function [keep_ix, uf_to_f, threshold] = removalFilter(x_101, varargin)

% FUNCTIONALITY:
%   Implements the bad level classification SVM to iteratively and remove
%   bad levels from an un-normalized, un-calibrated, un-filtered event.
% 
% INPUTS:
%   x_101: 101 x N_levels array of un-normalized, un-calibrated 101 dimensional feature
%          vectors for event to be filtered
%
% VARARGIN:
%   threshold: sets probability threshold for level to be called as
%              "bad". default is set to the logit(0), corresponding to the
%              SVM boundary. user input can range from 0 (must have P(good)
%              = 0 to be called "bad") to 1 (must have P(good) <= 1 to be
%              called "bad"). Closer to 1, the more levels will be called
%              "bad". default corresponds to value of ~0.45.
%
% OUTPUTS:
%   keep_ix: 1 x N_levels array of 1 (good levels to keep) and 0 (bad
%            levels to throw)
%   uf_to_f: 1 x N_levels_uf array of alignment of uf levels to new f
%            levels
%   threshold: returns P(good) required to be below to call level "bad"

%% Handle varargin
for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'threshold'
            user_input_threshold = varargin{cV + 1};
            if user_input_threshold >= 1
                fprintf('WARNING: user threshold is above 1, reverting to default.\n')
                clear user_input_threshold
            elseif user_input_threshold <= 0
                fprintf('WARNING: user threshold is below 0, reverting to default.\n')
                clear user_input_threshold
            end
    end
end

%% Generate calibration strand
% load prediction map
if ~exist('map', 'var')
    map = load('pore_model_6mer_variable_voltage.mat');
    map = map.model;
end

% load principal components
if ~exist('pc', 'var')
    pc = load('principal_components.mat');
    pc = pc.principal_components;
end

% generate random calibration sequence
if ~exist('refseq', 'var')
    alphabet = 'ACGT';
    refseq = alphabet(randi(4, 1, 5e5));
    clear alphabet
end

% generate prediction for calibration sequence
if ~exist('x_ref_3', 'var')
    [x_ref_3] = levelPredPipeline(refseq, map, false);
    x_ref_101 = pc * x_ref_3;
end

%% Filtering loop
% initialization for loop
exitflag = false;
x_current_101 = x_101;
scale = 1;
offset = 0;
ix_orig = 1:size(x_current_101, 2);
keep_ix = ones(1, size(x_current_101, 2));
iter = 0;

% parameters
min_levels_remaining = 10; % if less than this number of levels remaining,
end_trimming = 2; % buffer at start and end to prevent repetitive whittling; will be trimmed at end of process
foo = @(x, mu, cov, stiff) (1 / sqrt(((2 * pi) ^ size(x, 1)) * det(cov))) * exp((-1 / 2) * (x - mu)' * stiff * (x - mu)); % function calculating multivariate normal distribution level match score
max_iter = 10; % places a cap on how many iterations to go through before forcing completion

% load classifier materials
svm = load('svm_for_bad_level_filtering.mat');
svm = svm.svm;
params = load('logit_params_for_bad_level_filtering.mat');
params = params.params;
logit = load('logit_for_badlevel_filter.mat');
logit = logit.logit;

% define "good" threshold
% default: "good" are scores 0 and above (after passed through logit)
if exist('user_input_threshold', 'var')
    threshold = user_input_threshold;
else
    threshold = logit(0, params);
end

% do the loop
while ~exitflag
    
    iter = iter + 1;
    % make sure we have a minimum number of levels
    if size(x_current_101, 2) <= min_levels_remaining
        exitflag = true;
        continue
    end
    
    % determine trimmed indices
    trimmed_ix = (1 + end_trimming) : (size(x_current_101, 2) - end_trimming);
    
    % calibrate to reference
    % first, normalize measured levels by MEDIAN only
    x_current_101 = x_current_101 - repmat(nanmedian(x_current_101(:, trimmed_ix), 2), 1, size(x_current_101, 2));
    
    % second, calculate calibration of measured to reference
    meas_cal = reshape(x_current_101, 1, []);
    ref_cal = reshape(x_ref_101, 1, []);
    [scale_current, offset_current] = calibrate_by_iqr(meas_cal(~isnan(meas_cal)), ref_cal(~isnan(ref_cal)), false, [-inf inf]);
    
    % apply calibration to the measured levels
    x_current_101 = (scale_current .* x_current_101) + offset_current;
    
    % keep track of cummulative scale/offset
    scale = scale_current * scale;
    offset = (scale_current * offset) + offset_current;
    
    % reduce measurement to 3 pcs
    x_current_3 = pc' * x_current_101;
    
    % initialize x_bad
    x_bad = nan(12, size(x_current_3, 2));
    
    % fill x_bad
    % I, entries 1:3 are the N-1 level pcs
    x_bad(1:3, 2:end) = x_current_3(:, 1:(end - 1));
    % II, entries 4:6 are the N level pcs
    x_bad(4:6, 1:end) = x_current_3(:, 1:end);
    % III, entries 7:9 are the N+1 level pcs
    x_bad(7:9, 1:(end - 1)) = x_current_3(:, 2:end);
    % IV, entry 10 is the extreme value of x_101
    x_bad(10, :) = max(abs(x_current_101));
    % V, entry 11 is the residual variance relative to the 3 pc reduction
    residual_3pc = x_current_101 - (pc * x_current_3);
    var_residual_3pc = mean(residual_3pc .^ 2, 1);
    x_bad(11, :) = var_residual_3pc;
    % VI, entry 12 is the best score relative to map
    % loop over all levels
    % SLOW: only calculate on first run-through, then just track values
    % through
    if iter == 1
        map_matches = nan(1, size(x_current_3, 2));
        for cL = 1:size(x_current_3, 2)
            x = x_current_3(:, cL);
            scores = zeros(1, size(map.mean_3, 2));
            % loop over all map entries
            for cK = 1:size(map.mean_3, 2)
                % extract map values
                mu = map.mean_3(:, cK);
                cov = inv(map.stiffness_3{cK});
                stiff = map.stiffness_3{cK};
                % calculate score
                scores(cK) = log(foo(x, mu, cov, stiff));
            end
            % if mod(cL, 200) == 0
            %         fprintf('Done %d/%d\n', cL, size(x_current_3, 2))
            % end
            best_score = max(scores);
            x_bad(12, cL) = best_score;
            map_matches(cL) = best_score;
        end
    else
        x_bad(12, :) = map_matches_current;
    end
    
    % use classifier to determine label and score
    [label, score] = predict(svm.ClassificationSVM, x_bad');
    score = score(:, 2); % score > 0 means good, score < 0 means bad (score tells how far into "good" region example lies
    
    % use logit to convert scores to probabilities
    p_score = logit(score, params);
    
    % determine indices of good/bad levels
    is_good = ones(1, size(x_current_101, 2));
    is_good(p_score < threshold) = 0;
    
    % for edge consideration, keep first Ntrim and last Ntrim levels
    is_good(1:end_trimming) = 1;
    is_good((end - end_trimming + 1):end) = 1;
    
    % update kept levels and tracking indices
    keep_ix(ix_orig(is_good == 0)) = 0;
    ix_orig = ix_orig(is_good == 1);
    
    % update the current set of levels to look at
    x_current_101 = x_101(:, keep_ix == 1);
    map_matches_current = map_matches(keep_ix == 1);
    
    % check to see if we've removed anything new; if not, exit loop
    if sum(is_good) == length(is_good)
        exitflag = true;
    end
    
    % check to see if we are going to exceed the max_iter; if so, exit loop
    if iter >= max_iter
        exitflag = true;
    end
end

% chop off end trimming levels
keep_ix(1:end_trimming) = 0;
keep_ix((end - end_trimming + 1):end) = 0;

% create uf_to_f alignment
uf_to_f = cumsum(keep_ix);
uf_to_f(keep_ix == 0) = NaN;

end
