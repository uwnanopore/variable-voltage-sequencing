function [p, p_bad] = smartStepCounts(x_101, x_3, calibrated, calibration, varargin)

% FUNCTIONALITY:
%
% INPUTS:
%   x_101: (101 x n_levels) features array (should be unnormalized and uncalibrated)
%   x_3: (3 x n_levels) features array (should be unnormalized and uncalibrated)
%   calibrated: (boolean) tells whether passed features have yet been
%               calibrated
%   calibration: (1 x 2) array of [scale, offset] to use if passed
%                "calibrated" is true
%
% OUTPUTS:
%   p: (n_levels x max_sequential_bad x 12) array giving transition
%      probabilities into each level. 2nd dimension gives gap sizes, 3rd
%      dimension gives step sizes
%   p_bad: (n_levels x 1) array of bad probability of each level
%
% VARARGIN:
%   dumb_mode: turns on dumb_mode where we just use default step
%              probabilities
%   defaultpbad: allows user to provide the default p_bad
%   defaultpskip: allows user to provide the default p_skip
%   defaultpextension: allows user to provide the default p_extend
%   maxsequentialbad: provide the max consecutive bad levels allowed to be
%                     called
%   regularization: turn on (true) or off (false) the logit regularization
%                   by verification calling rates
%   svmtrust: how many step sizes out to trust the SVM classifiers to
%             provide useful information
%   pfalloff: fall off probability rate for extending steps after we no
%             longer trust the svms (0 = max fall off, 1 = no fall off)
%   priortostep: prior probability on how often steps should occur
%                relative to skips. Default is 0.5, in which case
%                we trust exactly what comes out of the svm ladder.
%                A value above 0.5 biases things for more steps
%   map: allows user to pass in a map to save on loading time

% set defaults
dumb_mode = false;
default_p_bad = 0.05;
default_p_skip = 0.125;
default_p_extension = 0.08;
max_sequential_bad = 3;
regularization = true;
svm_trust = 2;
p_fall_off = .075;
prior_to_step = .75;

% handle varargin
for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'dumbmode'
            dumb_mode = true;
        case 'defaultpbad'
            default_p_bad = varargin{cV + 1};
        case 'defaultpskip'
            default_p_skip = varargin{cV + 1};
            default_p_extension = default_p_skip;
        case 'defaultpextension'
            default_p_extension = varargin{cV + 1};
        case 'maxsequentialbad'
            max_sequential_bad = varargin{cV + 1};
        case 'regularization'
            regularization = varargin{cV + 1};
        case 'svmtrust'
            svm_trust = varargin{cV + 1};
        case 'pfalloff'
            p_fall_off = varargin{cV + 1};
        case 'priortostep'
            prior_to_step = varargin{cV + 1};
            % regularize the step prior from getting too extreme
            if (prior_to_step > 0.99)
                prior_to_step = 0.99;
            elseif prior_to_step < 0.01
                prior_to_step = 0.01;
            end
        case 'map'
            map = varargin{cV + 1};
    end
end

% initialize p, p_bad
% p is an (N_levels) x (max_sequential_bad) x (12) array
% the third dimension corresponds to step sizes ranging from 1 to 12
p = nan(size(x_3, 2), max_sequential_bad + 1, 12);
% p_bad is an (N_levels) x (1) array
p_bad = nan(size(x_3, 2), 1);

% if dumb mode, we just fill in with uniform step counts
if dumb_mode
    % fill p_bad with default
    p_bad(:, :) = default_p_bad;
    % fill 1st sheet with step default
    p(:, :, 1) = 1 - default_p_skip;
    % fill subsequent sheets with skip and skip extension defaults
    for cS = 2:12
        p(:, :, cS) = default_p_skip * (default_p_extension ^ (cS - 2));
    end
    % enforce normalization of p
    p = p ./ repmat(sum(p, 3), 1, 1, 12);
    return
end

% if not dumb mode...
% load necessary information for calibration
if ~exist('map', 'var')
    map = load('pore_model_6mer_variable_voltage.mat');
    map = map.model;
end
pc = load('principal_components.mat');
pc = pc.principal_components;

% if not passed in calibrated, we need to calibrate to a random sequence
% generated from the map
if ~calibrated
    
    % generate random calibration sequence
    alphabet = 'ACGT';
    random_sequence = alphabet(randi(4, 1, 5e5));
    [random_x_3] = levelPredPipeline(random_sequence, map, false);
    random_x_101 = pc * random_x_3;
    
    % get normalized features for calibration
    n_3 = pc' * ((pc * x_3) - normalizerAC(pc * x_3));
    n_101 = pc * n_3;
    
    % calculate calibration to the random sequence
    [calibration(1), calibration(2)] = calibrate_by_iqr(reshape(n_101, 1, []), reshape(random_x_101, 1, []), false, [-inf inf]);
end

% generate x_bad
x_bad = makeXBad(x_101, calibration, map, pc);

% generate x_skips for each of the possible spanning bads
x_skips = cell(1, max_sequential_bad + 1);
for cS = 0:max_sequential_bad
    x_skips{cS + 1} = makeXSkip(x_3, calibration, pc, cS);
end

% define extreme bounds on regularization to be implemented at all times
abs_min = 0.05;
abs_max = 1 - abs_min;

% load basic logit function
logit = load('logit_for_backstep_filter.mat'); 
logit = logit.logit;

% define regularized logit function
reglogit = @(score, params, min_p, max_p) (logit(score, params) .* ((logit(score, params) > min_p) & (logit(score, params) < max_p))) + ...
    (max_p .* (logit(score, params) >= max_p)) + ...
    (min_p .* (logit(score, params) <= min_p));

% Good vs Bad Classification
% load classifier materials
svm = load('svm_for_bad_level_filtering.mat');
svm_bad = svm.svm;
params = load('logit_params_for_bad_level_filtering.mat'); 
params_bad = params.params;
validation = load('validation_for_bad_level_filtering.mat');
val_bad = validation.validation;

% define regularization parameters for G/B classifier
if ~regularization
    min_p_GB = abs_min;
    max_p_GB = abs_max;
    
elseif regularization
    min_p_GB = max(1 - val_bad.good_v_bad.corr_rate_B, abs_min);
    max_p_GB = min(val_bad.good_v_bad.corr_rate_G, abs_max);
end

% make sure that we are not doing worse than coin flip
if (min_p_GB >= .5) || (max_p_GB <= .5)
    min_p_GB = .5 - 1e10;
    max_p_GB = .5 + 1e10;
end

% predict using classifier
[~, score] = predict(svm_bad.ClassificationSVM, x_bad');
score = score(:, 2);

% use logit to convert score to probability
p_score = reglogit(score, params_bad, min_p_GB, max_p_GB);

% fill in p_bad based on these results
p_bad = 1 - p_score;

% update the first and last entries in p_bad, as these are currently NaN
p_bad([1 length(p_bad)]) = default_p_bad;

% Step Size Classification
% load classifier materials
temp = load('classifiers_for_smart_step_counts.mat'); 
classifiers = temp.classifiers; clear temp

% initialize storage of marginal probabilities:
% p_marginal{gap_size, ii} = logit_ii(svm_score_ii)
% size is (gap_sizes) x (n_step_sizes - 1)
p_marginal = cell(size(p, 2), size(p, 3) - 1);

% loop over step sizes
for step_size = 1:size(p_marginal, 2)
    % define regularization for each step size classifier
    if ~regularization
        min_p_ss = abs_min;
        max_p_ss = abs_max;
    elseif regularization
        min_p_ss = max(1 - classifiers{step_size}.corr_rate_neg, abs_min);
        max_p_ss = min(classifiers{step_size}.corr_rate_pos, abs_max);
    end
    
    % make sure that we are not doing worse than coin flip
    if (min_p_ss >= .5) || (max_p_ss <= .5)
        min_p_ss = .5 - 1e-3;
        max_p_ss = .5 + 1e-3;
    end
    
    % loop over all spanning distances
    for gap_size = 1:size(x_skips, 2)
        x_skip = x_skips{gap_size};
        
        % predict using classifier
        [~, score] = predict(classifiers{step_size}.svm.ClassificationSVM, x_skip);
        score = score(:, 2);
        
        % use logit to convert score to probability
        p_score = reglogit(score, classifiers{step_size}.params, min_p_ss, max_p_ss);
        
        % fill in the marginal probability we have calculated
        p_marginal{gap_size, step_size} = p_score;
        
    end
end

% calculate the complementary probabilities at each step size
% q_marginal{gap_size, ii} = 1 - logit_ii(svm_score_ii)
q_marginal = cell(size(p, 2), size(p, 3) - 1);

% loop over step sizes
for step_size = 1:size(q_marginal, 2)
    % loop over spanning distances
    for gap_size = 1:size(q_marginal, 1)
        % take 1 - the marginal p
        q_marginal{gap_size, step_size} = 1 - p_marginal{gap_size, step_size};
    end
end

% we want to calculate the cummulate complementary probabilities
% q_marginal_cum{gap_size, ii} = Prod(from 1 to ii) of (1 - logit_ii(svm_score_ii))
q_marginal_cum = cell(size(p, 2), size(p, 3) - 1);
% loop over step sizes
for step_size = 1:size(q_marginal_cum, 2)
    % loop over spanning distances
    for gap_size = 1:size(q_marginal_cum, 1)
        % initialize each entry with ones
        q_marginal_cum{gap_size, step_size} = ones(size(q_marginal{gap_size, step_size}));
        % loop over all step sizes up to step_size
        for ii = 1:step_size
            q_marginal_cum{gap_size, step_size} = q_marginal_cum{gap_size, step_size} .* q_marginal{gap_size, ii};
        end
    end
end

% figure out the total probability at each step now
% p_tot{gap_size, ii} = L_ii(svm_score_ii) * prod(kk = 1 to ii - 1) [1 - L_kk(svm_score_kk)]
% p_tot{gap_size, ii} = p_marginal{gap_size, ii} .* q_marginal_cum{gap_size, ii - 1}
% p_tot should be (gap_sizes) x (n_step_sizes)
p_tot = cell(size(p, 2), size(p, 3));
% loop over step sizes
for step_size = 1:size(p_tot, 2)
    % loop over spanning distances
    for gap_size = 1:size(p_tot, 1)
        if step_size == 1
            % if we are in step size == 1, we don't need any q_marginals
            p_tot{gap_size, step_size} = p_marginal{gap_size, step_size};
        elseif (step_size > 1) && (step_size < size(p, 3))
            % if we are in step sizes from 1 to max - 1, we need both p and
            % q marginals
            p_tot{gap_size, step_size} = p_marginal{gap_size, step_size} .* q_marginal_cum{gap_size, step_size - 1};
        elseif step_size == size(p, 3)
            % if we are in terminal step size, we only want the q marginal
            p_tot{gap_size, step_size} = q_marginal_cum{gap_size, step_size - 1};
        end
    end
end

% set all nan entries in p_tot to the correct default values
% loop over step sizes
for step_size = 1:size(p_tot, 2)
    % loop over gap sizes
    for gap_size = 1:size(p_tot, 1)
        if step_size == 1
            % if single step, use default p_step
            p_tot{gap_size, step_size}(isnan(p_tot{gap_size, step_size})) = 1 - default_p_skip;
        elseif step_size > 1
            % if a sort of skip, use skip and extension penalties
            p_tot{gap_size, step_size}(isnan(p_tot{gap_size, step_size})) = default_p_skip * (default_p_extension ^ (step_size - 2));
        end
    end
end

% if we are out at a step size beyond where we trust the SVM
% classifiers, revert the marginal probability to a default value
% of p_fall_off
% loop over step sizes
for step_size = 1:size(p_tot, 2)
    % loop over gap sizes
    for gap_size = 1:size(p_tot, 1)
        % if we are at a step size larger than the trusted svms cover, we
        % set the probability of this step ii to be p_ii = p_fall_off * p_(ii - 1)
        if step_size > svm_trust
            p_tot{gap_size, step_size} = p_tot{gap_size, step_size - 1} .* p_fall_off;
        end
    end
end

% use p_tot to fill in the final p matrix
% loop over step sizes
for step_size = 1:size(p, 3)
    % loop over gap sizes
    for gap_size = 1:size(p, 2)
        % fill in all levels for this gap and step size
        p(:, gap_size, step_size) = p_tot{gap_size, step_size};
    end
end

% apply the prior to step
% factor is prior_prob / (1 - prior_prob)
prior_step_factor = prior_to_step / (1 - prior_to_step);
% we multiply this factor into the step probabilities, then renormalize
% later
p(:, :, 1) = prior_step_factor .* p(:, :, 1);

% make sure each level of p is normalized so sum of all step size
% probabilities leads to 1
% loop over gap sizes
for gap_size = 1:size(p, 2)
    % calculate total probability in all step sizes
    total_probability = sum(p(:, gap_size, :), 3);
    % divide all step probabilities by the total probability to normalize
    % the sum to 1
    p(:, gap_size, :) = p(:, gap_size, :) ./ repmat(total_probability, 1, 1, 12);
end

end

% ----------------------------------------------------------------------------------- %
%% function to calculate x_skip
function x_skip = makeXSkip(x_3, calibration, pc, spanning_length)

% normalize and calibrate x_3
n_101 = (pc * x_3) - normalizerAC(pc * x_3);
n_101 = (calibration(1) .* n_101) + calibration(2);
n_3 = pc' * n_101;

% initialize x_skip
x_skip = nan(size(n_3, 2), 6);

% fill in the "to" state
x_skip(:, 4:6) = n_3';

% fill in the "from" state
x_skip((2 + spanning_length):end, 1:3) = n_3(:, 1:(end - spanning_length - 1))';

end

%% function to calculate x_bad
function x_bad = makeXBad(x_101, calibration, map, pc)

% function calculating multivariate normal distribution level match score
foo = @(x, mu, cov, stiff) (1 / sqrt(((2 * pi) ^ size(x, 1)) * det(cov))) * exp((-1 / 2) * (x - mu)' * stiff * (x - mu));

% simple normalize, and calibrate
x_101 = x_101 - repmat(nanmedian(x_101, 2), 1, size(x_101, 2));
x_101 = (calibration(1) .* x_101) + calibration(2);

% reduce 101 to 3
x_3 = pc' * x_101;

% initialize x_bad
x_bad = nan(12, size(x_3, 2));

% fill x_bad
% I, entries 1:3 are the ii-1 level pcs
x_bad(1:3, 2:end) = x_3(:, 1:(end - 1));
% II, entries 4:6 are the ii level pcs
x_bad(4:6, 1:end) = x_3(:, 1:end);
% III, entries 7:9 are the ii+1 level pcs
x_bad(7:9, 1:(end - 1)) = x_3(:, 2:end);
% IV, entry 10 is the extreme value of x_101
x_bad(10, :) = max(abs(x_101));
% V, entry 11 is the residual variance relative to the 3 pc reduction
residual_3pc = x_101 - (pc * x_3);
var_residual_3pc = mean(residual_3pc .^ 2, 1);
x_bad(11, :) = var_residual_3pc;
% VI, entry 12 is the best score relative to map
for cL = 1:size(x_3, 2)
    x = x_3(:, cL);
    scores = zeros(1, size(map.mean_3, 2));
    for cK = 1:size(map.mean_3, 2)
        mu = map.mean_3(:, cK);
        cov = inv(map.stiffness_3{cK});
        stiff = map.stiffness_3{cK};
        scores(cK) = log(foo(x, mu, cov, stiff));
    end
    best_score = max(scores);
    x_bad(12, cL) = best_score;
end

end
