function [x_3_r, k_3_r, x_101_r, k_101_r, tf_to_r, npts_r, did_backstep_r] = ...
    reorderFilter(x_3_tf, k_3_tf, x_101_tf, k_101_tf, npts_tf, did_backstep_tf, varargin)

% FUNCTIONALITY:
%   Uses SVMs to evaluate Pstep/Pback/Pskip at each transition, then
%   reorders the levels to give the best overall probability
% VARARGIN:
%   maxiterations: (5) most times can loop implementation of reorderer
%                  before forcing exit
%   map: the kmer model to calibrate to
%   prior: sets a prior on the P_step / P_back / P_skip

%% Defaults and varargin
% defaults
% prin comps
pc = load('principal_components.mat');
pc = pc.principal_components;
max_iterations = 5;
prior = [.5 .1 .05]';

% varargin
for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'maxiterations'
            max_iterations = varargin{cV + 1};
        case 'map'
            map = varargin{cV + 1};
        case 'prior'
            prior = varargin{cV + 1};
            if size(prior, 1) < size(prior, 2)
                prior = prior';
            end
    end
end
 
%% Calibrate features
% Generate Calibration Strand
% load prediction map
if ~exist('map', 'var')
    map = load('pore_model_6mer_variable_voltage.mat');
    map = map.model;
end

% generate random calibration sequence
if ~exist('refseq', 'var')
    alphabet = 'ACGT';
    refseq = alphabet(randi(4, 1, 5e5));
    clear alphabet
end

% generate prediction for calibration sequence
n_3_ref = levelPredPipeline(refseq, map, false);
n_101_ref = pc * n_3_ref;

% Calibrate to reference
%   1, normalize measured levels
n_101 = (pc * x_3_tf) - normalizerAC(pc * x_3_tf);

%   2, calculate calibration of measured to reference
meas_cal = reshape(n_101, 1, []);
ref_cal = reshape(n_101_ref, 1, []);
[scale, offset] = calibrate_by_iqr(meas_cal(~isnan(meas_cal)), ref_cal(~isnan(ref_cal)), false, [-inf inf]);

%   3, apply calibration to measured levels
n_101 = (scale .* n_101) + offset;

%   4, reduce calibrated data back to 3 pcs
n_3 = pc' * n_101;

%% Looped implementation of reorderer
exitflag = false;
iter = 0;
new_ix = 1:size(n_3, 2);

while ~exitflag
    % update iteration
    iter = iter + 1;
    
    % new_ix is now previous_ix
    previous_ix = new_ix;
    
    % calculate SVM scores
    logPs = calculateStepProbabilitiesForReorderFilter(n_3(:, previous_ix), true, 'prior', prior);
    
    % calculate optimal reordering
    level_placements = calculateReordering(logPs');
    level_placements = level_placements';
    
    % remove gaps from level placements
    for cL = max(level_placements):-1:1
        if ~any(level_placements == cL)
            level_placements(level_placements > cL) = level_placements(level_placements > cL) - 1;
        end
    end
    
    % determine new indexing overall
    new_ix = level_placements(previous_ix);
    
    % exit loop if we've either made no changes, or reached the max
    % iterations
    if (sum(abs(new_ix - previous_ix)) == 0) || (iter >= max_iterations)
        exitflag = true;
    end
end

level_placements = new_ix';

% store reordered information
x_101_r = x_101_tf(:, level_placements');
k_101_r = k_101_tf(level_placements');
npts_r = npts_tf(level_placements');

x_3_r = x_3_tf(:, level_placements');
k_3_r = k_3_tf(level_placements');

did_backstep_r = did_backstep_tf(level_placements');

tf_to_r = level_placements';

end


