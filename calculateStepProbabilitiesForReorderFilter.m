function [logP] = calculateStepProbabilitiesForReorderFilter(X, regularization, varargin)

% FUNCTIONALITY:
%   Uses a system of svms and logits built in
%   buildRecombinationSVM.m in order to determine
%   step/back/skip probabilities for normalized, calibrated levels.
%   For use within the reorder filter.
%
% INPUTS:
%   X: (3 x N_levels) array of normalized PC 3D feature vectors. Features
%      should be calibrated to the current version of the 6mer map
%
% OUTPUTS:
%   logP: (N_levels - 1 x 3) array of log probabilities for different step
%         types for each transition between levels in X. 1st column gives
%         logP(Step), 2nd column gives logP(sKip), 3rd column gives
%         logP(Back)

% build N_transitions x 6 array out of X
prior = [1 1 1]'; % prior on S/K/B
for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'prior'
            prior = varargin{cV + 1};
            if size(prior, 1) < size(prior, 2)
                prior = prior';
            end
    end
end

X = [X(:, 1:(end - 1))' X(:, 2:end)'];

% load svms and logits
svm = load('svm_for_recombination_filter.mat');
svm = svm.svm;
logit = load('logit_for_backstep_filter.mat');
logit = logit.logit;
params = load('logit_params_for_recombination_filtering.mat');
params = params.params;
validation = load('validation_for_recombination_filter.mat');
validation = validation.validation;

% define extreme bounds on regularization level, to be implemented at all
% times
abs_min = 1e-2;
abs_max = 1 - abs_min;

% default regularization condition
if ~exist('regularization', 'var')
    regularization = true;
end

% define regularized logit function
if ~regularization
    min_p_SK = abs_min;
    max_p_SK = abs_max;
    
    min_p_SB = abs_min;
    max_p_SB = abs_max;
    
elseif regularization
    min_p_SK = max(1 - validation.step_v_skip.corr_rate_K, abs_min);
    max_p_SK = min(validation.step_v_skip.corr_rate_S, abs_max);
    
    min_p_SB = max(1 - validation.step_v_back.corr_rate_B, abs_min);
    max_p_SB = min(validation.step_v_back.corr_rate_S, abs_max);
end

reglogit = @(score, params, min_p, max_p) (logit(score, params) .* ((logit(score, params) > min_p) & (logit(score, params) < max_p))) + ...
    (max_p .* (logit(score, params) >= max_p)) + ...
    (min_p .* (logit(score, params) <= min_p));

% define regularized sigma function
sigma = @(score, params, min_p, max_p) reglogit(score, params, min_p, max_p) ./ (1 - reglogit(score, params, min_p, max_p));

% initialize output
logP = nan(3, size(X, 1));

% loop over transitions
for cT = 1:size(X, 1)
    X_ii = X(cT, :);
    
    % score X_ii on all svms
    [~, score_SK] = predict(svm.step_v_skip.ClassificationSVM, X_ii); score_SK = score_SK(2);
    [~, score_SB] = predict(svm.step_v_back.ClassificationSVM, X_ii); score_SB = score_SB(2);
    
    % generate H matrix
    H = [ ...
        1, zeros(1, 0), -1 * sigma(score_SK, params.step_v_skip, min_p_SK, max_p_SK), zeros(1, 1) ; ...
        1, zeros(1, 1), -1 * sigma(score_SB, params.step_v_back, min_p_SB, max_p_SB), zeros(1, 0) ; ...
        ones(1, 3)];
    
    % H * [p_step p_back p_skip]' = [0 0 1]', so
    p_i = inv(H) * [zeros(size(H, 1) - 1, 1) ; 1];
    
    % apply prior
    p_i = p_i .* prior;
    p_i = p_i ./ (sum(p_i));    
    
    % store the log prob
    logP(:, cT) = log(p_i);
end

end