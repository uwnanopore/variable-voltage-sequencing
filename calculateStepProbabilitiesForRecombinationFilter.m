function [P] = calculateStepProbabilitiesForRecombinationFilter(X, regularization)

% FUNCTIONALITY:
% function [P] = calculateStepProbabilitiesForRecombinationFilter(X, regularization)
%   Uses a system of svms and logits built in
%   buildRecombinationSVM.m in order to determine
%   step/back/hold/skip probabilities for normalized, calibrated levels.
%   For use within the recombination filter.
%
% INPUTS:
%   X: N_transitions x 6 array of the "transition features". Entries 1:3 in
%      each row are "level ii" normalized, calibrated features, entries 4:6
%      are "level ii + 1" normalized, calibrated features.
%   regularization: boolean saying whether (true) or not (false) to
%                   implement upper/lower bounds on the svm-based probabilities. 
%                   Defaults to true. 
%
% OUTPUTS:
%   P: 4 x N_transitions array of the step probabilities for each
%   transition. Entry 1 = p_step, 2 = p_back, 3 = p_hold, 4 = p_skip.
%

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
% abs_min = 1e-2;
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
    
    min_p_SH = abs_min;
    max_p_SH = abs_max;
    
elseif regularization
    min_p_SK = max(1 - validation.step_v_skip.corr_rate_K, abs_min);
    max_p_SK = min(validation.step_v_skip.corr_rate_S, abs_max);
    
    min_p_SB = max(1 - validation.step_v_back.corr_rate_B, abs_min);
    max_p_SB = min(validation.step_v_back.corr_rate_S, abs_max);
    
    min_p_SH = max(1 - validation.step_v_hold.corr_rate_H, abs_min);
    max_p_SH = min(validation.step_v_hold.corr_rate_S, abs_max);
    
end

reglogit = @(score, params, min_p, max_p) (logit(score, params) .* ((logit(score, params) > min_p) & (logit(score, params) < max_p))) + ...
    (max_p .* (logit(score, params) >= max_p)) + ...
    (min_p .* (logit(score, params) <= min_p));

% define (regularized) sigma function
sigma = @(score, params, min_p, max_p) reglogit(score, params, min_p, max_p) ./ (1 - reglogit(score, params, min_p, max_p));

% initialize output
P = nan(4, size(X, 1));

% loop over transitions
for cT = 1:size(X, 1)
    X_ii = X(cT, :);
    
    % score X_ii on all svms
    [~, score_SK] = predict(svm.step_v_skip.ClassificationSVM, X_ii); score_SK = score_SK(2);
    [~, score_SB] = predict(svm.step_v_back.ClassificationSVM, X_ii); score_SB = score_SB(2);
    [~, score_SH] = predict(svm.step_v_hold.ClassificationSVM, X_ii); score_SH = score_SH(2);

    % generate H matrix
    H = [ ...
        1, zeros(1, 0), -1 * sigma(score_SB, params.step_v_back, min_p_SB, max_p_SB), zeros(1, 2) ; ...
        1, zeros(1, 1), -1 * sigma(score_SH, params.step_v_hold, min_p_SH, max_p_SH), zeros(1, 1) ; ...
        1, zeros(1, 2), -1 * sigma(score_SK, params.step_v_skip, min_p_SK, max_p_SK), zeros(1, 0) ; ...
        ones(1, 4)];

    % H * [p_step p_back p_hold p_skip]' = [0 0 0 1]' so
    % H^-1 * [0 0 0 1]' = [p_step p_back p_hold p_bad]'
    p_i = inv(H) * [zeros(size(H, 1) - 1, 1) ; 1];
    
    P(:, cT) = p_i;
end

end


