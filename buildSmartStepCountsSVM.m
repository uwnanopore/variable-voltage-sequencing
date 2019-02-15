%%
% General Architecture
% Want probabilities of different step sizes from P_1 (step) to P_2 (single
% skip) up to P_12 (max skip, to any new kmer state)
% Stacked series of classifiers C_1 to C_11
% Each classifier C_i is characterized by a logit function L_i and an SVM
% score s_i. 
% Classifier C_i tells whether step was likely to be step of size i or step
% of size greater than i, given that was not a step of size less than i
% P(i | ~ <i) = l_i(s_i)
% ----------------------------
% P(1) = P(1)
% P(2) = P(2 | ~ <2) * P(~ <2)
% P(3) = P(3 | ~ <3) * P(~ <3)
% ...
% P(N) = P(N | ~ <N) * P(~ <N)
% ----------------------------
% P(1) = l_1(s_1)
% P(2) = l_2(s_2) * [1 - l_1(s_1)]
% P(3) = l_3(s_3) * [1 - l_2(s_2)] * [1 - l_1(s_1)]
% ...
% P(N) = l_N(s_N) * [1 - l_N-1(s_N-1)] * ... * [1 - l_1(s_1)]
% -----------------------------------------------------------
%
%

%% load examples
examples = load('smart_step_counts_training_examples.mat');
examples = examples.examples;

%% CLASSIFIERS
% determine some basic parameters for classifier training
% find minumum number examples extant of any step size
min_examples = min((cellfun(@numel, examples) / 6));
% fix validation size
validation_size = 1000;
% training size is min number of leftovers after forcing validation size
training_size = min_examples - validation_size;

% set rate of larger skips
extension_rate = 1/5;

% build ladder of classifiers
% initialize classifiers
classifiers = cell(1, 11);
for step_size = 1:11
    fprintf('Building Classifier %d/%d\n', step_size, 11)
    % initialize
    classifiers{step_size} = struct;
    % build example set for the step size in question
    % randomly permute the examples
    ix = randperm(size(examples{step_size}, 1));
    % grab training_size worth of random ix
    t_ix = ix(1:training_size);
    % grab validation_size worth of random ix
    v_ix = ix((training_size + 1):(training_size + validation_size));
    % build positive training examples
    t_examples_pos = examples{step_size}(t_ix, :);
    % build positive validation examples
    v_examples_pos = examples{step_size}(v_ix, :);
    
    % build example set for the step sizes greater than the step size in
    % question
    % we have 12 - step_size possible larger steps
    neg_example_multiplicity = ones(1, 12 - step_size);
    % each subsequent example multiplicity m_i should be extension_rate *
    % m_(i-1)
    exponent_values = 0:1:(size(neg_example_multiplicity, 2) - 1);
    neg_example_multiplicity = neg_example_multiplicity .* (extension_rate .^ exponent_values);
    % normalize neg_example_multiplicity so the sum of all example sizes is
    % training_size or validation_size
    t_neg_example_multiplicity = (training_size / sum(neg_example_multiplicity)) * neg_example_multiplicity;
    v_neg_example_multiplicity = (validation_size / sum(neg_example_multiplicity)) * neg_example_multiplicity;
    % ceil off all multiplicity values
    t_neg_example_multiplicity = ceil(t_neg_example_multiplicity);
    t_neg_example_multiplicity(t_neg_example_multiplicity < 1) = 1;
    v_neg_example_multiplicity = ceil(v_neg_example_multiplicity);
    v_neg_example_multiplicity(v_neg_example_multiplicity < 1) = 1;
    % reduce 1st entry by enough to restore correct total to training_size
    % or validation_size
    t_neg_example_multiplicity(1) = t_neg_example_multiplicity(1) - (sum(t_neg_example_multiplicity) - training_size);
    v_neg_example_multiplicity(1) = v_neg_example_multiplicity(1) - (sum(v_neg_example_multiplicity) - validation_size);
    % fill in the t_examples_neg and v_examples_neg with the correct number
    % of examples from each large step size
    % initialize
    t_examples_neg = zeros(0, 6);
    v_examples_neg = zeros(0, 6);
    % loop over all extended step sizes
    for cS = 1:size(neg_example_multiplicity, 2)
        % randomly permute the entries in this set of examples
        ix = randperm(size(examples{step_size + cS}, 1));
        % select out t_neg_example_multiplicity(cS) and v_neg_example_multiplicity(cS)
        t_ix = ix(1:t_neg_example_multiplicity(cS));
        v_ix = ix((t_neg_example_multiplicity(cS) + 1):(t_neg_example_multiplicity(cS) + v_neg_example_multiplicity(cS)));
        t_examples_neg = [t_examples_neg ; examples{step_size + cS}(t_ix, :)];
        v_examples_neg = [v_examples_neg ; examples{step_size + cS}(v_ix, :)];
    end
    
    % train classifier (quadratic polynomial kernel SVM) to distinguish
    % between step_size and >step_size
    X_train = [t_examples_pos ; t_examples_neg];
    y_train = [ones(size(t_examples_pos, 1), 1) ; -1 .* ones(size(t_examples_neg, 1), 1)];
    T = table(X_train, y_train);
    [classifiers{step_size}.svm, classifiers{step_size}.svm_validation] = trainClassifierQuadSVM(T);
    
    % logit regression on validation set to convert svm scores to
    % probabilities
    logit = @(x, params) glmval(params, x, 'logit');
    % use svm to predict labels and scores on the validation set
    [labels, scores] = predict(classifiers{step_size}.svm.ClassificationSVM, [v_examples_pos ; v_examples_neg]);
    % determination is second column of scores
    out = scores(:, 2);
    % target is trues for _pos and falses for _neg
    target = [true(size(v_examples_pos, 1), 1) ; false(size(v_examples_neg, 1), 1)];
    % logit parameters from global likelihood maximization
    classifiers{step_size}.params = glmfit(out, target, 'binomial');
    
    % evaluate accuracy on the validation set
    classifiers{step_size}.corr_rate_pos = sum(out(1:size(v_examples_pos, 1)) >= 0) / size(v_examples_pos, 1);
    classifiers{step_size}.corr_rate_neg = sum(out((size(v_examples_pos, 1) + 1):(size(v_examples_pos, 1) + size(v_examples_neg, 1))) <= 0) / size(v_examples_neg, 1);    
    
end

%% save?
save_on = false;
if save_on
    save('classifiers_for_smart_step_counts.mat', 'classifiers');
end

