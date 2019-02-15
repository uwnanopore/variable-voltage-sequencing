%% preliminaries
% structure of x_bad vectors:
%   entries 1-3: 3 pc values of previous level
%   entries 4-6: 3 pc values of current level
%   entries 7-9: 3 pc values of next level
%   entry 10   : max absolute value of the 101 (median-normalized) calibrated features for the level
%   entry 11   : mean square value of the residual from 3pc fit to 101 features
%   entry 12   : best log score relative to map

% load trainging data
examples = load('removal_svm_training_examples.mat');
examples = examples.examples;

% remove badly conditioned examples (i.e. containing NaNs, infs)
% goods
ill = zeros(1, size(examples.good, 2));
nan_count = sum(isnan(examples.good), 1);
ill(nan_count > 0) = 1;
inf_count = sum(abs(examples.good) > 1e20, 1);
ill(inf_count > 0) = 1;
examples.good(:, ill == 1) = [];

% bads
ill = zeros(1, size(examples.bad, 2));
nan_count = sum(isnan(examples.bad), 1);
ill(nan_count > 0) = 1;
inf_count = sum(abs(examples.bad) > 1e20, 1);
ill(inf_count > 0) = 1;
examples.bad(:, ill == 1) = [];

clear nan_count inf_count ill

%% divide up into training and validation sets
% how big to make T and V sets
training_size_bad = 800;
training_size_good = training_size_bad * 1;

validation_size_bad = 400;
validation_size_good = validation_size_bad * 1;

% get GOOD examples
ix = randperm(size(examples.good, 2));
training_ix = ix(1:training_size_good);
validation_ix = ix((training_size_good + 1):(training_size_good + validation_size_good));
training_examples.good = examples.good(:, training_ix);
validation_examples.good = examples.good(:, validation_ix);

% get BAD examples
ix = randperm(size(examples.bad, 2));
training_ix = ix(1:training_size_bad);
validation_ix = ix((training_size_bad + 1):(training_size_bad + validation_size_bad));
training_examples.bad = examples.bad(:, training_ix);
validation_examples.bad = examples.bad(:, validation_ix);

clear ix training_ix validation_ix *_size_*

%% train SVMs

% build table
X_train = [training_examples.good' ; training_examples.bad'];
y_train = [1 .* ones(size(training_examples.good, 2), 1) ; -1 .* ones(size(training_examples.bad, 2), 1)];
T = table(X_train, y_train);

% train svm
[svm, val] = trainRemovalClassifierQuadSVM(T);

clear X_train y_train T

%% logit fitting to convert svm scores to probabilities

% define function
logit = @(x, params) glmval(params, x, 'logit');

% predict labels and scores on validation set
[labels, scores] = predict(svm.ClassificationSVM, [validation_examples.good' ; validation_examples.bad']);
scores = scores(:, 2);
target = [true(size(validation_examples.good, 2), 1) ; false(size(validation_examples.bad, 2), 1)];

% maximum likelihood estimator
params = glmfit(scores, target, 'binomial');

%% save outputs
% save('svm_for_bad_level_filtering', 'svm');
% save('logit_params_for_bad_level_filtering', 'params');
% save('logit_for_badlevel_filtering', 'logit');
