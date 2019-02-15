function [trained_classifier, validation_accuracy] = trainClassifierQuadSVM(training_data)

%% TRAINING

% extract predictors (X) and responses (y)
input_table = training_data;

% split matrices in input table into vectors
input_table.X_train_1 = input_table.X_train(:, 1);
input_table.X_train_2 = input_table.X_train(:, 2);
input_table.X_train_3 = input_table.X_train(:, 3);
input_table.X_train_4 = input_table.X_train(:, 4);
input_table.X_train_5 = input_table.X_train(:, 5);
input_table.X_train_6 = input_table.X_train(:, 6);

predictor_names = {'X_train_1', 'X_train_2', 'X_train_3', 'X_train_4', 'X_train_5', 'X_train_6'};
predictors = input_table(:, predictor_names);
response = input_table.y_train;
is_categorical_predictor = [false, false, false, false, false, false];

% train a classifier
svm_classifier = fitcsvm( ...
    predictors, ...
    response, ...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 2, ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', [-1; 1]);

% create the result struct with predict function
splitMatricesInTableFcn = @(t) [t(:, setdiff(t.Properties.VariableNames, {'X_train'})), array2table(table2array(t(:, {'X_train'})), 'VariableNames', {'X_train_1', 'X_train_2', 'X_train_3', 'X_train_4', 'X_train_5', 'X_train_6'})];
extractPredictorsFromTableFcn = @(t) t(:, predictor_names);
predictorExtractionFcn = @(x) extractPredictorsFromTableFcn(splitMatricesInTableFcn(x));
svmPredictFcn = @(x) predict(svm_classifier, x);
trained_classifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% add additional fields to the result struct
trained_classifier.RequiredVariables = {'X_train'};
trained_classifier.ClassificationSVM = svm_classifier;
trained_classifier.About = 'This struct is a trained model, based on training in classification learner R2017a.';
trained_classifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats( e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

%% VALIDATION

% perform cross validation
partitionedModel = crossval(trained_classifier.ClassificationSVM, 'KFold', 5);

% compute validation predictions
[validation_predictions, validation_scores] = kfoldPredict(partitionedModel);

% compute the validation accuracy
validation_accuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

end