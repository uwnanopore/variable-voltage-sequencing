%% load pcs
pc = load('principal_components.mat');
pc = pc.principal_components;

%% load training data
training_data = load('training_data_for_recombination_filter.mat');
training_data = training_data.training_data;

%% extract from training data
x_101 = training_data.x_101;
x_3 = training_data.x_3;
w = training_data.weight;

%% calculate level averages
mu_101 = zeros(101, size(x_101, 2));
mu_3 = zeros(3, size(x_101, 2));
cov_3 = cell(1, size(x_101, 2));
bad_cond = zeros(1, size(cov_3, 2));

for cL = 1:size(x_101, 2)
    if size(x_101{cL}, 2) == 0
        mu_101(:, cL) = zeros(101, 1);
        mu_3(:, cL) = zeros(3, 1);
        cov_3{cL} = zeros(3, 3);
        bad_cond(cL) = 1;
        continue
    end
    mu_101(:, cL) = sum(x_101{cL} .* repmat(w{cL}, 101, 1), 2) / sum(w{cL});
    mu_3(:, cL) = sum(x_3{cL} .* repmat(w{cL}, 3, 1), 2) / sum(w{cL});
    cov_3{cL} = weightedcov(x_3{cL}', w{cL}' ./ sum(w{cL}));
    if rcond(cov_3{cL}) < 1e-10
        bad_cond(cL) = 1;
    end
end

%% fill in ill-conditioned covariances with 90th percentile covariance
good_C = cov_3(bad_cond == 0);
dets = zeros(1, size(good_C, 2));
for cC = 1:size(good_C, 2)
    dets(cC) = det(good_C{cC});
end
dets90 = prctile(dets, 90);
[~, dets90loc] = min(abs(dets - dets90));
C90 = good_C{dets90loc};
cov_3(bad_cond == 1) = {C90};

%% decide which levels to include in training set
checkl = load('check_levels_for_recombination_training.mat'); 
checkl = checkl.checkl;
usel = load('use_levels_for_recombination_training.mat');
usel = usel.usel;

%% calibrate consensus to current 6mer map
map = load('pore_model_6mer_variable_voltage.mat');
map = map.model;
phixseq = fastaread('phix.fasta');
phixseq = phixseq.Sequence;
[x_pred, ~] = levelPredPipeline(phixseq, map, false);
x_pred = x_pred(:, region_ix);

cal_meas = pc * mu_3;
cal_meas(:, usel == 0) = [];
cal_meas = reshape(cal_meas, 1, []);

cal_ref = pc * x_pred;
cal_ref(:, usel == 0) = [];
cal_ref = reshape(cal_ref, 1, []);

[xcal, ycal] = fitprep(cal_meas, cal_ref);
foo = fit(xcal, ycal, 'poly1');
scale = foo.p1;
offset = foo.p2;

calibrated_features = (scale .* (pc * mu_3)) + offset;
calibrated_features = pc' * calibrated_features;
calibrated_covariances = cov_3;
for cC = 1:length(calibrated_covariances); calibrated_covariances{cC} = (scale ^ 2) .* calibrated_covariances{cC}; end

X = calibrated_features;
X101 = (scale .* mu_101) + offset;
C = calibrated_covariances;

%% put things in MEASURED order
X = fliplr(X);
C = fliplr(C);
usel = fliplr(usel);

%% generate example transitions
ex_m7 = zeros(0, 6); % back 7 = B*Bext*Bext*Bext*Bext*Bext*Bext
ex_m6 = zeros(0, 6); % back 6 = B*Bext*Bext*Bext*Bext*Bext
ex_m5 = zeros(0, 6); % back 5 = B*Bext*Bext*Bext*Bext
ex_m4 = zeros(0, 6); % back 4 = B*Bext*Bext*Bext
ex_m3 = zeros(0, 6); % back 3 = B*Bext*Bext
ex_m2 = zeros(0, 6); % back 2 = B*Bext
ex_m1 = zeros(0, 6); % back 1 = B
ex_p0 = zeros(0, 6); % hold = H
ex_p1 = zeros(0, 6); % step = S
ex_p2 = zeros(0, 6); % skip 1 = K
ex_p3 = zeros(0, 6); % skip 2 = K*Kext
ex_p4 = zeros(0, 6); % skip 3 = K*Kext*Kext
ex_p5 = zeros(0, 6); % skip 4 = K*Kext*Kext*Kext
ex_p6 = zeros(0, 6); % skip 5 = K*Kext*Kext*Kext*Kext
ex_p7 = zeros(0, 6); % skip 6 = K*Kext*Kext*Kext*Kext*Kext
ex_p8 = zeros(0, 6); % skip 7 = K*Kext*Kext*Kext*Kext*Kext*Kext

% loop over levels, filling in step examples
for cL = 1:size(X, 2)
    % Back7 (-7)
    level_ii_index = cL;
    level_jj_index = cL - 7;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m7(end + 1, :) = example;
    end
    
    % Back6 (-6)
    level_ii_index = cL;
    level_jj_index = cL - 6;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m6(end + 1, :) = example;
    end
    
    % Back5 (-5)
    level_ii_index = cL;
    level_jj_index = cL - 5;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m5(end + 1, :) = example;
    end
    
    % Back4 (-4)
    level_ii_index = cL;
    level_jj_index = cL - 4;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m4(end + 1, :) = example;
    end
    
    % Back3 (-3)
    level_ii_index = cL;
    level_jj_index = cL - 3;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m3(end + 1, :) = example;
    end
    
    % Back2 (-2)
    level_ii_index = cL;
    level_jj_index = cL - 2;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m2(end + 1, :) = example;
    end
    
    % Back1 (-1)
    level_ii_index = cL;
    level_jj_index = cL - 1;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_m1(end + 1, :) = example;
    end
    
    % Hold (+0) ***special case: look too easily classifiable w/o noise, so
    % take random samples based on level covariances
    level_ii_index = cL;
    level_jj_index = cL;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example_ii = mvnrnd(X(:, level_ii_index)', C{level_ii_index});
        example_jj = mvnrnd(X(:, level_jj_index)', C{level_jj_index});
        example = [example_ii, example_jj];
        ex_p0(end + 1, :) = example;
    end
    
    % Step (+1)
    level_ii_index = cL;
    level_jj_index = cL + 1;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p1(end + 1, :) = example;
    end
    
    % sKip1 (+2)
    level_ii_index = cL;
    level_jj_index = cL + 2;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p2(end + 1, :) = example;
    end
    
    % sKip2 (+3)
    level_ii_index = cL;
    level_jj_index = cL + 3;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p3(end + 1, :) = example;
    end
    
    % sKip3 (+4)
    level_ii_index = cL;
    level_jj_index = cL + 4;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p4(end + 1, :) = example;
    end
    
    % sKip4 (+5)
    level_ii_index = cL;
    level_jj_index = cL + 5;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p5(end + 1, :) = example;
    end
    
    % sKip5 (+6)
    level_ii_index = cL;
    level_jj_index = cL + 6;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p6(end + 1, :) = example;
    end
    
    % sKip6 (+7)
    level_ii_index = cL;
    level_jj_index = cL + 7;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p7(end + 1, :) = example;
    end
    
    % sKip7 (+8)
    level_ii_index = cL;
    level_jj_index = cL + 8;
    if (level_jj_index <= size(X, 2)) && (level_jj_index > 0) && (usel(level_ii_index) == 1) && (usel(level_jj_index) == 1)
        example = [X(:, level_ii_index)', X(:, level_jj_index)'];
        ex_p8(end + 1, :) = example;
    end
end

%% set extension rates for K and B
rate_Kext = 1/8; % ==> 1/8 skips is skip2, 1/64 is skip3, ...
rate_Bext = 1/8; % ==> 1/8 backs is back2, 1/64 is back3, ...

%% determine training/validation sizes
training_size = 3000;
validation_size = 1000;

% training...
% backsteps
t_B7 = (rate_Bext ^ 6) * training_size; t_B7 = round(t_B7);
t_B6 = (rate_Bext ^ 5) * training_size; t_B6 = round(t_B6);
t_B5 = (rate_Bext ^ 4) * training_size; t_B5 = round(t_B5);
t_B4 = (rate_Bext ^ 3) * training_size; t_B4 = round(t_B4);
t_B3 = (rate_Bext ^ 2) * training_size; t_B3 = round(t_B3);
t_B2 = (rate_Bext ^ 1) * training_size; t_B2 = round(t_B2);
t_B1 = training_size - sum([t_B7 t_B6 t_B5 t_B4 t_B3 t_B2]);

% skips 
t_K7 = (rate_Kext ^ 6) * training_size; t_K7 = round(t_K7);
t_K6 = (rate_Kext ^ 5) * training_size; t_K6 = round(t_K6);
t_K5 = (rate_Kext ^ 4) * training_size; t_K5 = round(t_K5);
t_K4 = (rate_Kext ^ 3) * training_size; t_K4 = round(t_K4);
t_K3 = (rate_Kext ^ 2) * training_size; t_K3 = round(t_K3);
t_K2 = (rate_Kext ^ 1) * training_size; t_K2 = round(t_K2);
t_K1 = training_size - sum([t_K7 t_K6 t_K5 t_K4 t_K3 t_K2]);

% holds 
t_H = training_size;

% validation...
% backsteps
v_B7 = (rate_Bext ^ 6) * validation_size; v_B7 = round(v_B7);
v_B6 = (rate_Bext ^ 5) * validation_size; v_B6 = round(v_B6);
v_B5 = (rate_Bext ^ 4) * validation_size; v_B5 = round(v_B5);
v_B4 = (rate_Bext ^ 3) * validation_size; v_B4 = round(v_B4);
v_B3 = (rate_Bext ^ 2) * validation_size; v_B3 = round(v_B3);
v_B2 = (rate_Bext ^ 1) * validation_size; v_B2 = round(v_B2);
v_B1 = validation_size - sum([v_B7 v_B6 v_B5 v_B4 v_B3 v_B2]);

% skips 
v_K7 = (rate_Kext ^ 6) * validation_size; v_K7 = round(v_K7);
v_K6 = (rate_Kext ^ 5) * validation_size; v_K6 = round(v_K6);
v_K5 = (rate_Kext ^ 4) * validation_size; v_K5 = round(v_K5);
v_K4 = (rate_Kext ^ 3) * validation_size; v_K4 = round(v_K4);
v_K3 = (rate_Kext ^ 2) * validation_size; v_K3 = round(v_K3);
v_K2 = (rate_Kext ^ 1) * validation_size; v_K2 = round(v_K2);
v_K1 = validation_size - sum([v_K7 v_K6 v_K5 v_K4 v_K3 v_K2]);

% holds 
v_H = validation_size;

%% Step vs sKip
% build training/validation sets
% Steps
ix = randperm(size(ex_p1, 1));
t_ix = ix(1:training_size);
v_ix = ix((training_size + 1):(training_size + validation_size));
t_ex_p1 = ex_p1(t_ix, :);
v_ex_p1 = ex_p1(v_ix, :);

% Skip7
ix = randperm(size(ex_p8, 1));
t_ix = ix(1:t_K7);
v_ix = ix((t_K7 + 1):(t_K7 + v_K7));
t_ex_p8 = ex_p8(t_ix, :);
v_ex_p8 = ex_p8(v_ix, :);

% Skip6
ix = randperm(size(ex_p7, 1));
t_ix = ix(1:t_K6);
v_ix = ix((t_K6 + 1):(t_K6 + v_K6));
t_ex_p7 = ex_p7(t_ix, :);
v_ex_p7 = ex_p7(v_ix, :);

% Skip5
ix = randperm(size(ex_p6, 1));
t_ix = ix(1:t_K5);
v_ix = ix((t_K5 + 1):(t_K5 + v_K5));
t_ex_p6 = ex_p6(t_ix, :);
v_ex_p6 = ex_p6(v_ix, :);

% Skip4
ix = randperm(size(ex_p5, 1));
t_ix = ix(1:t_K4);
v_ix = ix((t_K4 + 1):(t_K4 + v_K4));
t_ex_p5 = ex_p5(t_ix, :);
v_ex_p5 = ex_p5(v_ix, :);

% Skip3
ix = randperm(size(ex_p4, 1));
t_ix = ix(1:t_K3);
v_ix = ix((t_K3 + 1):(t_K3 + v_K3));
t_ex_p4 = ex_p4(t_ix, :);
v_ex_p4 = ex_p4(v_ix, :);

% Skip2
ix = randperm(size(ex_p3, 1));
t_ix = ix(1:t_K2);
v_ix = ix((t_K2 + 1):(t_K2 + v_K2));
t_ex_p3 = ex_p3(t_ix, :);
v_ex_p3 = ex_p3(v_ix, :);

% Skip1
ix = randperm(size(ex_p2, 1));
t_ix = ix(1:t_K1);
v_ix = ix((t_K1 + 1):(t_K1 + v_K1));
t_ex_p2 = ex_p2(t_ix, :);
v_ex_p2 = ex_p2(v_ix, :);

% all skips
t_ex_step = t_ex_p1;
v_ex_step = v_ex_p1;
t_ex_skip = [t_ex_p8 ; t_ex_p7 ; t_ex_p6 ; t_ex_p5 ; t_ex_p4 ; t_ex_p3 ; t_ex_p2];
v_ex_skip = [v_ex_p8 ; v_ex_p7 ; v_ex_p6 ; v_ex_p5 ; v_ex_p4 ; v_ex_p3 ; v_ex_p2];

% train SVM (Quadratic Polynomial Kernel)
X_train = [t_ex_step ; t_ex_skip];
y_train = [ones(size(t_ex_step, 1), 1) ; -1 .* ones(size(t_ex_skip, 1), 1)];
T = table(X_train, y_train);
[svm.step_v_skip, val.step_v_skip] = trainClassifierQuadSVM(T);

% logit on validation set
logit = @(x, params) glmval(params, x, 'logit');
[labels, scores] = predict(svm.step_v_skip.ClassificationSVM, [v_ex_step ; v_ex_skip]);
out = scores(:, 2);
target = [true(size(v_ex_step, 1), 1) ; false(size(v_ex_skip, 1), 1)];
params.step_v_skip = glmfit(out, target, 'binomial');

% evaluate accuracy on validation set
validation.step_v_skip.corr_rate_S = sum(out(1:size(v_ex_step, 1)) >= 0) / size(v_ex_step, 1);
validation.step_v_skip.corr_rate_K = sum(out((size(v_ex_step, 1) + 1):(size(v_ex_step, 1) + size(v_ex_skip, 1))) <= 0) / size(v_ex_skip, 1);

%% Step vs Back
% build training/validation sets
% Steps 
ix = randperm(size(ex_p1, 1));
t_ix = ix(1:training_size);
v_ix = ix((training_size + 1):(training_size + validation_size));
t_ex_p1 = ex_p1(t_ix, :);
v_ex_p1 = ex_p1(v_ix, :);

% Back7
ix = randperm(size(ex_m7, 1));
t_ix = ix(1:t_B7);
v_ix = ix((t_B7 + 1):(t_B7 + v_B7));
t_ex_m7 = ex_m7(t_ix, :);
v_ex_m7 = ex_m7(v_ix, :);

% Back6
ix = randperm(size(ex_m6, 1));
t_ix = ix(1:t_B6);
v_ix = ix((t_B6 + 1):(t_B6 + v_B6));
t_ex_m6 = ex_m6(t_ix, :);
v_ex_m6 = ex_m6(v_ix, :);

% Back5
ix = randperm(size(ex_m5, 1));
t_ix = ix(1:t_B5);
v_ix = ix((t_B5 + 1):(t_B5 + v_B5));
t_ex_m5 = ex_m5(t_ix, :);
v_ex_m5 = ex_m5(v_ix, :);

% Back4
ix = randperm(size(ex_m4, 1));
t_ix = ix(1:t_B4);
v_ix = ix((t_B4 + 1):(t_B4 + v_B4));
t_ex_m4 = ex_m4(t_ix, :);
v_ex_m4 = ex_m4(v_ix, :);

% Back3
ix = randperm(size(ex_m3, 1));
t_ix = ix(1:t_B3);
v_ix = ix((t_B3 + 1):(t_B3 + v_B3));
t_ex_m3 = ex_m3(t_ix, :);
v_ex_m3 = ex_m3(v_ix, :);

% Back2
ix = randperm(size(ex_m2, 1));
t_ix = ix(1:t_B2);
v_ix = ix((t_B2 + 1):(t_B2 + v_B2));
t_ex_m2 = ex_m2(t_ix, :);
v_ex_m2 = ex_m2(v_ix, :);

% Back1
ix = randperm(size(ex_m1, 1));
t_ix = ix(1:t_B1);
v_ix = ix((t_B1 + 1):(t_B1 + v_B1));
t_ex_m1 = ex_m1(t_ix, :);
v_ex_m1 = ex_m1(v_ix, :);

% all backs
t_ex_step = t_ex_p1;
v_ex_step = v_ex_p1;
t_ex_back = [t_ex_m7 ; t_ex_m6 ; t_ex_m5 ; t_ex_m4 ; t_ex_m3 ; t_ex_m2 ; t_ex_m1];
v_ex_back = [v_ex_m7 ; v_ex_m6 ; v_ex_m5 ; v_ex_m4 ; v_ex_m3 ; v_ex_m2 ; v_ex_m1];

% train SVM (Quadratic Polynomial Kernel)
X_train = [t_ex_step ; t_ex_back];
y_train = [ones(size(t_ex_step, 1), 1) ; -1 .* ones(size(t_ex_back, 1), 1)];
T = table(X_train, y_train);
[svm.step_v_back, val.step_v_back] = trainClassifierQuadSVM(T);

% logit on validation set
logit = @(x, params) glmval(params, x, 'logit');
[labels, scores] = predict(svm.step_v_back.ClassificationSVM, [v_ex_step ; v_ex_back]);
out = scores(:, 2);
target = [true(size(v_ex_step, 1), 1) ; false(size(v_ex_back, 1), 1)];
params.step_v_back = glmfit(out, target, 'binomial');

% evaluate accuracy on validation set
validation.step_v_back.corr_rate_S = sum(out(1:size(v_ex_step, 1)) >= 0) / size(v_ex_step, 1);
validation.step_v_back.corr_rate_B = sum(out((size(v_ex_step, 1) + 1):(size(v_ex_step, 1) + size(v_ex_back, 1))) <= 0) / size(v_ex_back, 1);

%% Step vs Hold
% build training/validation sets
% Steps
ix = randperm(size(ex_p1, 1));
t_ix = ix(1:training_size);
v_ix = ix((training_size + 1):(training_size + validation_size));
t_ex_p1 = ex_p1(t_ix, :);
v_ex_p1 = ex_p1(v_ix, :);

% Holds
ix = randperm(size(ex_p0, 1));
t_ix = ix(1:t_H);
v_ix = ix((t_H + 1):(t_H + v_H));
t_ex_p0 = ex_p0(t_ix, :);
v_ex_p0 = ex_p0(v_ix, :);

% all holds
t_ex_step = t_ex_p1;
v_ex_step = v_ex_p1;
t_ex_hold = t_ex_p0;
v_ex_hold = v_ex_p0;

% train SVM (Quadratic Polynomial Kernel)
X_train = [t_ex_step ; t_ex_hold];
y_train = [ones(size(t_ex_step, 1), 1) ; -1 .* ones(size(t_ex_hold, 1), 1)];
T = table(X_train, y_train);
[svm.step_v_hold, val.step_v_hold] = trainClassifierQuadSVM(T);

% logit on validation set
logit = @(x, params) glmval(params, x, 'logit');
[labels, scores] = predict(svm.step_v_hold.ClassificationSVM, [v_ex_step ; v_ex_hold]);
out = scores(:, 2);
target = [true(size(v_ex_step, 1), 1) ; false(size(v_ex_hold, 1), 1)];
params.step_v_hold = glmfit(out, target, 'binomial');

% evaluate accuracy on validation set
validation.step_v_hold.corr_rate_S = sum(out(1:size(v_ex_step, 1)) >= 0) / size(v_ex_step, 1);
validation.step_v_hold.corr_rate_H = sum(out((size(v_ex_step, 1) + 1):(size(v_ex_step, 1) + size(v_ex_hold, 1))) <= 0) / size(v_ex_hold, 1);

%% Add in some info to svms
svm.step_v_skip.n_skip = [t_K7 t_K6 t_K5 t_K4 t_K3 t_K2 t_K1];
svm.step_v_skip.name_skip = [8 7 6 5 4 3 2];
svm.step_v_skip.file_of_origin = 'buildRecombinationSVM.m';

svm.step_v_back.n_back = [t_B7 t_B6 t_B5 t_B4 t_B3 t_B2 t_B1];
svm.step_v_back.name_back = [-7 -6 -5 -4 -3 -2 -1];
svm.step_v_back.file_of_origin = 'buildRecombinationSVM.m';

svm.step_v_hold.n_hold = [t_H];
svm.step_v_hold.name_hold = [0];
svm.step_v_hold.file_of_origin = 'buildRecombinationSVM.m';

%% Save Products
% save('svm_for_recombination_filter.mat', 'svm');
% save('logit_for_backstep_filter.mat', 'logit');
% save('logit_params_for_recombination_filtering.mat', 'params');
% save('validation_for_recombination_filter.mat', 'validation');

