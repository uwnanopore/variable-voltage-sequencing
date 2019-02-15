function [f_to_tf, x_101_tf, k_101_tf, x_3_tf, k_3_tf, c_3_tf, npts_tf, did_backstep_f, did_backstep_tf] = recombinationFilter(x_101, k_101, x_3, k_3, npts, varargin)

%% Defaults and Varargin
% defaults
% enzyme step probabilities
% p_step, p_skip, p_skip+, p_back, p_back+, p_hold
enzyme_step_probabilities = [0.75, 0.05, 0.05, 0.125, 0.125, 0.075];

% self alignment lookout distance
lookback = 8;

% principal components
pc = load('principal_components.mat');
pc = pc.principal_components;

% varargin
for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'enzymestepprobabilities'
            enzyme_step_probabilities = varargin{cV + 1};
        case 'principlecomponents'
            pc = varargin{cV + 1};
        case 'calibrationmap'
            map = varargin{cV + 1};
        case 'referencesequence'
            refseq = varargin{cV + 1};
        case 'lookback'
            lookback = varargin{cV + 1};
    end
end

%% Normalize enzyme step probabilities
total_enzyme_step_probability = ...
    enzyme_step_probabilities(1) + ...
    (enzyme_step_probabilities(2) / (1 - enzyme_step_probabilities(3))) + ...
    (enzyme_step_probabilities(4) / (1 - enzyme_step_probabilities(5))) + ...
    enzyme_step_probabilities(6);

%% Load contiguity-based step type classifiers
svm = load('svm_for_recombination_filter.mat');
svm = svm.svm;
params = load('logit_params_for_recombination_filtering.mat');
params = params.params;
logit = load('logit_for_backstep_filter.mat');
logit = logit.logit;

%% Generate Calibration Strand
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

%% Initialize placeholder variables for iterative self alignment
n_levels = size(x_101, 2);

temp_x_101 = x_101;
temp_k_101 = k_101;
temp_x_3 = x_3;
temp_k_3 = k_3;
temp_npts = npts;

total_self_alignment = 1:n_levels;
temp_self_alignment = ones(1, n_levels);

%% Iterative self-alignment filtering

% self align and recombine until self alignment makes no changes
while length(temp_self_alignment) ~= max(temp_self_alignment)
    
    % CALIBRATE TO REFERENCE -- only for step probabilities (will do
    % self-alignment in un-calibrated, un-normalized space)
    %   1, normalize measured levels by normalizerAC
    temp_n_101 = temp_x_101 - normalizerAC(temp_x_101);
    
    %   2, calculate calibration of measured to reference
    meas_cal = reshape(temp_n_101, 1, []);
    ref_cal = reshape(n_101_ref, 1, []);
    [temp_scale, temp_offset] = calibrate_by_iqr(meas_cal(~isnan(meas_cal)), ref_cal(~isnan(ref_cal)), false, [-inf inf]);
    
    %   3, apply calibration to the measured levels
    temp_n_101 = (temp_scale .* temp_n_101) + temp_offset;
    
    %   4, reduce calibrated data back to 3 pcs
    temp_n_3 = pc' * temp_n_101;
    
    % CALCULATE STEP-BY-STEP TRANSITION PENALTIES
    %   1, build transition vectors out of the normalized, calibrated 3 pcs
    %   data
    transition_vectors = nan(n_levels - 1, 6);
    transition_vectors(:, 1:3) = temp_n_3(:, 1:(end - 1))';
    transition_vectors(:, 4:6) = temp_n_3(:, 2:end)';
    
    %   2, use step/back/hold/skip svms to generate probabilities for each
    %   transition
    P = calculateStepProbabilitiesForRecombinationFilter(transition_vectors, true);
    
    %   3, tack on empty (NaN) row to the start of P to make right size
    P0 = nan(4, 1);
    P = [P0 P];
    
    %   4, build step-by-step transition probabilities
    individual_step_probabilities = nan(size(temp_x_101, 2), 6);
    individual_step_probabilities(:, 1) = P(1, :)'; % column 1 is 'step'
    individual_step_probabilities(:, 2) = P(4, :)'; % column 2 is 'skip'
    individual_step_probabilities(:, 3) = individual_step_probabilities(:, 2); % column 3 is 'skip extend', set to match skip
    individual_step_probabilities(:, 4) = P(2, :)'; % column 4 is 'back'
    individual_step_probabilities(:, 5) = individual_step_probabilities(:, 4); % column 5 is 'back extend', set to match back
    individual_step_probabilities(:, 6) = P(3, :)'; % column 6 is 'hold'
    
    %   5, fold in enzyme stepping priors
    individual_step_probabilities = individual_step_probabilities .* repmat(enzyme_step_probabilities, size(individual_step_probabilities, 1), 1);    
 
    % CALCULATE SELF ALIGNMENT PENALTIES FROM LEVEL DURATION DATA
    self_alignment_penalty = -0.5 * size(temp_x_3, 1);
    
    % DO SELF ALIGNMENT
    temp_self_alignment = selfAlignAC( ...
        temp_x_3, temp_k_3, ...
        individual_step_probabilities, ...
        self_alignment_penalty, ...
        lookback);  
    
    % initialize updated arrays
    new_x_3 = zeros(size(temp_x_3, 1), max(temp_self_alignment));
    new_k_3 = cell(1, max(temp_self_alignment));
    new_x_101 = zeros(size(temp_x_101, 1), max(temp_self_alignment));
    new_k_101 = cell(1, max(temp_self_alignment));
    new_npts = zeros(1, max(temp_self_alignment));
    
    % walk through the new levels
    for cL = 1:size(new_x_3, 2)
        
        % find original levels corresponding to each filtered level
        this_level_locations = find(temp_self_alignment == cL);
        
        % compute and record the stats of the new level
        [new_x_3(:, cL), new_k_3{cL}] = multivarSampleStats(temp_x_3(:, this_level_locations), temp_k_3(this_level_locations), true);
        [new_x_101(:, cL), new_k_101{cL}] = multivarSampleStats(temp_x_101(:, this_level_locations), temp_k_101(this_level_locations), true);
        new_npts(cL) = sum(temp_npts(this_level_locations));
    end
    
    % keep track of the alignment from initial _f to final _tf
    for cSL = 1:length(temp_self_alignment)
        total_self_alignment(total_self_alignment == cSL) = temp_self_alignment(cSL);
    end
    
    % update the data we are going to self align for the next round
    temp_x_3 = new_x_3;
    temp_k_3 = new_k_3;
    temp_x_101 = new_x_101;
    temp_k_101 = new_k_101;
    temp_npts = new_npts;
    n_levels = size(temp_x_3, 2);
end

%% finalize outputs
f_to_tf = total_self_alignment;
x_101_tf = temp_x_101;
k_101_tf = temp_k_101;
x_3_tf = temp_x_3;
k_3_tf = temp_k_3;
c_3_tf = cell(1, size(k_3_tf, 2));
for cC = 1:size(c_3_tf, 2); c_3_tf{cC} = inv(k_3_tf{cC}); end
npts_tf = temp_npts;
did_backstep_f = f_to_tf(1:end-1) > f_to_tf(2:end);
did_backstep_f = [did_backstep_f 0];
did_backstep_tf = true(1, size(x_101_tf, 2));
for cl = 1:size(x_101_tf,2)
    did_backstep_tf(cl) = any(did_backstep_f(f_to_tf == cl));
end
    
end
    
    
    
    
    
    
    
    
    
    
    