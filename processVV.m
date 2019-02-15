function [read] = processVV(current_data, voltage_data, is_raw)

%% set default
% if is_raw is not passed in, assume to be false as this is the case for
% all data deposited here. in this case, will jump ahead past flicker
% filtering straight to level finding
if nargin < 3
    is_raw = false;
end

%% Load the processing parameters
params = load('processing_parameters_vv.mat');
params = params.processing_parameters;

%% initialize output struct
read = struct;

%% Raw data through flicker filter
if is_raw
    % only use on RAW data (i.e. is_raw = true). in this data deposition, all
    % provided current/voltage data is post-flicker-filtering, meaning this
    % step will be skipped. this node will only be used if generating new data
    % using variable voltage + hel308 tga helicase
    
    % store raw data in the read
    read.raw_current_data = current_data;
    read.raw_voltage_data = voltage_data;
    
    % remove the fractional cycle at end of the data and any cycles with
    % voltages outside a certain range
    downsample_rate = pp.general.cycle_period;
    n_truncated = floor(length(current_data) / downsample_rate) * downsample_rate;
    truncated_raw_current_data = current_data(1:n_truncated);
    truncated_raw_voltage_data = voltage_data(1:n_truncated);
    reshaped_raw_current_data = reshape(truncated_raw_current_data, downsample_rate, []);
    reshaped_raw_voltage_data = reshape(truncated_raw_voltage_data, downsample_rate, []);
    
    % remove any voltage cycles where voltage was out of the accepted range
    voltage_range = pp.general.voltage_range;
    bad_voltage_cycles = any(reshaped_raw_voltage_data < voltage_range(1) - 5 | reshaped_raw_voltage_data > voltage_range(2) + 5);
    reshaped_raw_voltage_data(:, bad_voltage_cycles) = [];
    reshaped_raw_current_data(:, bad_voltage_cycles) = [];
    
    % make downsampled data
    downsampled_current_data = mean(reshaped_raw_current_data);
    
    % find and remove flicker cycles
    flicker_cycles = filterFlickers( ...
        downsampled_current_data, ...
        'pflicker', params.filterFlickers.pflicker, ...
        'windows', params.filterFlickers.windows);
    reshaped_raw_current_data(:, flicker_cycles) = [];
    reshaped_raw_voltage_data(:, flicker_cycles) = [];
    
    % store flicker filtered data in the read
    read.current_data = reshaped_raw_current_data(:)';
    read.voltage_data = reshaped_raw_voltage_data(:)';
    
else
    % if is_raw is false, means that the input current/voltage data has
    % already been treated for flicker filtering (as is the case for all
    % data deposited here). in this case, directly store the input data
    % into the read and proceed to level finding.
    
    read.current_data = current_data;
    read.voltage_data = voltage_data;
end

%% Find levels at enzyme steps
% load principal components for level finding
pc_levelfind = load('principal_components_for_level_finding.mat');
pc_levelfind = pc_levelfind.principal_components_for_level_finding;
pc_levelfind = pc_levelfind(:, 1:5);

% do level finding
[transitions] = findLevelsPC( ...
    read.current_data, ...
    read.voltage_data, ...
    params.general.cycle_period, ...
    pc_levelfind, ...
    'sensitivity', params.levelFinding.sensitivity, ...
    'minlevellength', params.levelFinding.minimum_level_length);

% remove too-short transitions
transitions_length = [inf diff(transitions)];
too_short_transitions = transitions_length < params.levelFinding.minimum_level_length;
transitions(too_short_transitions) = [];

% store transitions results
read.transitions = transitions;

%% Capacitance compensation
[compensated_data, levels, ~, ~, ~, ~, ~, ~, ~, Xraw_3, ~, Kraw_3, ~, ~, ~] = capacitanceCompensation( ...
    read.voltage_data, ...
    read.current_data, ...
    read.transitions, ...
    params.general.cycle_period, ...
    params.capacitanceCompensation.low_voltage, ...
    params.capacitanceCompensation.high_voltage, ...
    params.general.sampling_frequency);

% store results in read
read.compensated_data = compensated_data;
read.levels = levels;
read.x3_uf = Xraw_3;
read.k3_uf = Kraw_3;

%% Make unfiltered features
[x101_uf, k101_uf, npts_uf, ~, ~, feature_type] = makeUnfilteredFeatures( ...
    read.levels);

% store results in read
read.x101_uf = x101_uf;
read.k101_uf = k101_uf;
read.npts_uf = npts_uf;
read.feature_type = feature_type;

%% Removal filtering
[keep_ix, uf_to_f, ~] = removalFilter( ...
    read.x101_uf, ...
    'threshold', params.removalFilter.good_level_threshold);

% apply the filter and store results in read
read.x3_f = read.x3_uf(:, keep_ix == 1);
read.k3_f = read.k3_uf(keep_ix == 1);
read.x101_f = read.x101_uf(:, keep_ix == 1);
read.k101_f = read.k101_uf(keep_ix == 1);
read.npts_f = read.npts_uf(keep_ix == 1);
read.uf_to_f = uf_to_f;

%% Recombination filtering
[f_to_tf, x101_tf, k101_tf, x3_tf, k3_tf, ~, npts_tf, did_backstep_f, did_backstep_tf] = recombinationFilter( ...
    read.x101_f, ...
    read.k101_f, ...
    read.x3_f, ...
    read.k3_f, ...
    read.npts_f, ...
    'enzymestepprobabilities', params.recombinationFilter.enzyme_step_probabilities, ...
    'lookback', params.recombinationFilter.backstep_lookback);

% store results in read
read.did_backstep_f = did_backstep_f;
read.x3_tf = x3_tf;
read.k3_tf = k3_tf;
read.x101_tf = x101_tf;
read.k101_tf = k101_tf;
read.npts_tf = npts_tf;
read.f_to_tf = f_to_tf;
read.did_backstep_tf = did_backstep_tf;

%% Reordering filtering
[x3_r, k3_r, x101_r, k101_r, tf_to_r, npts_r, did_backstep_r] = reorderFilter( ...
    read.x3_tf, ...
    read.k3_tf, ...
    read.x101_tf, ...
    read.k101_tf, ...
    read.npts_tf, ...
    read.did_backstep_tf, ...
    'maxiterations', params.reorderFilter.max_iterations, ...
    'prior', params.reorderFilter.prior);

% store results in read
read.x3_r = x3_r;
read.k3_r = k3_r;
read.x101_r = x101_r;
read.k101_r = k101_r;
read.npts_r = npts_r;
read.tf_to_r = tf_to_r;
read.did_backstep_r = did_backstep_r;

end


