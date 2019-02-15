processing_parameters = struct;
% general
processing_parameters.general.sampling_frequency = 50000;
processing_parameters.general.voltage_cycle_frequency = 200;
processing_parameters.general.cycle_period = 250;
processing_parameters.general.voltage_range = [100, 200];
% flicker filtering
processing_parameters.filterFlickers.pflicker = [3, 3];
processing_parameters.filterFlickers.windows = [2, 5];
% change point detection
processing_parameters.levelFinding.minimum_level_length = 250;
processing_parameters.levelFinding.sensitivity = 4;
% capacitance compensation
processing_parameters.capacitanceCompensation.low_voltage = 95;
processing_parameters.capacitanceCompensation.high_voltage = 205;
% removal filter
processing_parameters.removalFilter.good_level_threshold = 0.35;
% recombination filter
processing_parameters.recombinationFilter.enzyme_step_probabilities_names = 'step | skip | skip extend | back | back extend | hold';
processing_parameters.recombinationFilter.enzyme_step_probabilities = [0.75, 0.05, 0.05, 0.125, 0.125, 0.075];
processing_parameters.recombinationFilter.backstep_lookback = 16;
% reordering filter
processing_parameters.reorderFilter.prior = [.5, .1, .05]';
processing_parameters.reorderFilter.max_iterations = 5;

% saving
saveon = false;
if saveon
    save('C:/Users/nanopore/Desktop/DataDeposition/processing_parameters_vv.mat', 'processing_parameters');
end