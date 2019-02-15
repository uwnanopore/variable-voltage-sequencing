function [c_data,levels,up_levels,down_levels,levels_nocc,up_levels_nocc,down_levels_nocc,phase,phase_index, Xraw_3, Craw_3, Kraw_3, well_conditioned, pc, evaluation_voltage] = ...
    capacitanceCompensation(v_data,i_data,transitions,ac_period,low_voltage,high_voltage,sample_frequency)

% turn off bad conditioning warnings
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:illConditionedMatrix');

% load default pcs and evaluation voltages for cap comp 3pc calculation
pc = load('principal_components.mat');
pc = pc.principal_components;

spacing = load('equal_x_spacing.mat');
spacing = spacing.spacing;
evaluation_voltage = spacing.v';

% calculate phase for event based on voltage data
phase = angle(sum(v_data.*exp(2*pi*1i/ac_period*(1:length(v_data)))));
phase = phase*ac_period/(2*pi);

% constrain the phase to be in domain [0,2*ac_period)
if phase < 0
    phase = phase + ac_period;
end

% round phase to nearest integer
phase = round(phase);

%find phase value for each point in data time series
phase_index = 1:length(v_data);
phase_index = phase_index - phase;
phase_index = mod(phase_index,ac_period);
phase_index = phase_index+1;

% initialize levels output in event_result
levels = cell(1,length(transitions)-1);
up_levels = cell(1,length(transitions)-1);
down_levels = cell(1,length(transitions)-1);
levels_nocc = cell(1,length(transitions)-1);
up_levels_nocc = cell(1,length(transitions)-1);
down_levels_nocc = cell(1,length(transitions)-1);
c_data = nan(1,length(i_data));

% initialize pca output
Xraw_3 = nan(3, length(transitions) - 1);
Craw_3 = cell(1, length(transitions) - 1);
Kraw_3 = cell(1, length(transitions) - 1);
well_conditioned = nan(1, length(transitions) - 1);

% loop through all found levels and conduct compensation
temp_iterator = 0;
insufficient_data = [];
for cL = 1:(length(transitions)-1)
    temp_i_data = i_data(transitions(cL)+1:(transitions(cL+1)));
    temp_v_data = v_data(transitions(cL)+1:(transitions(cL+1)));
    temp_phase_data = phase_index(transitions(cL)+1:(transitions(cL+1)));
    
    % grab only data where voltage sweep is in range
    in_sweep = (temp_v_data >= low_voltage) & (temp_v_data <= high_voltage);
    temp_v_data(~in_sweep) = nan;
    temp_i_data(~in_sweep) = nan;

    % examine to make sure the in-sweep region is at least one period long
    % if not, print warning
    if sum(~isnan(temp_v_data)) < ac_period
        insufficient_data = [insufficient_data cL];
        fprintf('WARNING: level %d of %d does not have a full cycle in bounds.\n',cL,length(transitions)-1);
        continue
    end
    
    % do cap comp
    [temp_c_data, iv_data, temp_Xraw_3, temp_Craw_3, temp_Kraw_3, temp_well_conditioned] = ...
        doCapacitanceCompensation(temp_v_data, temp_i_data, temp_phase_data, 'f_samp', sample_frequency, 'f_ac', sample_frequency / ac_period, 'pc', pc, 'evaluation', evaluation_voltage);
    
    % fill in event_result output
    levels{cL} = iv_data.levels;
    up_levels{cL} = iv_data.up_levels;
    down_levels{cL} = iv_data.down_levels;
    levels_nocc{cL} = iv_data.levels_nocc;
    up_levels_nocc{cL} = iv_data.up_levels_nocc;
    down_levels_nocc{cL} = iv_data.down_levels_nocc;
    
    % fill in pca output
    Xraw_3(:, cL) = temp_Xraw_3;
    Craw_3{cL} = temp_Craw_3;
    Kraw_3{cL} = temp_Kraw_3;
    well_conditioned(cL) = temp_well_conditioned;
    
    % fill in event_result.c_data for this level
    c_data((temp_iterator+1):(temp_iterator+length(temp_c_data))) = temp_c_data;
    temp_iterator = temp_iterator + length(temp_c_data);
end
    
% remove entries with insufficient data
levels(insufficient_data) = [];
up_levels(insufficient_data) = [];
down_levels(insufficient_data) = [];
levels_nocc(insufficient_data) = [];
up_levels_nocc(insufficient_data) = [];
down_levels_nocc(insufficient_data) = [];

Xraw_3(:, insufficient_data) = [];
Craw_3(insufficient_data) = [];
Kraw_3(insufficient_data) = [];
well_conditioned(insufficient_data) = [];

% correct ill-conditioned covariance/stiffness measurements by setting them
% to the Covariance with the 90th percentile determinant
good_C = Craw_3(well_conditioned == 1);
dets = zeros(1, size(good_C, 2));
for cC = 1:size(good_C, 2)
    dets(cC) = det(good_C{cC});
end
dets90 = prctile(dets, 90);
[~, dets90loc] = min(abs(dets - dets90));
C90 = good_C{dets90loc};
Craw_3(well_conditioned == 0) = {C90};
K90 = inv(C90);
Kraw_3(well_conditioned == 0) = {K90};

% turn back on bad conditioning warnings
warning('on', 'MATLAB:nearlySingularMatrix');
warning('on', 'MATLAB:singularMatrix');

end

function [compData, ivData, Xraw_3, Craw_3, Kraw_3, well_conditioned] = doCapacitanceCompensation(v_data, i_data, phase_ix, varargin)

% FUNCTIONALITY:
%   subroutine to implement the cap comp algorithm on individual enzyme
%   step data.
%   [compData, ivData, Xraw_3, Craw_3, Kraw_3, well_conditioned] = capacitanceCompensationNew(v_data, i_data, phase_ix, varargin)
%   Takes in raw time series data for a single level (pre capacitive
%   correction). Computes a correction function for the bilayer capacitance, 
%   and applies this correction to the current data. Additionally,
%   calculates the 3pc decomposition for each indiviual voltage sweep in
%   the level (post cap-comp) and the covariance/stiffness of these 3pc
%   measurements.
%   
% INPUTS:
%   v_data: 1 x Npts array of time series voltage values for level
%   i_data: 1 x Npts array of time series current values for level
%   phase_ix: 1 x Npts array of phase indices foro level time series
%   
% VARARGIN:
%   f_samp: sampling frequency, default = 5e4
%   f_ac: ac cycle frequency, default = 200
%   pc: principle components to use in decomposition to 3 pcs. default is
%       normalized AC princomps from MatlabCode
%   evaluation: 101 x 1 array of voltage points (in mV) at which to
%               evaluate interpolation for the 101 feature vectors. default is from
%               equal_x_spacing_160915 in ACnew
% 
% OUTPUTS:
%   compData: 1 x Npts array of cap-comped current time series
%   ivData: structure containing detailed info for (reg, up,
%           down)_levels(cc, nocc)
%   Xraw_3: 3 x 1 array with mean value of 3pcs for all cycles
%   Craw_3: 3 x 3 array with covariance of 3pcs calculated from all cycles
%   Kraw_3: 3 x 3 array with stiffness of 3pcs calculated from all cycles
%   well_conditioned: true if we had enough cycles to accurately estimate
%                     Craw (need 6 or more), false if we didn't have enough
%   

%% defaults and varargin
f_samp = 5e4;
f_ac = 200;

for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'f_samp'
            f_samp = varargin{cV + 1};
        case 'f_ac'
            f_ac = varargin{cV + 1};
        case 'pc'
            pc = varargin{cV + 1};
        case 'evaluation'
            evaluation_voltage = varargin{cV + 1};
    end
end

% calculate period T
T = round(f_samp / f_ac);

% grab starting length of voltage data
initialLength = length(v_data);

% initialize ivData structure
ivData = struct;

% initialize ivData levels structure
ivData.levels = struct;
ivData.levels.numPts = nan(1, round(T / 2));
ivData.levels.avgV = nan(1, round(T / 2));
ivData.levels.stdV = nan(1, round(T / 2));
ivData.levels.semV = nan(1, round(T / 2));
ivData.levels.avgI = nan(1, round(T / 2));
ivData.levels.stdI = nan(1, round(T / 2));
ivData.levels.semI = nan(1, round(T / 2));

% initialize up_levels, down_levels, levels_nocc, up_levels_nocc, and
% down_levels_nocc to match levels
ivData.up_levels = ivData.levels;
ivData.down_levels = ivData.levels;
ivData.levels_nocc = ivData.levels;
ivData.up_levels_nocc = ivData.levels;
ivData.down_levels_nocc = ivData.levels;

%% check to make sure have enough non NaN entries
% if we don't have at least a full period of non-nan entries, nan out the
% outputs and continue
if sum(~isnan(v_data)) < T
    compData = nan(1, initialLength);
    ivData = nan;
    Xraw_3 = nan(3, 1);
    Craw_3 = {nan(3, 3)};
    Kraw_3 = {nan(3, 3)};
    well_conditioned = nan;
    return
end

%% check to make sure level is not on DC
if std(v_data, 'omitnan') < 5
    compData = nan(1, initialLength);
    ivData = nan;
    Xraw_3 = nan(3, 1);
    Craw_3 = {nan(3, 3)};
    Kraw_3 = {nan(3, 3)};
    well_conditioned = nan;
    return
end

%% make into full cycles by appending NaNs to front and back
v_data = [nan * (1:1:phase_ix(1) - 1) v_data nan * (phase_ix(end) + 1:1:T)];
i_data = [nan * (1:1:phase_ix(1) - 1) i_data nan * (phase_ix(end) + 1:1:T)];

% find how many cycles we have
cyc = round(length(v_data) / T);

% reshape v_data into a matrix where each column is the data for one
% complete cycle
v_mat = reshape(v_data, T, cyc);
% do the same for i_data
i_mat = reshape(i_data, T, cyc);

% create long matrices where each row now corresponds to a single voltage,
% as opposed to a single phase point (each voltage point corresponds to two
% phase points, one for up, one for down)
v_mat_long = [v_mat(1:round(T / 2), :) , flipud(v_mat(round(T / 2) + 1:T, :))];
i_mat_long = [i_mat(1:round(T / 2), :) , flipud(i_mat(round(T / 2) + 1:T, :))];

%% fill in what you can so far
% numPts is a 1 x num_voltages array of the number of non-nan measurements
% at each voltage value
ivData.levels.numPts = fliplr(sum(~isnan(v_mat_long), 2)');

% avgV is the 1 x num_voltages array of average voltage value at each
% voltage point
ivData.levels.avgV = fliplr(median(v_mat_long, 2, 'omitnan')');

% stdV is the 1 x num_voltages array of voltage standard deviations at each
% voltage point
ivData.levels.stdV = fliplr(std(v_mat_long, 0, 2, 'omitnan')');

% semV is the 1 x num_voltages array of standard error of the mean at each
% voltage point
ivData.levels.semV = ivData.levels.stdV ./ sqrt(ivData.levels.numPts);

% Fields here correspond to fields for .levels, but only for up-slope
% regions
ivData.up_levels.numPts = fliplr(sum(~isnan(v_mat(round(T / 2) + 1:T,:)), 2)');
ivData.up_levels.avgV = fliplr(median(v_mat(round(T / 2) + 1:T,:), 2, 'omitnan')');
ivData.up_levels.stdV = fliplr(std(v_mat(round(T / 2) + 1:T,:), 0, 2, 'omitnan')');
ivData.up_levels.semV = ivData.up_levels.stdV ./ sqrt(ivData.up_levels.numPts);

% Same for down-slope regions
ivData.down_levels.numPts = fliplr(sum(~isnan(v_mat(1:round(T / 2),:)), 2)');
ivData.down_levels.avgV = fliplr(median(v_mat(1:round(T / 2), :), 2, 'omitnan')');
ivData.down_levels.stdV = fliplr(std(v_mat(1:round(T / 2), :), 0, 2, 'omitnan')');
ivData.down_levels.semV = ivData.down_levels.stdV ./ sqrt(ivData.down_levels.numPts);

% _nocc fields have same voltage information, so just copy it directly
ivData.levels_nocc = ivData.levels;
ivData.up_levels_nocc = ivData.up_levels;
ivData.down_levels_nocc = ivData.down_levels;

% can also fill in current values for _nocc fields
% .avgI is median current by voltage point
ivData.levels_nocc.avgI = fliplr(median(i_mat_long, 2, 'omitnan')');
% standard deviation of current by voltage point
ivData.levels_nocc.stdI = fliplr(std(i_mat_long, 0, 2, 'omitnan')');
% standard error of mean of current by voltage point
ivData.levels_nocc.semI = ivData.levels_nocc.stdI ./ sqrt(ivData.levels_nocc.numPts);

% same for only up slopes
ivData.up_levels_nocc.avgI = fliplr(median(i_mat(round(T / 2) + 1:T, :), 2, 'omitnan')');
ivData.up_levels_nocc.stdI = fliplr(std(i_mat(round(T / 2) + 1:T, :), 0, 2, 'omitnan')');
ivData.up_levels_nocc.semI = ivData.up_levels_nocc.stdI ./ sqrt(ivData.up_levels_nocc.numPts);

% same for only down slopes
ivData.down_levels_nocc.avgI = fliplr(median(i_mat(1:round(T / 2), :), 2, 'omitnan')');
ivData.down_levels_nocc.stdI = fliplr(std(i_mat(1:round(T / 2), :), 0, 2, 'omitnan')');
ivData.down_levels_nocc.semI = ivData.down_levels_nocc.stdI ./ sqrt(ivData.down_levels_nocc.numPts);

%% compute matrix averages
% average voltage by voltage point is median along the 2nd dim of long
% matrix
v_avg_all = median(v_mat_long, 2, 'omitnan')';

% up slope entries are contained in the right half of long matrix
i_avg_up = median(i_mat(round(T / 2) + 1:T, :), 2, 'omitnan')';
% down slope entries are contained in the left half of long matrix
i_avg_do = median(i_mat(1:round(T / 2), :), 2, 'omitnan')';

%% calculate the correction function from the averaged data
% grab voltage, upslope current, and downslope cur7rent, and put all three
% in the same ordering
V = fliplr(v_avg_all);
i_up = i_avg_up;
i_do = fliplr(i_avg_do);

% H_V is the difference between down and up slopes
H_V = i_do - i_up;

% vq1, vq2, vq3 are the first, second, and third voltage quartiles
vq1 = median(V(V < median(V, 'omitnan')), 'omitnan');
vq2 = median(V, 'omitnan');
vq3 = median(V(V > median(V, 'omitnan')), 'omitnan');

% q0ix, q1ix, q2ix, q3ix, q4ix are the indices of the 0th, 1st, 2nd, 3rd,
% 4th voltage quartiles
q0ix = 1;
[~, q1ix] = min(abs(V - vq1));
[~, q2ix] = min(abs(V - vq2));
[~, q3ix] = min(abs(V - vq3));
q4ix = length(V);

% calculate overall offset using quadratic fit to the middle section of H_V
[X, Y] = fitPrep2(V(q1ix:q3ix), H_V(q1ix:q3ix));
% enforce vertex of fit to occur at voltage midpoint
vertex_pt = vq2;
Xmat = zeros(size(X, 1), 2);
Xmat(:, 1) = (X - vertex_pt) .^ 2;
Xmat(:, 2) = ones(size(X, 1), 1);
out_fit = inv((Xmat' * Xmat)) * Xmat' * Y;
m = out_fit(2);

% initialize C_up and C_do--the corrector functions for up and down slope
C_up = zeros(1, length(V));
C_do = zeros(1, length(V));

% calculate left half of C_up and C_do
C_up(q0ix:q2ix) = H_V(q0ix:q2ix) - (m / 2);
C_do(q0ix:q2ix) = m / 2;

% calculate right half of C_up and C_do
C_up(q2ix + 1:q4ix) = m / 2;
C_do(q2ix + 1:q4ix) = H_V(q2ix + 1:q4ix) - (m / 2);

%% apply the correction function to all of the data
% put correction function into the right order to be added into current
% matrix directly
C_do = -fliplr(C_do);
corrector = [C_do C_up]';
corrector = repmat(corrector, 1, size(i_mat, 2));

% generate comp_mat (matrix of compensated data) by directly adding
% correction matrix into the data matrix
comp_mat = i_mat + corrector;

% reshape comp_mat to go by voltage points rather than phase index
comp_mat_long = [comp_mat(1:round(T / 2), :) , flipud(comp_mat(round(T / 2) + 1:T, :))];

%% HERE IS WHERE TO DO PCA STUFF
comp_mat_for_pca = flipud(comp_mat_long);
v_mat_for_pca = flipud(v_mat_long);
pc_mat = zeros(101, size(comp_mat_for_pca, 2));
for cI = 1:size(pc_mat, 2)
    if sum(isnan(comp_mat_for_pca(:, cI)) > 0)
        pc_mat(:, cI) = nan(101, 1);
        continue
    end
    temp_i_avg = comp_mat_for_pca(:, cI);
    temp_v_avg = v_mat_for_pca(:, cI);
    temp_g_avg = temp_i_avg ./ temp_v_avg;
    [temp_v_fit, temp_g_fit] = fitPrep2(temp_v_avg, temp_g_avg);
    interpolation = interp1q(temp_v_fit, temp_g_fit, evaluation_voltage);
    pc_mat(:, cI) = interpolation;
end

[~, badix] = filterSpikes(nanmean(pc_mat,2));
goods = ~badix;
pc_mat(badix, :) = [];
pc_mat = ((pc(goods,:)'*pc(goods,:))\pc(goods,:)') * pc_mat;
pc_mat(:, sum(isnan(pc_mat), 1) > 0) = [];
Xraw_3 = mean(pc_mat, 2);
Craw_3 = cov(pc_mat');
Kraw_3 = inv(Craw_3);

if size(pc_mat, 2) < 6
    well_conditioned = false;
else
    well_conditioned = true;
end

%% Filling in outputs
% can now fill in avg, std, and sem for currents in compensated structures
ivData.levels.avgI = fliplr(median(comp_mat_long,2,'omitnan')');
ivData.levels.stdI = fliplr(std(comp_mat_long,0,2,'omitnan')');
ivData.levels.semI = ivData.levels.stdI./sqrt(ivData.levels.numPts);

% upslopes are the bottom half of the comp_mat
ivData.up_levels.avgI = fliplr(median(comp_mat(round(T/2)+1:T,:),2,'omitnan')');
ivData.up_levels.stdI = fliplr(std(comp_mat(round(T/2)+1:T,:),0,2,'omitnan')');
ivData.up_levels.semI = ivData.up_levels.stdI./sqrt(ivData.up_levels.numPts);

% downslopes are the top half of the comp_mat
ivData.down_levels.avgI = fliplr(median(comp_mat(1:round(T/2),:),2,'omitnan')');
ivData.down_levels.stdI = fliplr(std(comp_mat(1:round(T/2),:),0,2,'omitnan')');
ivData.down_levels.semI = ivData.down_levels.stdI./sqrt(ivData.down_levels.numPts);

%%
% this entire block only serves to change any STDs or SEMs that == 0 to
% NaNs
% these 0 STD/SEM entries are result of single meaurements of given point
ivData.levels.stdV(ivData.levels.stdV==0) = nan;
ivData.levels.semV(ivData.levels.semV==0) = nan;
ivData.levels.stdI(ivData.levels.stdI==0) = nan;
ivData.levels.semI(ivData.levels.semI==0) = nan;

ivData.up_levels.stdV(ivData.up_levels.stdV==0) = nan;
ivData.up_levels.semV(ivData.up_levels.semV==0) = nan;
ivData.up_levels.stdI(ivData.up_levels.stdI==0) = nan;
ivData.up_levels.semI(ivData.up_levels.semI==0) = nan;

ivData.down_levels.stdV(ivData.down_levels.stdV==0) = nan;
ivData.down_levels.semV(ivData.down_levels.semV==0) = nan;
ivData.down_levels.stdI(ivData.down_levels.stdI==0) = nan;
ivData.down_levels.semI(ivData.down_levels.semI==0) = nan;

ivData.levels_nocc.stdV(ivData.levels_nocc.stdV==0) = nan;
ivData.levels_nocc.semV(ivData.levels_nocc.semV==0) = nan;
ivData.levels_nocc.stdI(ivData.levels_nocc.stdI==0) = nan;
ivData.levels_nocc.semI(ivData.levels_nocc.semI==0) = nan;

ivData.up_levels_nocc.stdV(ivData.up_levels_nocc.stdV==0) = nan;
ivData.up_levels_nocc.semV(ivData.up_levels_nocc.semV==0) = nan;
ivData.up_levels_nocc.stdI(ivData.up_levels_nocc.stdI==0) = nan;
ivData.up_levels_nocc.semI(ivData.up_levels_nocc.semI==0) = nan;

ivData.down_levels_nocc.stdV(ivData.down_levels_nocc.stdV==0) = nan;
ivData.down_levels_nocc.semV(ivData.down_levels_nocc.semV==0) = nan;
ivData.down_levels_nocc.stdI(ivData.down_levels_nocc.stdI==0) = nan;
ivData.down_levels_nocc.semI(ivData.down_levels_nocc.semI==0) = nan;

%%
% reshape comp_mat into a 1xNpts array for output as compData
compData = reshape(comp_mat,1,size(i_mat,1)*size(i_mat,2));

% remove NaN entries from compData
compData(isnan(compData)) = [];

end
    










































