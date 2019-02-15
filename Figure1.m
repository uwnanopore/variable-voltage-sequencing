%% panel A: no data

%% panel B: raw data with error modes

% load the data we need
jumps = load('figure1_jumps_data.mat');
jumps = jumps.jumps;

reduced = load('figure1_reduced_data.mat');
reduced = reduced.reduced;

% choose the indices in reduced we want to display
ix_red = 1707100:1712150;

% extract the current data of interest
i_data = reduced.data(ix_red);

% calculate conductance data from currents
g_data = i_data ./ 180;

% grab all found level locations and level means
level_locs = cell2mat(jumps.reducedStart);

% trim level locations by which are in the region of interest
level_locs(level_locs < min(ix_red)) = [];

% i_level(level_locs > max(ix_red)) = [];
level_locs(level_locs > max(ix_red)) = [];

% tack on left and right bounds to level locs
level_locs = [min(ix_red), level_locs, max(ix_red)];

% plotting options
data_lw = .4;
data_colr = .55 .* ones(1, 3); %[.65 .65 .65];
downsample_rate = 2;

level_lw = 1.2;
level_colr = [0 0 0];

scale_lw = 3;
scale_colr = [0 0 0];

highlight_lw = 1.6;
highlight_box_h = .012;
highlight_colr = [0 1 0];
highlight_edge_colr = 'none';
highlight_transparency = 0.5;

misstep_marker = 'pentagram';
misstep_colr = [1 0 0];
misstep_markersize = 3;
marker_y_spacing = 0.012;

% PLOTTING
f1b = figure(6483); clf; hold on;

% plot highlights of where we have degen levels
left_bound = level_locs(13);
right_bound = level_locs(29);
box_y = 0.396;
rectangle('Position', [left_bound, box_y - (highlight_box_h / 2), right_bound - left_bound + 1, highlight_box_h], 'FaceColor', [highlight_colr, highlight_transparency], 'EdgeColor', highlight_edge_colr);

left_bound = level_locs(12);
right_bound = level_locs(20);
box_y = 0.324;
rectangle('Position', [left_bound, box_y - (highlight_box_h / 2), right_bound - left_bound + 1, highlight_box_h], 'FaceColor', [highlight_colr, highlight_transparency], 'EdgeColor', highlight_edge_colr);

left_bound = level_locs(14);
right_bound = level_locs(26);
box_y = 0.430;
rectangle('Position', [left_bound, box_y - (highlight_box_h / 2), right_bound - left_bound + 1, highlight_box_h], 'FaceColor', [highlight_colr, highlight_transparency], 'EdgeColor', highlight_edge_colr);

% plot the data
plot(downsampleinmatlab(ix_red, downsample_rate), downsampleinmatlab(g_data, downsample_rate), 'color', data_colr, 'linewidth', data_lw);

% plot the levels and plot markers at misstep levels
misstep_levels = [6 7 8 10 16 17];
marker_on_top = [1 0 1 0 1 0];
for cL = 1:(length(level_locs) - 1)
    left_bound = level_locs(cL);
    right_bound = level_locs(cL + 1);
    level_median = median(reduced.data(left_bound:right_bound)) ./ 180;
    plot([left_bound, right_bound], level_median .* [1 1], 'color', level_colr, 'linewidth', level_lw);
    if any(cL == misstep_levels)
        cLL = find(cL == misstep_levels);
        marker_x = (left_bound + right_bound) / 2;
        if marker_on_top(cLL) == 1
            marker_y = level_median + marker_y_spacing;
        else
            marker_y = level_median - marker_y_spacing;
        end
        plot(marker_x, marker_y, 'marker', misstep_marker, 'markersize', misstep_markersize, 'markeredgecolor', misstep_colr, 'markerfacecolor', misstep_colr);
    end
end

% plot a time scale bar
plot([1707500 1708500], [.48 .48], 'color', scale_colr, 'linewidth', scale_lw);

% set up the canvas
xlim([ix_red(1), ix_red(end)]);
ylim([0.24, 0.51]);
set(gcf, 'renderer', 'painters');
set(gca, 'XTick', [], 'YTick', [.25:.05:.50], 'YTickLabel', [.25:.05:.50]);
xlabel('Time');
ylabel('Conductance (nS)');
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 2, 'FontSize', 12);
f1b.PaperUnits = 'inches';
f1b.PaperPosition = [0 0 4 2.5];

%% panel C: constant voltage sampling
% load AC and DC maps
if ~exist('acmap', 'var') || ~exist('dcmap', 'var')
    acmap = load('pore_model_6mer_variable_voltage.mat');
    acmap = acmap.model;
    dcmap = load('pore_model_6mer_constant_voltage.mat');
    dcmap = dcmap.model;
end

% sequence to find candidate for segment
alphabet = 'ACGT';
seq = 'ACGACACTCGATAG';
[x_dc, ~, e_dc] = levelPredPipelineDC(seq, dcmap);

% calibrate the data
cal_to_g_scale = 6;
cal_to_g_offset = 0.39;

x_dc = (x_dc .* cal_to_g_scale) + cal_to_g_offset;
e_dc = (e_dc .* cal_to_g_scale);

% PLOT THE DC FIGURE
f1c = figure(10242); clf; hold on;

% plotting options
% data_colr = [0 0 0];
data_colr = [218 165 32] / 255;
data_colr_ind = [.8 .1 .1];
data_colr_dep = [.1 .1 .8];
data_marker_ind = '.';
data_markersize_ind = 18;

data_marker_dep = '.';
data_markersize_dep = 18;

stairs_style = ':';
stairs_colr = .5 .* ones(1, 3);
stairs_lw = 1;

linker_style = stairs_style;
linker_colr = stairs_colr;
linker_lw = stairs_lw;

x_border = 0.01;
y_border = 0.025;

% plot the stairs behind the measurements
for cL = 1:length(x_dc)
    plot([cL - .5, cL + .5], x_dc(cL) .* [1 1], 'linestyle', stairs_style, 'color', stairs_colr, 'linewidth', stairs_lw);
end

% plot linkers between the stairs
for cL = 1:(length(x_dc) - 1)
    plot([cL + .5, cL + .5], [x_dc(cL) x_dc(cL + 1)], 'linestyle', linker_style, 'color', linker_colr, 'linewidth', linker_lw);
end

% plot the data
for cL = 1:length(x_dc)
    if mod(cL, 2) == 1
        plot(cL, x_dc(cL), 'color', data_colr_ind, 'marker', data_marker_ind, 'markersize', data_markersize_ind);
    elseif mod(cL, 2) == 0
        plot(cL, x_dc(cL), 'color', data_colr_dep, 'marker', data_marker_dep, 'markersize', data_markersize_dep, 'markerfacecolor', [1 1 1]);
    end
end

% plot the sequence above
for cS = 1:length(seq)
    text(cS * 2 - 5.5, 0.65, seq(cS), 'fontweight', 'bold', 'color', [0 0 0], 'fontname', 'courier', 'fontsize', 12);
end

% set up the canvas
xlim([4.5 + x_border, 14.5 - x_border]);
ylim([.35 - y_border, 0.57 + y_border]);
set(gca, 'XTick', [7 8 9 10 11 12 13 14 15 16] - 2, 'XTickLabels', {'', '1', '', '2', '', '3', '', '4', '', '5'});
set(gca, 'YTick', [0.35:.05:0.60], 'YTickLabels', {'.35', '', '.45', '', '.55', ''});
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 2, 'FontSize', 12);
set(gcf, 'renderer', 'painters');
xlabel('DNA Position (nt)');
ylabel('Conductance (nS)');
f1c.PaperUnits = 'inches';
f1c.PaperPosition = [0 0 4 2.5];

%% Panel D: no data

%% Panel E: DNA position shift as a function of voltage
% load the data
stretchdata = load('stretching_fitting_result.mat');
stretchdata = stretchdata.stretching;

f1e = figure(12021); clf; hold on;

% plot shift data
shadedErrorBar(stretchdata.mV, -stretchdata.shift, stretchdata.dshift, {'-', 'color', [0 0 0], 'linewidth', 1.5}, 1);

% plot reference lines
plot([180 180], [-1.2 .2], '--', 'color', [0 0 0]);
plot([100 200], [0 0], '--', 'color', [0 0 0]);

% set up the canvas
xlim([100 200]);
ylim([-1.2 0.2]);
set(gca, 'XTick', [100 150 200]); 
set(gca, 'YTick', [-1.0 -0.5 0]);
ylabel('DNA shift (nt from 180 mV)');
xlabel('Voltage (mV)');
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 2);
set(gca, 'FontSize', 14);
set(gcf, 'renderer', 'painters');
f1e.PaperUnits = 'inches';
f1e.PaperPosition = [0 0 4 2.5];


%% Panel F: variable voltage sampling
% get ac data
[x_ac, k_ac] = levelPredPipeline(seq, acmap, true);
f1f = figure(5239); clf; hold on;
% load equal x spacing info
spacing = load('equal_x_spacing.mat');
spacing = spacing.spacing;
% calibrate ac data
x_ac = (x_ac .* cal_to_g_scale) + cal_to_g_offset;
for cK = 1:length(k_ac); k_ac{cK} = k_ac{cK} ./ (cal_to_g_scale ^ 2); end
% downsample ac features
n_ac_features = 8;
[x_ac_ds, k_ac_ds] = downsampleFeatures(x_ac, k_ac, n_ac_features);
spacing.x_ds = downsampleinmatlab(spacing.x, floor(101 / n_ac_features));
spacing.v_ds = downsampleinmatlab(spacing.v, floor(101 / n_ac_features));

% plotting options
dc_data_colr = [218 165 32] / 255;
dc_data_colr_ind = [.8 .1 .1];
dc_data_colr_dep = [.1 .1 .8];

dc_data_marker_ind = '.';
dc_data_markersize_ind = 18;

dc_data_marker_dep = '.';
dc_data_markersize_dep = 18;

data_colr_ind = [.8 .1 .1];
data_colr_dep = [.1 .1 .8];
data_ls = '-';
data_lw = 1.25;

x_border = x_border;
y_border = y_border;

% plot the dc data
for cL = 1:length(x_dc)
    if mod(cL, 2) == 1
        plot(cL, x_dc(cL), 'color', dc_data_colr_ind, 'marker', dc_data_marker_ind, 'markersize', dc_data_markersize_ind);
    elseif mod(cL, 2) == 0
        plot(cL, x_dc(cL), 'color', dc_data_colr_dep, 'marker', dc_data_marker_dep, 'markersize', dc_data_markersize_dep, 'markerfacecolor', [1 1 1]);
    end
end

% plot the ac data
for cL = 1:size(x_ac_ds, 2)
    % decide which markers to use
    if mod(cL, 2) == 1
        lw = data_lw;
        ls = data_ls;
        colr = data_colr_ind;
    elseif mod(cL, 2) == 0
        lw = data_lw;
        ls = data_ls;
        colr = data_colr_dep;
    end
    plot((-2 * spacing.x_ds(:)) + cL, x_ac_ds(:, cL), 'linestyle', ls, 'linewidth', lw, 'color', colr);
end

% set up the canvas
xlim([4.5 + x_border, 14.5 - x_border]);
ylim([.35 - y_border, 0.57 + y_border]);
set(gca, 'XTick', [7 8 9 10 11 12 13 14 15 16] - 2, 'XTickLabels', {'', '1', '', '2', '', '3', '', '4', '', '5'});
set(gca, 'YTick', [0.35:.05:0.60], 'YTickLabels', {'.35', '', '.45', '', '.55', ''});
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 2, 'FontSize', 12);
set(gcf, 'renderer', 'painters');
xlabel('DNA Position (nt)');
ylabel('Conductance (nS)');
f1f.PaperUnits = 'inches';
f1f.PaperPosition = [0 0 4 2.5];





