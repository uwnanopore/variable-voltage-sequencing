function [ff1, ff2, gg1, gg2] = Figure2(mappass)

%% load principle components
global pc
pc = load('principal_components.mat');
pc = pc.principal_components;
global map
if exist('mappass', 'var')
    map = mappass;
else
    map = load('pore_model_6mer_variable_voltage.mat');
    map = map.model;
end

%% define plotting options
global smart_skip_threshold
smart_skip_threshold = .49;

global base_colr
base_colr = 0 .* ones(1, 3); % base color for levels
global base_lw
base_lw = 2; % linewidth for levels

global inst_colr
inst_colr = .0 .* ones(1, 3); % color for aligned instances in corrected plots
global inst_lw
inst_lw = 1; % linewidth of aligned instances in corrected plots
global inst_ls
inst_ls = '--'; % linestyle of aligned instances in corrected plots

global vert_shift
vert_shift = 30; % offset between dc and ac
global under_shift
under_shift = -1; % how far to stay under dc level
global over_shift
over_shift = 1; % how far to go over ac level

trans_colr = 0 .* ones(1, 3); % color of transition lines
global trans_lw
trans_lw = 1.25; % linewidth for transition lines
global trans_shift
trans_shift = 0; % offset for transition lines

global L
L = 0.75; % size scale of markers
global LX_factor
LX_factor = .8; % augmentation of markers along x-direction
global LY_factor
LY_factor = 3.5; % augmentation of markers along y-direction
global LY_factor_square
LY_factor_square = 2.3; % different scale factors for hold markers so they don't look huge

global marker_vertloc
marker_vertloc = 69; % where to place markers vertically
global back_colr
global back_edge_colr
back_colr = [.7 .2 .2]; back_edge_colr = [1 1 1]; % color marking backsteps
global forw_colr
global forw_edge_colr
forw_colr = [.2 .7 .2]; forw_edge_colr = [1 1 1]; % color marking forward steps
global skip_colr
skip_colr = [.2 .2 .7]; skip_edge_colr = [1 1 1]; % color marking skips
global hold_colr
hold_colr = [.8 .5 .2]; hold_edge_colr = [1 1 1]; % color marking holds
global edge_width
edge_width = 2.1; % how thick to make edge of transition markers
rect_edge_width = 1; % how thick tomake edge of hold markers
global line_factor
line_factor = 1.4; % darkening factor for transition lines

global n_random_draws
n_random_draws = 100; % how many random draws from multivariate gaussian to use in estimating errors

global rescale_factor
rescale_factor = 180;

% first event
cE = 1;
roi = 146:170; 

dilate_x_thistime = 1;
[ff1, ff2] = plotPlayButton(cE, roi, 101, 10101, dilate_x_thistime, 1);

% second event
cE = 2;
roi = 58:88; 

dilate_x_thistime = 1.25;
[gg1, gg2] = plotPlayButton(cE, roi, 202, 20202, dilate_x_thistime, 0);

end


%% define a function which does all the the plotting given an event ID and a region of interest
function [ff1, ff2] = plotPlayButton(cE, roi, fignum1, fignum2, dilate_x, calshift)

if ~exist('dilate_x', 'var')
    dilate_x = 1;
end

% define global variables
global pc
global map
global n_random_draws
global rescale_factor
global smart_skip_threshold

global trans_shift
global vert_shift
global under_shift
global over_shift
global trans_lw

global back_colr
global hold_colr
global forw_colr
global skip_colr

global edge_width
global back_edge_colr
global forw_edge_colr

global marker
global marker_vertloc

global L
global LX_factor
final_lx_factor = LX_factor * dilate_x;
global LY_factor
global LY_factor_square

global line_factor

global base_lw
global base_colr

% define markers
marker = struct;
marker.p1.x = [0 0 (sqrt(3)/2)*L 0] - ((sqrt(3)/4)*L);
marker.p1.x = final_lx_factor * marker.p1.x;
marker.p1.y = [0 L L/2 0] - (L/2);
marker.p1.y = LY_factor * marker.p1.y;

marker.m1.x = -marker.p1.x;
marker.m1.y = marker.p1.y;

marker.p2.x = [0 0 (2*sqrt(3)/8)*L (2*sqrt(3)/8)*L (6*sqrt(3)/8)*L (2*sqrt(3)/8)*L (2*sqrt(3)/8)*L 0]-(2*sqrt(3)/8)*L;
marker.p2.x = final_lx_factor * marker.p2.x;
marker.p2.y = [0 L 6*L/8 L L/2 0 2*L/8 0]-L/2;
marker.p2.y = LY_factor * marker.p2.y;

marker.m2.x = -marker.p2.x;
marker.m2.y = marker.p2.y;

marker.p3.x = [0 0 (sqrt(3)/2)*L (sqrt(3)/2)*L sqrt(3)*L sqrt(3)*L (3*sqrt(3)/2)*L sqrt(3)*L sqrt(3)*L (sqrt(3)/2)*L (sqrt(3)/2)*L 0] - (3*sqrt(3)/4)*L;
marker.p3.x = final_lx_factor * marker.p3.x;
marker.p3.y = [0 L L/2 L L/2 L L/2 0 L/2 0 L/2 0] - L/2;
marker.p3.y = LY_factor * marker.p3.y;

marker.m3.x = -marker.p3.x;
marker.m3.y = marker.p3.y;

marker.p0.x = [-L/2 L/2 0 -L/2];
marker.p0.y = [(sqrt(3)/2)*L (sqrt(3)/2)*L 0 (sqrt(3)/2)*L] - (sqrt(3)/4)*L;

% ----------------------------------- $

% load event
event_result = load(['figure2_eventdata_' num2str(cE) '.mat']);
er = event_result.event_result; clear event_result

% calculate rescaling to picoamps
rescale_mean = mean(reshape(normalizerAC(er.features_r), 1, []));

% extract features and calibrate
yy_f_101 = (er.scale .* ((pc * er.Xraw_3_f) - normalizerAC(pc * er.Xraw_3_f))) + er.offset;
yy_corr_101 = (er.scale .* ((pc * er.Xraw_3_r) - normalizerAC(pc * er.Xraw_3_r))) + er.offset;

yy_f_3 = pc' * yy_f_101;
yy_corr_3 = pc' * yy_corr_101;

% extract errors and calibrate
ee_f_3 = er.Craw_3_f;
for ii = 1:length(ee_f_3)
    ee_f_3{ii} = ee_f_3{ii} .* (er.scale ^ 2);
end

ee_corr_3 = er.Craw_3_r;
for ii = 1:length(ee_corr_3)
    ee_corr_3{ii} = ee_corr_3{ii} .* (er.scale ^ 2);
end

% figure out the 101-dim errors from mvnrnd
ee_f_101 = zeros(size(yy_f_101));
for ii = 1:size(ee_f_101, 2)
    random_draws = mvnrnd(yy_f_3(:, ii)', ee_f_3{ii}, n_random_draws);
    random_draws = pc * (random_draws');
    ee_f_101(:, ii) = std(random_draws, 0, 2);
end

ee_corr_101 = zeros(size(yy_corr_101));
for ii = 1:size(ee_corr_101, 2)
    random_draws = mvnrnd(yy_corr_3(:, ii)', ee_corr_3{ii}, n_random_draws);
    random_draws = pc * (random_draws');
    ee_corr_101(:, ii) = std(random_draws, 0, 2);
end

% rescale features and errors into picoamps
yy_f_101 = (yy_f_101 + rescale_mean) .* rescale_factor;
ee_f_101 = ee_f_101 .* rescale_factor;

yy_corr_101 = (yy_corr_101 + rescale_mean) .* rescale_factor;
yy_corr_101 = yy_corr_101 + calshift;
ee_corr_101 = ee_corr_101 .* rescale_factor;

% get alignment of f to r
f_to_r = er.tf_to_r(er.f_to_tf);
f_to_r_base = er.tf_to_r(er.f_to_tf);

% figure out smart step counts
[p_smart, pbad_smart] = smartStepCounts(er.features_r, er.Xraw_3_r, true, [er.scale, er.offset], 'svmtrust', 2, 'map', map);
residual_skips = sum(p_smart(:, 1, 2:12), 3);
residual_skips = residual_skips >= smart_skip_threshold;

% build in the inferred skips into the f_to_r
residual_skips = find(residual_skips == 1)';

for ii = numel(residual_skips):-1:1
    f_to_r(f_to_r >= residual_skips(ii)) = f_to_r(f_to_r >= residual_skips(ii)) + 1;
end

% set the x-spacing
xx_f_101 = (1 .* er.details.unfiltered.evaluation_x) + (1:size(yy_f_101, 2));
xx_r_101 = (2 .* er.details.unfiltered.evaluation_x) + f_to_r;
xx_corr_101 = (2 .* er.details.unfiltered.evaluation_x) + (1:size(yy_corr_101, 2));

% figure out transition types
transition_types = diff(f_to_r);

% label which steps are where
m3_steps = zeros(1, length(transition_types)); m3_steps(transition_types <= -3) = 1;
m2_steps = zeros(1, length(transition_types)); m2_steps(transition_types == -2) = 1;
m1_steps = zeros(1, length(transition_types)); m1_steps(transition_types == -1) = 1;
p0_steps = zeros(1, length(transition_types)); p0_steps(transition_types == 0) = 1;
p1_steps = zeros(1, length(transition_types)); p1_steps(transition_types == 1) = 1;
p2_steps = zeros(1, length(transition_types)); p2_steps(transition_types == 2) = 1;
p3_steps = zeros(1, length(transition_types)); p3_steps(transition_types >= 3) = 1;

% extract dc information
yy_f_dc = yy_f_101(90, :);
ee_f_dc = ee_f_101(90, :);

% plotting
ff1 = figure(fignum1); clf; hold on;

% build lines marking transitions
lines_x = [roi(2:end) ; roi(2:end)] + trans_shift;
lines_y = [min([yy_f_dc(roi(2:end)) ; yy_f_dc(roi(1:(end - 1)))]) + vert_shift + under_shift ; max([yy_f_dc(roi(2:end)) ; yy_f_dc(roi(1:(end - 1)))]) + over_shift];

% plot the lines
plot(lines_x(:, m3_steps(roi(1:(end - 1))) == 1), lines_y(:, m3_steps(roi(1:(end - 1))) == 1), 'color', back_colr / line_factor, 'linewidth', trans_lw);
plot(lines_x(:, m2_steps(roi(1:(end - 1))) == 1), lines_y(:, m2_steps(roi(1:(end - 1))) == 1), 'color', back_colr / line_factor, 'linewidth', trans_lw);
plot(lines_x(:, m1_steps(roi(1:(end - 1))) == 1), lines_y(:, m1_steps(roi(1:(end - 1))) == 1), 'color', back_colr / line_factor, 'linewidth', trans_lw);
plot(lines_x(:, p0_steps(roi(1:(end - 1))) == 1), lines_y(:, p0_steps(roi(1:(end - 1))) == 1), 'color', hold_colr / line_factor, 'linewidth', trans_lw);
plot(lines_x(:, p1_steps(roi(1:(end - 1))) == 1), lines_y(:, p1_steps(roi(1:(end - 1))) == 1), 'color', forw_colr / line_factor, 'linewidth', trans_lw);
plot(lines_x(:, p2_steps(roi(1:(end - 1))) == 1), lines_y(:, p2_steps(roi(1:(end - 1))) == 1), 'color', skip_colr / line_factor, 'linewidth', trans_lw);
plot(lines_x(:, p3_steps(roi(1:(end - 1))) == 1), lines_y(:, p3_steps(roi(1:(end - 1))) == 1), 'color', skip_colr / line_factor, 'linewidth', trans_lw);

% put in the markers
for cM = roi(1:(end - 1))
    if m3_steps(cM) == 1
        p = patch(cM + 1 + marker.m3.x + trans_shift, marker_vertloc + marker.m3.y, back_colr);
        p.EdgeColor = back_edge_colr;
        p.LineWidth = edge_width;
    elseif m2_steps(cM) == 1
        p = patch(cM + 1 + marker.m2.x + trans_shift, marker_vertloc + marker.m2.y, back_colr);
        p.EdgeColor = back_edge_colr;
        p.LineWidth = edge_width;
    elseif m1_steps(cM) == 1
        p = patch(cM + 1 + marker.m1.x + trans_shift, marker_vertloc + marker.m1.y, back_colr);
        p.EdgeColor = back_edge_colr;
        p.LineWidth = edge_width;
    elseif p0_steps(cM) == 1
        rectangle('position', [cM + 1 - final_lx_factor * (L / 2) + trans_shift, marker_vertloc - LY_factor_square * (L / 2), final_lx_factor * (L / 3), LY_factor_square * L], 'facecolor', hold_colr, 'edgecolor', [1 1 1]);%'none');
        rectangle('position', [cM + 1 - final_lx_factor * (L / 6) + trans_shift, marker_vertloc - LY_factor_square * (L / 2), final_lx_factor * (L / 3), LY_factor_square * L], 'facecolor', [1 1 1], 'edgecolor', [1 1 1]);%'none');
        rectangle('position', [cM + 1 + final_lx_factor * (L / 6) + trans_shift, marker_vertloc - LY_factor_square * (L / 2), final_lx_factor * (L / 3), LY_factor_square * L], 'facecolor', hold_colr, 'edgecolor', [1 1 1]);%'none');
    elseif p1_steps(cM) == 1
        p = patch(cM + 1 + marker.p1.x + trans_shift, marker_vertloc + marker.p1.y, forw_colr);
        p.EdgeColor = forw_edge_colr;
        p.LineWidth = edge_width;
    elseif p2_steps(cM) == 1
        p = patch(cM + 1 + marker.p2.x + trans_shift, marker_vertloc + marker.p2.y, skip_colr);
        p.EdgeColor = forw_edge_colr;
        p.LineWidth = edge_width;
    elseif p3_steps(cM) == 1
        p = patch(cM + 1 + marker.p3.x + trans_shift, marker_vertloc + marker.p3.y, skip_colr);
        p.EdgeColor = forw_edge_colr;
        p.LineWidth = edge_width;
    end
end

% plotting f_ac
for cL = roi
    shadedErrorBar(xx_f_101(:, cL), yy_f_101(:, cL), ee_f_101(:, cL), {'-', 'color', base_colr, 'linewidth', base_lw}, 1);
end

% plottng just the means again for illustrator
for cL = roi
    plot(xx_f_101(:, cL), yy_f_101(:, cL), '-', 'color', base_colr, 'linewidth', base_lw);
end

% plotting f_dc
for cL = roi
    shadedErrorBar([cL cL + 1], yy_f_dc(cL) .* [1 1] + vert_shift, ee_f_dc(cL) .* [1 1], {'-', 'color', base_colr, 'linewidth', base_lw}, 1);
end

% plot vertical lines on stairs
plot([roi(2:end) ; roi(2:end)], [yy_f_dc(roi(1:(end - 1))) ; yy_f_dc(roi(2:end))] + vert_shift, 'color', base_colr, 'linewidth', base_lw);

% set axis properties
set(gcf, 'renderer', 'painters');

set(gca, 'XLim', [roi(1) + .3, roi(end) + .8]);
set(gca, 'YLim', [.23 * rescale_factor, (.37 * rescale_factor) + vert_shift]);
set(gca, 'YTick', rescale_factor .* [.24, .28, .32, .36, ((.24 * rescale_factor) + vert_shift) / rescale_factor, ((.28 * rescale_factor) + vert_shift) / rescale_factor, ((.32 * rescale_factor) + vert_shift) / rescale_factor, ((.36 * rescale_factor) + vert_shift) / rescale_factor]);
set(gca, 'YTickLabel', [.24 .28 .32 .36 .24 .28 .32 .36]);
set(gca, 'XTick', [(roi(1) + 5):6:(roi(end) + 4)], 'XTickLabel', [6:6:length(roi)]);
set(gca, 'XMinorTick', 'on');
ylabel('Conductance (nS)');
xlabel('Observed Step Number');

set(gca, 'LineWidth', 3, 'FontSize', 24, 'XColor', [0 0 0], 'YColor', [0 0 0]);
ff1.PaperUnits = 'inches';
ff1.PaperPosition = [0 0 8 8];

% error-corrected plot
ff2 = figure(fignum2); clf; hold on;

usedsofar = [];
for cL = roi
    if ~any(usedsofar == f_to_r_base(cL))
        shadedErrorBar(xx_r_101(:, cL), yy_corr_101(:, f_to_r_base(cL)), ee_corr_101(:, f_to_r_base(cL)), {'-', 'color', base_colr, 'linewidth', base_lw}, 1);
        usedsofar = [usedsofar f_to_r_base(cL)];
    else
        continue
    end
end

% for illustrator, plot again (just the means)
for cL = roi
    plot(xx_r_101(:, cL), yy_corr_101(:, f_to_r_base(cL)), '-', 'color', base_colr, 'linewidth', base_lw);
end

% set axis properties
set(gcf, 'renderer', 'painters');

set(gca, 'XLim', [xx_r_101(101, roi(1)),  xx_r_101(1, roi(end))]);
set(gca, 'YLim', [.23 * rescale_factor, .37 * rescale_factor]);
set(gca, 'YTick', rescale_factor .* [.24 .28 .32 .36], 'YTickLabel', [.24 .28 .32 .36]);
set(gca, 'XMinorTick', 'on');
set(gca, 'XTick', [(f_to_r(roi(1)) + 5):6:(f_to_r(roi(end)) + 4)], 'XTickLabel', [6:6:length((f_to_r(roi)))] ./ 2);
ylabel('Conductance (nS)');
xlabel('Corrected DNA Position (nt)');

set(gca, 'LineWidth', 3, 'FontSize', 24, 'XColor', [0 0 0], 'YColor', [0 0 0]);
ff2.PaperUnits = 'inches';
ff2.PaperPosition = [0 0 8 5];

end