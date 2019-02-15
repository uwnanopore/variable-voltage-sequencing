% general stuff
% define alphabet
alphabet = 'ACGT';
roundlen = 3;
fontsize = 24;
fontcolor = [0 0 0];
conpos = [0 0 7 7];
bins = [0.49:0.03:1.05];
randfacecolor = [.7 .1 .1];
accfacecolor = [.1 .1 .7];
normalization = 'count';
minacc = .54;
yshift = 1;
histpos = [0 0 16 5];

%% AC
% load result
load('variable_voltage_sequencing_results.mat');
ac = variable_voltage_sequencing_results; clear variable_voltage_sequencing_results
accs = ac.accuracy;
als = ac.alignment;
bigal = ac.big_alignment;
acrandaccs = ac.random_accuracy;

% initialize confusion matrix
accon = zeros(5);

% make bigal numeric
nbigal = zeros(size(bigal));
nbigal(bigal == 'A') = 1;
nbigal(bigal == 'C') = 2;
nbigal(bigal == 'G') = 3;
nbigal(bigal == 'T') = 4;
nbigal(bigal == ' ') = 5;

% fillup the confusion matrix
for cB = 1:size(bigal, 2)
    accon(nbigal(1, cB), nbigal(2, cB)) = accon(nbigal(1, cB), nbigal(2, cB)) + 1;
end

% initialize the normalized confusion matrix
acncon = zeros(5);

% fill the normalized matrix
acncon(1:4, :) = accon(1:4, :) ./ sum(accon(1:4, :), 2);
acncon(5, :) = accon(5, :) ./ sum(accon, 1);

% round off acncon
acncon = round(acncon, roundlen);
acncon = 100 .* acncon;

% plotting
fac = figure(23902); clf; %hold on;
hac = heatmap({'A', 'C', 'G', 'T', 'del'}, {'A', 'C', 'G', 'T', 'ins'}, acncon, 'ColorLimits', [0 100]);

% plot options
set(gcf, 'renderer', 'painters');
hac.FontSize = fontsize;
hac.FontColor = fontcolor;
fac.PaperUnits = 'inches';
fac.PaperPosition = conpos;

% read accuracy histogram
fach = figure(30942); clf; hold on;
hrandac = histogram(acrandaccs, bins, 'normalization', normalization, 'facecolor', randfacecolor);
haccac = histogram(ac.accuracy, bins, 'normalization', normalization, 'facecolor', accfacecolor);
xlim([minacc 1]);
maxcount = max([hrandac.Values, haccac.Values]);
ylim([0 maxcount + yshift]);

set(gcf, 'renderer', 'painters');
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 2.5);
set(gca, 'XTick', [.55:.05:1], 'XTickLabel', [.55:.05:1] .* 100);
set(gca, 'YTick', [10 30 50]);
set(gca, 'FontSize', fontsize);
xlabel('Accuracy (%)');
ylabel('Counts');
fach.PaperUnits = 'inches';
fach.PaperPosition = histpos;

%% DC
% load result
load('constant_voltage_sequencing_results.mat');
dc = constant_voltage_sequencing_results; clear constant_voltage_sequencing_results;
accs = dc.accuracy;
als = dc.alignment;
bigal = cell2mat(als);
dcrandaccs = dc.random_accuracy;

% initialize confusion matrix
dccon = zeros(5);

% make bigal numeric
nbigal = zeros(size(bigal));
nbigal(bigal == 'A') = 1;
nbigal(bigal == 'C') = 2;
nbigal(bigal == 'G') = 3;
nbigal(bigal == 'T') = 4;
nbigal(bigal == ' ') = 5;

% fillup the confusion matrix
for cB = 1:size(bigal, 2)
    dccon(nbigal(1, cB), nbigal(2, cB)) = dccon(nbigal(1, cB), nbigal(2, cB)) + 1;
end

% initialize the normalized confusion matrix
dcncon = zeros(5);

% fill the normalized matrix
dcncon(1:4, :) = dccon(1:4, :) ./ sum(dccon(1:4, :), 2);
dcncon(5, :) = dccon(5, :) ./ sum(dccon, 1);

% round off acncon
dcncon = round(dcncon, roundlen);
dcncon = 100 .* dcncon;

% plotting
fdc = figure(53345); clf; %hold on;
hdc = heatmap({'A', 'C', 'G', 'T', 'del'}, {'A', 'C', 'G', 'T', 'ins'}, dcncon, 'ColorLimits', [0 100]);

% plot options
set(gcf, 'renderer', 'painters');
hdc.FontSize = fontsize;
hdc.FontColor = fontcolor;
fdc.PaperUnits = 'inches';
fdc.PaperPosition = conpos;

% read accuracy histogram
fdch = figure(66732); clf; hold on;
hranddc = histogram(dcrandaccs, bins, 'normalization', normalization, 'facecolor', randfacecolor);
haccdc = histogram(dc.accuracy, bins, 'normalization', normalization, 'facecolor', accfacecolor);
xlim([minacc 1]);
maxcount = max([hranddc.Values, haccdc.Values]);
ylim([0 maxcount + yshift]);

set(gcf, 'renderer', 'painters');
set(gca, 'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 2.5);
set(gca, 'XTick', [.55:.05:1], 'XTickLabel', [.55:.05:1] .* 100);
% set(gca, 'YTick', [10 30 50]);
set(gca, 'YTick', [5 10 15]);
set(gca, 'FontSize', fontsize);
xlabel('Accuracy (%)');
ylabel('Counts');
fdch.PaperUnits = 'inches';
fdch.PaperPosition = histpos;

