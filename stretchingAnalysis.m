%% 
% load the consensus
if ~exist('con', 'var')
    load('stretching_sequence_consensus.mat');
    con = consensus; clear consensus
end

% pull out the individually measured, interpolated, calibrated conductance
% values
g_con = con.guf_int_cal;

% extract just the relevant levels from the consensus
g_con = g_con(:, 45:231);
% define our region of interest
roi = 5:2:175;
% define our voltage bins
bins = 1:1:101;

% calculate the all-level averages
g_all = zeros(101, size(g_con, 2));
for cL = 1:size(g_con, 2)
    if size(g_con{cL}, 2) == 0
        fprintf('Warning, %d is empty!\n', cL)
        continue
    end
    g_all(:, cL) = sum(g_con{cL}, 2) ./ size(g_con{cL}, 2);
end

%%
% calculate the nt stretch for all odd levels
fprintf('Stretch All\n')
analysis.odd = zeros(1, length(bins) - 1);
analysis.scale = zeros(1, length(bins) - 1);
analysis.offset = zeros(1, length(bins) - 1);
for ii = 1:length(analysis.odd)
    if mod(ii, 10) == 0
        fprintf('\tBin %d/%d\n', ii, length(analysis.odd))
    end
    [analysis.odd(ii), cal_calc] = calculateNtShift(g_all(:, roi), bins(ii), bins(ii + 1));
    analysis.scale(ii) = cal_calc(1);
    analysis.offset(ii) = cal_calc(2);
end

%% 
% bootstraping by subsections
boot_levels = 16;
boot_overlap = 8;

% initialize storage
analysis.boots = cell(1, 0);

% do the bootstrap
fprintf('Bootstrapping\n')
exitflag = false;
iter = 0;
while ~exitflag
    iter = iter + 1;
    roi_boot_ix = (1:boot_levels) + ((iter - 1) * boot_overlap);
    if max(roi_boot_ix) > length(roi)
        exitflag = true;
        continue
    end
    roi_boot = roi(roi_boot_ix);
    analysis.boots{end + 1} = zeros(1, length(bins) - 1);
    fprintf('\tBoot %d\n', iter)
    for ii = 1:length(analysis.boots{end})
        if mod(ii, 10) == 0
            fprintf('\t\tBin %d/%d\n', ii, size(analysis.boots{end}, 2))
        end
        analysis.boots{end}(ii) = calculateNtShift(g_all(:, roi_boot), bins(ii), bins(ii + 1), [analysis.scale(ii), analysis.offset(ii)]);
    end
end

analysis.enzyme = 'Hel308';
analysis.roi = roi;
analysis.bootlevels = 16;
analysis.bootoverlap = 8;

%%
% save the product
save('stretching_analysis_result.mat', 'analysis', '-v7.3');



