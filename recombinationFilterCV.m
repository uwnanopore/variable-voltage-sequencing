% removeToggles() performs a self-alignment to combine repeated levels,
% saving the most reliable (longest-duration) level in each combined level
%
% removeToggles(datafolders) -- adds toggle-filtered (TF) fields to the
% jumps data structures
%
% removeToggles( [ medians, stdevs ], errs, durations ) -- returns the
% self-alignment and consensus levels (selected based on durations)

function filteredread = recombinationFilterCV(read, varargin)


error_floor = [0 100];%0.0005;
psa = -4; %XXX
lookback = 10;
% step_probabilities = [step | skip | skip+ | back | back+ | hold | bad | bad+ | optional p_sa/p_slip]
step_probabilities =  [0.6 1e-2 1e-2 .27 .2 .13 0 0];

for ca = 1:length(varargin)
    if ~ischar(varargin{ca})
        continue
    end
    switch upper(char(varargin{ca}))
        case 'PSA'
            psa = varargin{ca+1};
        case 'LOOKBACK'
            lookback = varargin{ca+1};
        case 'STEPPROBABILITIES'
            step_probabilities = varargin{ca+1};
        case 'ERRORFLOOR'
            error_floor = varargin{ca+1};
    end
end

cov_floor = diag(error_floor.^2);






% do the actual self-alignment

filteredFeatures = read.x_f;
filteredStiffs = read.k_f;
filteredNumPts = read.npts_f;
total_self_alignment = 1:size(filteredFeatures,2);


selfalignment = ones(1, size(filteredFeatures,2));


while length(selfalignment) ~= max(selfalignment)
    
    
    
    stiffs_array = cell(size(filteredStiffs));
    for cL = 1:length(filteredStiffs)
        covmat = inv(filteredStiffs{cL});
        covmat = covmat + cov_floor;
        stiffs_array{cL} = inv(covmat);
    end
    
    
    selfalignment = alignLevels(filteredFeatures, stiffs_array, [], [], 'stepprobabilities', [step_probabilities nan], 'selfalignment', psa, lookback);
    
    newFilteredFeatures = zeros(size(filteredFeatures, 1), max(selfalignment));
    newFilteredStiffs = cell(1, max(selfalignment));
    newFilteredNumPts = zeros(1, max(selfalignment));
    
    for cL = 1:size(newFilteredFeatures, 2)
        this_level_locations = find(selfalignment == cL);
        [newFilteredFeatures(:, cL), total_stiff] = multivarSampleStats(filteredFeatures(:, this_level_locations), filteredStiffs(this_level_locations));
        newFilteredStiffs{cL} = total_stiff;
        newFilteredNumPts(cL) = sum(filteredNumPts(this_level_locations));
    end
    
    for csl = 1:length(selfalignment)
        total_self_alignment(total_self_alignment == csl) = selfalignment(csl);
    end
    
    filteredFeatures = newFilteredFeatures;
    filteredStiffs = newFilteredStiffs;
    filteredNumPts = newFilteredNumPts;
end

read.f_to_tf = total_self_alignment;
read.x_tf = filteredFeatures;
read.k_tf = filteredStiffs;

filteredread = read;
