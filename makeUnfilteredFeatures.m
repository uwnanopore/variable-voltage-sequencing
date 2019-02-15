function [features_uf, stiffnesses_uf, npts_uf, evaluation_voltage, evaluation_x, feature_type] = makeUnfilteredFeatures(levels)

% describe feature type
feature_type = 'Features are conductances, evaluated at equal x spacing from x(110mV) to x(190mV) with a total of 101 points via linear interpolation to IV averages. Equal spacing from equal_x_spacing.mat. Openstate calibration not applied. Errors obtained from standard error of mean of binned values';

% load in evaluation voltages and positions for equally spaced points in
% position
spacing = load('equal_x_spacing.mat');
spacing = spacing.spacing;
evaluation_voltage = spacing.v';
evaluation_x = spacing.x';

% set number of features and levels in the event, initialize features,
% stiffnesses, and number of points
n_features = 101;
n_levels = length(levels);

features_uf = nan(n_features, n_levels);
stiffnesses_uf = cell(1, n_levels);
for cL = 1:n_levels
    stiffnesses_uf{cL} = nan(n_features, n_features);
end
npts_uf = nan(1, n_levels);

% loop through all levels in event
for cL = 1:n_levels
    % pull out avg and sem current and voltage from cap comp results
    i_avg = levels{cL}.avgI;
    i_err = levels{cL}.semI;
    v_avg = levels{cL}.avgV;
    
    % will use linear interpolation to evaluate at feature voltages
    % (feature positions). first need to put information into column vector
    % format and avoid nans
    temp = v_avg;
    [v_fit_avg, i_fit_avg] = fitPrep2(temp, i_avg);
    [v_fit_err, i_fit_err] = fitPrep2(temp, i_err);
    clear temp    
    
    % evaluate at all evaluation voltages using linear interpolation
    if sum(~isnan(i_fit_avg) & ~isnan(v_fit_avg)) > 10
        i_avg = interp1q(v_fit_avg, i_fit_avg, evaluation_voltage);
    else
        i_avg = nan .* evaluation_voltage;
        fprintf('WARNING: level %d/%d has too many NaN errors\n', cL, n_levels)
    end
    
    % only do for i_err if > 10 points with real error values, otherwise
    % fill with NaNs and warn
    if sum(~isnan(i_fit_err) & ~isnan(v_fit_err)) > 10
        i_err = interp1q(v_fit_err, i_fit_err, evaluation_voltage);
    else
        i_err = nan .* evaluation_voltage;
        fprintf('WARNING: level %d/%d has too many NaN errors\n', cL, n_levels)
    end

    % calculate conductance (G = I / V) and fill in features
    g_avg = i_avg ./ evaluation_voltage;
    g_err = i_err ./ evaluation_voltage;
    features_uf(:, cL) = g_avg';
    
    % calculate error floor by fitting quadratic curve to v vs g data
    [error_floor_x, error_floor_y] = fitPrep2(levels{cL}.avgV, levels{cL}.avgI ./ levels{cL}.avgV);
    if sum(~isnan(error_floor_x) & ~isnan(error_floor_y)) > 10
        [error_floor_params, ~, error_floor_entropy] = quadraticFitWithFisherInfo(error_floor_x, error_floor_y);
    else
        error_floor_params = nan;
        error_floor_entropy = nan;
    end
    
    if sum(isnan(error_floor_params)) < 1
        error_floor = prctile(abs(polyval(error_floor_params(1:3), error_floor_x) - error_floor_y), 68);
    else
        error_floor = nan;
    end
    
    % apply error floor to conductance errors
    g_err((g_err < error_floor) | (isnan(g_err))) = error_floor;
    
    % use conductance errors to fill in stiffnesses
    stiffnesses_uf{cL} = diag(g_err .^ -2);
    
    % fill in npts_uf
    npts_uf(cL) = sum(levels{cL}.numPts);
end

end


