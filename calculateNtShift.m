function [nt_shift, calibration] = calculateNtShift(G, vix1, vix2, calibration)

% default values
if ~exist('increment', 'var')
    increment = 0.002;
end

% calculate the x values at which to calculate the spline
points_per_level = 1 / increment;
xx = linspace(1, size(G, 2), (size(G, 2) - 1) * points_per_level + 1);

% calculate the spline
y1 = G(vix1, :);
y2 = G(vix2, :);

yy1 = spline(1:size(G, 2), y1, xx);
yy2 = spline(1:size(G, 2), y2, xx);

% define which shifts we are goint to test
max_left_shift = -0.5;
max_right_shift = 0.5;

shift_points = (max_left_shift * 2 * points_per_level):1:(max_right_shift * 2 * points_per_level);

% calculate the calibration between yy1 and yy2
if ~exist('calibration', 'var')
    [scale, offset] = calibrate_by_iqr(yy1, yy2, false, [-inf inf]);
else
    scale = calibration(1);
    offset = calibration(2);
end

% apply calibration
yy1 = (scale .* yy1) + offset;
calibration = [scale, offset];

% calculate scores at each shift position
scores = nan(1, length(shift_points));
for cS = 1:length(scores)
    scores(cS) = calculateResidual(yy1, yy2, shift_points(cS));
end

% find where the best score occured
[~, optimal_ix] = min(scores);

% convert the optimal point to a shift measurement
nt_shift = shift_points(optimal_ix) / (points_per_level);

end

% ----------------------------------------------- %
function score = calculateResidual(yy1, yy2, shift)

% figure out the points to take from yy1 and yy2
if shift < 0
    left1 = 1;
    right1 = length(yy1) + shift;
    ix1 = left1:right1;
    
    left2 = 1 - shift;
    right2 = length(yy2);
    ix2 = left2:right2;
    
elseif shift == 0
    ix1 = 1:length(yy1);
    ix2 = 1:length(yy2);
    
elseif shift > 0
    left1 = 1 + shift;
    right1 = length(yy1);
    ix1 = left1:right1;
    
    left2 = 1;
    right2 = length(yy2) - shift;
    ix2 = left2:right2;
end

% grab the points from the spline
yy1 = yy1(ix1);
yy2 = yy2(ix2);

% calculate the sum square error
sse = sum((yy2 - yy1) .^ 2);

% weight the sse by the number of shared points
score = sse / numel(yy1);

end