function [read] = applyCalibration(read, scale, offset)

% purpose:
%   applies a given scale/offset calibration to a read's _r features and
%   stiffnesses, making it ready to sequence
% inputs:
%   read:
%       a read struct, processed up through reordering filter
%   scale:
%       multiplicative scale calibration to apply
%   offset:
%       additive offset calibration to apply
%

% store calibration used
read.calibration_scale = scale;
read.calibration_offset = offset;

% load principle components
pc = load('principal_components');
pc = pc.principal_components;

% apply calibration to features
n_101 = (pc * read.x3_r) - normalizerAC(pc * read.x3_r);
ncal_101 = (scale .* n_101) + offset;
ncal_3 = pc' * ncal_101;

% apply calibration to stiffnesses
k_3 = read.k3_r;
for cK = 1:numel(k_3); k_3{cK} = k_3{cK} ./ (scale .^ 2); end

% store calibrated results
read.x3_cal = ncal_3;
read.k3_cal = k_3;

end