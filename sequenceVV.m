function [sequence] = sequenceVV(read, varargin)

% purpose:
%   function serves as a wrapper, taking a processed variable voltage read
%   (all filtering steps conducted, and calibration applied) and determines
%   the maximum likelihood sequence to have generated it, subject to
%   allowed state-to-state transitions and a pore model of
%   sequence-to-signal.
% inputs:
%   read:
%       a fully processed variable voltage read, as found in
%       /variableVoltageReads/read_*_section_*.mat. must contain the fields
%       x101_r, x3_r, x3_cal, k3_cal, calibration_scale, and
%       calibration_offset
% outputs:
%   sequence:
%       maximum likelihood sequence generating the observed variable
%       voltage levels, given in 5' -> 3' order
% varargin:
%   map:
%       user input of the pore model to use in decoding and in
%       smartStepCounts
%   transitioninfo:
%       user input of the allowed state-to-state transitions and their
%       associated penalties for decoding
%

% handle varargin
for cV = 1:length(varargin)
    if ~ischar(varargin{cV})
        continue
    end
    switch lower(varargin{cV})
        case 'map'
            map = varargin{cV + 1};
        case 'transitioninfo'
            transition_info = varargin{cV + 1};
    end
end

% load structs needed for step counts and sequencing
if ~exist('map', 'var')
    map = load('pore_model_6mer_variable_voltage.mat');
    map = map.model;
end
if ~exist('transition_info', 'var')
    transition_info = load('transition_info_hel308_6mer.mat');
    transition_info = transition_info.transition_info;
end

% calculate smart step counts
x101_for_stepcounts = read.x101_r;
x3_for_stepcounts = read.x3_r;
scale_for_stepcounts = read.calibration_scale;
offset_for_stepcounts = read.calibration_offset;
[p_steps, p_bads] = smartStepCounts( ...
    x101_for_stepcounts, ...
    x3_for_stepcounts, ...
    true, ...
    [scale_for_stepcounts, offset_for_stepcounts], ...
    'map', map);

% sequencing
sequence = calculateSequenceVV( ...
    read.x3_cal, ...
    read.k3_cal, ...
    p_steps, ...
    p_bads, ...
    'map', map, ...
    'backsteps', read.did_backstep_r', ...
    'scorecutoff', -10, ...
    'transitioninfo', transition_info, ...
    'quiet');

end