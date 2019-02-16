function [transitions, features, errors, stiffnesses] = findLevels_CVsub(data, varargin)
% recursive maximum-likelihood level finder
% includes change point information criterion bias to accurately find
% levels on all time scales.
%
%
% INPUTS
%
% data              a time series of data we want to find transitions in
%
% 
%
% OPTIONS
%
% sensitivity       the multiplier for the CPIC bias. higher numbers lead
%                   to a less sensitive level finder. natural value is 1.
%                   numbers less than 0 don't make sense.
%
% minlevellength    the shortest level length we will search for. needs to
%                   be at least 2.
%

% 
% OUTPUTS
% 
% transitions       the indices of transition points
%
% features          a 2x(num levels) array of 
%                   [medians; standard deviations] for each level
%
% errors            a 2x(num levels) array of the uncertainties in the 
%                   features
%
% stiffnesses       a 1x(num levels) cell array of the fisher information
%                   matrices for each level


original_data_mapping = 1:length(data);
data_length = length(data);
original_data_mapping(isnan(data)) = [];
data(isnan(data)) = [];

%declare global variables
global x xsq low_n_cpic high_n_cpic min_level_length cpic_multiplier

%the function defining the cpic (bias estimator for introducing
%transitions)
load('cpic_fits_CV.mat', 'cpic_fits');
low_n_cpic = cpic_fits.lessthan1000;
high_n_cpic = cpic_fits.greaterthan1000;

%define the number determining level finder sensitivity
cpic_multiplier = 9;

%the shortest level length allowed.
min_level_length = 2;

%Take optional arguments
for ca = 1:length(varargin)
    switch upper(char(varargin{ca}))
        case 'SENSITIVITY'
            cpic_multiplier = varargin{ca+1};
        case 'MINLEVELLENGTH'
            min_level_length = varargin{ca+1};
    end
end

%cumulative sum of data and squared data, used to rapidly calculate the
%sum-squared-error in different regions
x = [0 cumsum(data)];
xsq = [0 cumsum(data.^2)];

%call the recursive level finding function

transitions = sort(findTransitions(2, length(data)))-1;

%delete transitions at the beginning and end which are too short
transitions(transitions <= min_level_length ...
                    | transitions > length(data) - min_level_length) = [];

%Append 0 and the end
transitions = [0 transitions length(data)];


%preallocate the arrays containing calculated data about each level
features = zeros(2, length(transitions)-1);
errors = zeros(2, length(transitions)-1);
stiffnesses = cell(1, length(transitions)-1);

%loop through all the levels
for ct = 2:length(transitions)
    
    %features are median and noise
    features(:,ct-1) = [median( data( (transitions(ct-1)+1):transitions(ct) ) );
                           std( data( (transitions(ct-1)+1):transitions(ct) ) )];
    
    %standard errors in the features
    errors(:,ct-1) = [features(2,ct-1)/sqrt(transitions(ct)-transitions(ct-1)-1);
                      features(2,ct-1)/sqrt( 2*(transitions(ct)-transitions(ct-1)-1))];
    
    %stiffness is estimated as diagonal
    stiffnesses{ct-1} = diag(errors(:,ct-1).^-2);
    
end


transitions = [0 original_data_mapping(transitions(2:end-1)) data_length];

end




function transitions = findTransitions(left, right)
%Recursive function that calls itself to continually divide and examine
%regions of the data.

%declare global variables
global x xsq low_n_cpic high_n_cpic min_level_length cpic_multiplier

%by default, we have no transitions. the empty array will be returned if
%nothing is found.
transitions = [];

%the number of points total, to the left of the transition point, and to
%the right of the transition point. N_L and N_R are vectors, since we
%simultaneously calculate the likelihood for all candidate points.
N_T = right-left+1;
N_L = min_level_length:N_T-min_level_length;
N_R = N_T - N_L;

%if the even is under twice the minimum level length, there are no
%candidate points and we end by returning an empty array
if isempty(N_L)
    return
end

%calculate the mean of the left, the right, and the entire interval. THE
%TRANSITION POINT IS THE LAST POINT OF THE LEFT PARTITION.
x_mean_L = (x(left + N_L - 1) - x(left-1))./N_L;
x_mean_R = (x(right) - x(right - N_R))./N_R;
x_mean_T = (x(right) - x(left-1))/N_T;

%calculate the mean of the square of the left, the right, and the entire
%interval.
xsq_mean_L = (xsq(left + N_L - 1) - xsq(left-1))./N_L;
xsq_mean_R = (xsq(right) - xsq(right - N_R))./N_R;
xsq_mean_T = (xsq(right) - xsq(left-1))/N_T;

%calculate the variance <dx^2> = <x^2> - <x>^2
var_L = max(xsq_mean_L - x_mean_L.^2, 0.0003);
var_R = max(xsq_mean_R - x_mean_R.^2, 0.0003);
var_T = max(xsq_mean_T - x_mean_T.^2, 0.0003);

%calculate the CPIC penalty from the loaded functions. This penalty was
%calculated by Monte Carlo a la LaMont/Wiggins 2016, and serves to prevent
%the overfitting bias in the likelihood maximizing model.'
if N_T >= 1e6
    p_CPIC = high_n_cpic(1e6);
elseif N_T > 1000 && N_T < 1e6
    p_CPIC = high_n_cpic(N_T);
else
    p_CPIC = low_n_cpic(N_T);
end

%The total CPIC is the negative log likelihood plus the CPIC penalty times
%a multiplier that determines the algorithm's sensitivity. A higher
%multiplier means a less-sensitive level finder. Totally naïve level
%finding should see cpic_multiplier = 1.
CPIC = 0.5*(N_L.*reallog(var_L) + N_R.*reallog(var_R) - N_T*reallog(var_T)) ...
    + cpic_multiplier*p_CPIC;

%find the best possible transition point
[minCPIC, wheremin] = min(CPIC);

%correct the index from the location in the edge-chopped interval to the
%location in the arrays x and xsq.
min_index = wheremin+min_level_length+left-2;

%clear dangling vectors that shouldn't hang around as we descend into
%recursion
clear N_T N_L N_R x_mean_L x_mean_R x_mean_T ...
    xsq_mean_L xsq_mean_R xsq_mean_T var_L var_R var_T CPIC p_CPIC

%CPIC < 0 means that encoding is better if we make the partition
if minCPIC < 0
    
    %return the found transition, together with any found by repeating the
    %search on the intervals left and right of the transition.
    transitions = [ min_index ...
        findTransitions(left, min_index) ...
        findTransitions(min_index+1, right)];
    
end


end
