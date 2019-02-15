function [transitions] = findLevelsPC(data, vdata, T, basis_fns, varargin)
%arbitrary basis function level finder

%declare global variables
global XB Xsq BB low_n_cpic high_n_cpic min_level_length cpic_multiplier num_bf period

%the function defining the cpic (bias estimator for introducing
%transitions)
cpic_fits = load('cpic_fits.mat');
cpic_fits = cpic_fits.cpic_fits;

low_n_cpic = cpic_fits.lessthan1000;
high_n_cpic = cpic_fits.greaterthan1000;

%define the number determining level finder sensitivity
cpic_multiplier = 4;

%the shortest level length allowed.
min_level_length = T;

%the number of basis functions
num_bf = size(basis_fns, 2);


%Take optional arguments
for ca = 1:length(varargin)
    switch upper(char(varargin{ca}))
        case 'SENSITIVITY'
            cpic_multiplier = varargin{ca+1};
        case 'MINLEVELLENGTH'
            if varargin{ca+1} < num_bf + 1
                error(['cannot set min level length below number of model '...
                    'parameters = ' num2str(num_bf + 1) '.'])
            end
            min_level_length = varargin{ca+1};
    end
end

% %find the phase
fourier_component = sum( exp(1i*2*pi/T*(1:min(1e4,length(vdata))))...
                                            .*vdata(1:min(1e4,length(vdata))));
if abs(fourier_component) < 1
    error(['no AC voltage applied in first 10,000 points of ' ...
                                'file- could not determine phase.'])
end
phase = round((angle(fourier_component)+pi)/(2*pi/T));

%delete the spikes
pts_to_delete = repmat(((phase-T):(T/2):length(data)+T/2)', 1, 21);
pts_to_delete = pts_to_delete + repmat(0:20, size(pts_to_delete, 1), 1);
pts_to_delete = reshape(pts_to_delete, 1, []);
pts_to_delete(pts_to_delete < 1 | pts_to_delete > length(data)) = [];
original_indices = 1:length(data);
L_original = length(data);
data(pts_to_delete) = [];
original_indices(pts_to_delete) = [];
T = T - 2*21;
period = T;
phase = phase-sum(pts_to_delete < phase);


%call the recursive level finding function

transition_list = [];
for seg = 0:1e6:length(data)
    start_point = seg+1;
    end_point = min(length(data),seg + 1.1e6);
    data_segment = data(start_point:end_point);
    
    x = [0 data_segment];
    t = start_point-1:end_point;
    
    %define basis functions
    
    phased_bf = repmat(circshift(basis_fns', mod(phase - seg, size(basis_fns, 1)), 2), 1, ceil(length(t)/size(basis_fns, 1)));
    bf = phased_bf(:,1:length(t));
    
    
    %define necessary cumulates
    Xsq = cumsum(x.^2); %cumulate of x^2
    XB = cumsum( repmat(x, size(bf, 1), 1).*bf , 2); %vector of cumulates of x*b_i
    BB = zeros(num_bf, length(x)); %matrix of cumulates of b_i*b_j, the cumulate ...
                                            %for each entry of the matrix is a row of BB
    for cbf = 1:num_bf
        BB( ((cbf-1)*num_bf+1):(cbf*num_bf), :) = cumsum(repmat(bf(cbf,:), ...
                                                            size(bf, 1), 1).*bf, 2);
    end
    transition_list = [transition_list, findTransitions(2, length(x)) + start_point - 1];
end

transitions = sort(unique(transition_list))-1;

%Append 0 and the end
transitions = [0 original_indices(transitions) L_original];

end




function transitions = findTransitions(left, right)
%Recursive function that calls itself to continually divide and examine
%regions of the data.

%declare global variables
global Xsq XB BB low_n_cpic high_n_cpic min_level_length cpic_multiplier num_bf period

%by default, we have no transitions. the empty array will be returned if
%nothing is found.
transitions = [];

%the number of points total, to the left of the transition point, and to
%the right of the transition point. N_L and N_R are vectors, since we
%simultaneously calculate the likelihood for all candidate points.
N_T =right-left+1;
N_L = min_level_length:N_T-min_level_length;
N_R = fliplr(N_L);

%if the even is under twice the minimum level length, there are no
%candidate points and we end by returning an empty array
if isempty(N_L)
    return
end

%calculate the mean of the left, the right, and the entire interval. THE
%TRANSITION POINT IS THE LAST POINT OF THE LEFT PARTITION.

B_matrix_L = BB(:,left+N_L-1)            - BB(:,left-1)*ones(size(N_L));
B_matrix_R = BB(:,right)*ones(size(N_R)) - BB(:,right-N_R);
B_matrix_T = BB(:,right)                 - BB(:,left-1);

XB_L =       XB(:,left+N_L-1)            - XB(:,left-1)*ones(size(N_L));
XB_R =       XB(:,right)*ones(size(N_R)) - XB(:,right-N_R);
XB_T =       XB(:,right)                 - XB(:,left-1);

Xsq_L =      Xsq(left+N_L-1)             - Xsq(left-1);
Xsq_R =      Xsq(right)                  - Xsq(right-N_R);
Xsq_T =      Xsq(right)                  - Xsq(left-1);


var_L = zeros(size(N_L));
var_R = zeros(size(N_R));


for ct = 1:period:length(N_L)
    var_L(ct) = (Xsq_L(ct) - ( XB_L(:, ct)'/reshape(B_matrix_L(:,ct), num_bf, num_bf) )...
                                                                   *XB_L(:, ct) )/N_L(ct);
    var_R(ct) = (Xsq_R(ct) - ( XB_R(:, ct)'/reshape(B_matrix_R(:,ct), num_bf, num_bf) )...
                                                                   *XB_R(:, ct) )/N_R(ct);
end
var_T = (Xsq_T - (XB_T'/reshape(B_matrix_T, num_bf, num_bf) )*XB_T)/N_T;

var_L(var_L <= 0) = nan;
var_R(var_R <= 0) = nan;

%calculate the CPIC penalty from the loaded functions. This penalty was
%calculated by Monte Carlo a la LaMont/Wiggins 2016, and serves to prevent
%the overfitting bias in the likelihood maximizing model.
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
CPIC = 0.5*(N_L.*log(var_L) + N_R.*log(var_R) - N_T*log(var_T)) ...
    + cpic_multiplier*p_CPIC;

%find the best possible transition point
[~, wheremin] = nanmin(CPIC);


%calculate the summed variance from the model of each side
for ct = max(wheremin-period,1):min(wheremin+period,length(N_L))
    var_L(ct) = (Xsq_L(ct) - ( XB_L(:, ct)'/reshape(B_matrix_L(:,ct), num_bf, num_bf) )...
                                                                   *XB_L(:, ct) )/N_L(ct);
    var_R(ct) = (Xsq_R(ct) - ( XB_R(:, ct)'/reshape(B_matrix_R(:,ct), num_bf, num_bf) )...
                                                                   *XB_R(:, ct) )/N_R(ct);
end
var_L(var_L <= 0) = nan;
var_R(var_R <= 0) = nan;

%calculate the CPIC
CPIC = 0.5*(N_L.*log(var_L) + N_R.*log(var_R) - N_T*log(var_T)) ...
    + cpic_multiplier*p_CPIC;
[minCPIC, wheremin] = nanmin(CPIC);

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