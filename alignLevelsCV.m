function [alignment, best_score, score_matrix, alignment_matrix, ...
    cumulate_score] = alignLevelsCV(x_meas, K_meas, x_ref, K_ref, varargin)
%self alignment:
% [alignment, aligned_levels, best_score, score_matrix, alignment_matrix, cumulate_score] = alignLevels([1 2 3 4 5 4 5 6 7 8], num2cell(0.01^-2*ones(1, 10), [], [], 'stepcounts', [1000 100 1 200 1 200 0 0 350], 'selfalignment')
%reference alignment:
% [alignment, aligned_levels, best_score, score_matrix, alignment_matrix, cumulate_score] = alignLevels([1 2 3 4 5 4 5 6 7 8], num2cell(0.01^-2*ones(1, 10), [1 2 3 4 5 6 7 8], num2cell(0.01^-2*ones(1, 8), 'stepcounts', [1000 100 1 200 1 200 0 0 0])


step_penalties = [-0.2528 ...  
                  -2.0928 ...
                  -2.0928 ...
                  -3.8128 ...
                  -5.9128 ...
                  -2.9128 ...
                  -5.1328 ...
                  -1.8496 ...
                  -8.5100
                  ];
              
pfbad = -1000000*ones(1, size(x_meas, 1));
use_periodic_boundaries = false;
modes = 1;
lookback = 10;
slip_locations = zeros(1, size(x_meas, 2));
sa_lookback = -1; %-1 for no self alignment
alignment_type = 'R';
%psa = -50; %reasonable estimate to push self-alignments lower 
           %than well-matched backsteps or holds 
psa = -14; %extreme estimate to push self-alignments lower than all
%             %backsteps and holds
reorder = false;
isdiag = false;
isdc = false;
isprincomp = false;
for ca = 1:length(varargin)
    if ischar(varargin{ca})
    switch upper(varargin{ca})
        case 'STEPCOUNTS'
            step_counts = varargin{ca+1};
        case 'STEPPROBABILITIES'
            step_probabilities = varargin{ca+1};
        case 'STEPPENALTIES'
            step_penalties = varargin{ca+1};
        case 'MODES'
            modes = varargin{ca+1};
        case 'PERIODIC'
            use_periodic_boundaries = true;
        case 'SCOREMATRIX'
            score_matrix = varargin{ca+1};
        case 'LOOKBACK'
            lookback = varargin{ca+1};
        case 'SLIPLOCS'
            slip_locations(varargin{ca+1}) = 1;
        case 'SELFALIGNMENT'
            sa_lookback = min(7, size(x_meas, 2)-1);%length(x_meas)-1;
            alignment_type = 'S';
            x_ref = x_meas;
            K_ref = K_meas;
            if isnumeric(varargin{ca+1})
                psa = varargin{ca+1};
            end
            if isnumeric(varargin{ca+1}) && isnumeric(varargin{ca+2})
                sa_lookback = min(varargin{ca+2}, size(x_meas,2)-1);
                lookback = sa_lookback;
            end

            step_penalties(:,end) = psa;

        case 'PFBAD'
            pfbad = varargin{ca+1};
        case 'SELFALIGNMENTLOOKBACK'
            sa_lookback = varargin{ca+1};
        case 'REORDER'
            reorder = true;
        case 'DIAG'
            isdiag = true;
        case 'DC'
            isdc = true;
        case 'PRINCOMP'
            isprincomp = true;
    end
    end
end


if isempty(x_meas) || isempty(x_ref)
    alignment = [];
    best_score = 0;
    score_matrix = [];
    alignment_matrix = [];
    cumulate_score = [];
    return;
end

if ~exist('score_matrix', 'var')
    if isdiag
        score_matrix = smatrixACdiag(x_meas, K_meas, x_ref, K_ref);
    elseif isdc
        score_matrix = smatrix1d(x_meas, K_meas, x_ref, K_ref);
    elseif isprincomp
        score_matrix = smatrixACpar(x_meas, K_meas, x_ref, K_ref);
    else
        score_matrix = smatrixACBF(x_meas, K_meas, x_ref, K_ref, 'selfalignment', sa_lookback, 'pfbad', pfbad(1));
    end
end

if ~exist('step_probabilities', 'var') && exist('step_counts', 'var')
    step_probabilities = 0*step_counts;
    for cm = 1:length(modes)
        step_probabilities(cm,[1 2 4 6 7 9]) = step_counts(cm,[1 2 4 6 7 9])/sum(step_counts(cm,[1 2 4 6 7 9]));
        step_probabilities(cm, 3) = step_counts(cm, 3)/step_counts(cm, 2);
        step_probabilities(cm, 5) = step_counts(cm, 5)/step_counts(cm, 4);
        step_probabilities(cm, 8) = step_counts(cm, 8)/step_counts(cm, 7);
    end

end

if exist('step_probabilities', 'var')
    step_penalties = log(step_probabilities);
    step_penalties(:,8) = step_penalties(:,7);
    if alignment_type == 'S'
        step_penalties(:,9) = psa;
    end
end

if ismac && ~exist('alignmentC.mexmaci64', 'file')
    mex alignmentC.c
end

if size(step_penalties, 2) == 9
    if size(step_penalties,1) == length(modes)
        step_penalties = repmat(step_penalties, 1, size(score_matrix, 1));
    else
        step_penalties_unformatted = step_penalties;
        step_penalties = zeros(length(modes), 9*size(score_matrix,1));
        for ii = 1:length(modes):size(step_penalties_unformatted, 1)
            step_penalties(:,(9*(ii-1)+1):(9*ii)) = step_penalties_unformatted(ii:ii+length(modes)-1,:);
        end
    end
end

[alignment, bad_levels, cumulate_score, alignment_matrix] = ...
    alignmentC(score_matrix, step_penalties, modes, lookback, ...
    alignment_type, slip_locations, use_periodic_boundaries);

if alignment_type == 'S'

    alignment = alignment - sa_lookback + (1:length(alignment)) - 1;
    
    for cl = length(alignment):-1:1

        while alignment(cl) ~= alignment(alignment(cl))
            alignment(cl) = alignment(alignment(cl));
        end
    end
    
    largest_label = 0;
    for cl = 1:length(alignment)
        if alignment(cl) > largest_label + 1
            rest_of_alignment = alignment(cl:end);
            rest_of_alignment(rest_of_alignment == alignment(cl)) = largest_label + 1;
            alignment(cl:end) = rest_of_alignment;
        end
        
        largest_label = max(alignment(cl), largest_label);
        
    end
    
    if reorder
        alignment = renumberAlignment(alignment, step_penalties([1 2 4 6]));
    end
    
end


best_score = cumulate_score(end);
bad_levels = logical(bad_levels);



alignment(bad_levels) = nan;
                
                
end