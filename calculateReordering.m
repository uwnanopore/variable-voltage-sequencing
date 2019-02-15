function [level_placements] = calculateReordering(svm_scores)

% FUNCTIONALITY:
%   Solves for best possible rearrangement of levels given the svm-assessed
%   likelihoods of each transition being one of several possible step
%   types. Default implementation allows skips, backsteps, and steps.
%   Enforces that each level location must be visited exactly once. Returns
%   the placement 1:N_levels of the levels.
%
% DEFAULT BEHAVIOR:
%   The default state types will be 4: 
%   Step (1 S), 
%   Back (2 B), 
%   sKip|Step or sKip (3 K|SK),
%   sKip|Back (4 K|B)
%   The default allowed transitions force each level to be visited exactly
%   once, determining only the most like true ordering of levels
%   chronologically. Skips from skips are now allowed.
%
%   transition_matrix: (N_states) x (N_states) array. N_states may be
%                      different from N_stepsizes. Rows represent the
%                      "from" states, columns represent the "to" states.
%                      Entries should be 1s or 0s specifing if transition
%                      named by that entry is allowed
%
% INPUTS:
%   svm_scores: (N_levels - 1) x (N_stepsizes) array. Each row should
%               be the calculated SVM log transition probabilities for
%               that level transition being each of the different step
%               sizes. Default svm_scores column 1 = step, column 2 = 
%               skip, column 3 = backstep                   
%
% OUTPUTS: 
%   level_placements: (N_levels) x 1 array, with the best location
%                     placements for each of the levels. Re-indexed to
%                     begin always from 1 and count up to N_levels
%

% Set up default transition matrix
transition_matrix = [ ...
    1 0 1 0 ; ...
    0 0 0 1 ; ...
    1 1 1 0 ; ...
    1 0 1 0 ];
T = transition_matrix; clear transition_matrix

% Make sure svm_scores is oriented correctly
if size(svm_scores, 2) > size(svm_scores, 1)
    warning('SVM scores matrix seems to be oriented incorrectly, rotating.')
end
P = svm_scores; clear svm_scores

% Define some useful figures
n_states = size(T, 1);
n_stepsizes = size(P, 2);
n_transitions = size(P, 1);
n_levels = size(P, 1) + 1;

% Encode origin state information into T
T(T == 0) = nan;
T = T + ([0 : (n_states - 1)]');
T(isnan(T)) = n_states + 1;

% Initialize the Score Matrix (S) and Traceback Matrix (B)
S = zeros(n_transitions, n_states);
B = zeros(n_transitions, n_states);

% Fix first row of S and B
S(1, :) = [P(1, 1), P(1, 3), P(1, 2), -Inf];
B(1, :) = [NaN, NaN, NaN, NaN];

% Loop moves down through all transitions, following maximal paths
for cT = 2:n_transitions
    % define a scores vector, giving the score held at the start point of
    % each allowed pathway (entries 1:n_states) and disallowing the
    % disallowed pathways (entry n_states + 1)
    scores_ii = [S(cT - 1, :) -Inf];
    % fill a version of the transition matrix with the score for each
    % allowed transition, and -Inf for disallowed transitions
    T_ii = scores_ii(T);    
    % find value and row of maximal score in each column of T_ii
    [max_val, max_loc] = max(T_ii);
    % define base values to add into new row of S
    base_val = [P(cT, 1), P(cT, 3), P(cT, 2), P(cT, 2)];
    % fill new row in S
    S(cT, :) = base_val + max_val;
    % fill in traceback pointing back to where each new entry came from
    B(cT, :) = max_loc;
end

% Initialize storage for best path
path = zeros(n_transitions, 1);

% Finding max of bottom row of S tells us where the best path ends and
% where to start our traceback
[~, best_loc] = max(S(end, :));

% Store this end point in our path, and generate a pointer back to where
% the path came from
path(end) = best_loc;
pointer = B(end, best_loc);

% Follow back through traceback matrix B to reconstruct best path
for cT = (n_transitions - 1):-1:1
    path(cT) = pointer;
    pointer = B(cT, pointer);
end

% Convert path entries into transition types
transition_translator = [1, -1, 2, 2];
steps = transition_translator(path);

% Need to solve for location of each level. "steps" tells us the difference
% between each pair of adjacent levels (level5 - level4 = steps(4)). Be
% default, set level1 = 1. 

% Build matrix to encode transition differences
X = sparse(diag(-1 .* ones(1, n_transitions)));
X = sparse([[zeros(1, n_transitions) ; X] zeros(n_levels, 1)]);
X = sparse(X + diag(ones(1, n_levels)));

% Determine location of each level by inverting X
level_placements = X \ ([1 steps]');

% Re-index level_placements to start at 1
level_placements = level_placements + (1 - min(level_placements));

end