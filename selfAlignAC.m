function [alignment] = selfAlignAC(X, K, T, PSA, lookback)

% FUNCTIONALITY:
%
% INPUTS:
%   X: [N_features x N_levels] array of level feature values
%   K: [1 x N_levels] cell array of [N_features x N_features] arrays of level
%      feature stiffnesses
%   T: Either a [1 x 6] array or [N_levels x 6] array of transition
%      probabilies. If dim(1) = 1, then these are assumed to apply to all
%      transitions equally. If dim(1) = N_levels, then these are
%      interpreted as transition penalties for each transition. In this
%      case, the first row corresponds to the "null" transition into the
%      1st level. Column meanings are                   
%      [step, skip, skip+, back, back+, hold]
%   PSA: Either a [1 x 1] array or an [N_levels x 1] array of log "create
%        new state" penalties. If only 1, then applies to all states,
%        otherwise if N_levels, then applies to states individually
%   lookback: Number of levels to look back for self alignment. Default is
%             10 levels.
%
% OUTPUTS:
%   alignment: [1 x N_levels] array detailing the input levels -> 
%              filtered levels alignment
%

%% Assign defaults, check dimensions
% extract important dimensions
N_levels = size(X, 2);

% set default lookback if not specified
if ~exist('lookback', 'var')
    lookback = 10;
end

% repmat PSA if singlular value
if size(PSA, 1) == 1
    PSA = repmat(PSA, N_levels, 1);
end

% repmat T if singular row
if size(T, 1) == 1
    T = repmat(T, N_levels, 1);
end

%% calculate score matrix
S = smatrixACBF(X, K, X, K, 'selfalignment', lookback);

%% pre-calculate all transitions (back2+, skip2+ up to lookback)
% calculate multi-step/back penalties
Tback = repmat(T(:, 4), 1, lookback) .* (repmat(T(:, 5), 1, lookback) .^ repmat([0:1:(lookback - 1)], N_levels, 1));
Tskip = repmat(T(:, 2), 1, lookback) .* (repmat(T(:, 3), 1, lookback) .^ repmat([0:1:(lookback - 1)], N_levels, 1));
% rebuild a T matrix [backLB, ..., back1, hold, step, skip1, ..., skipLB]
T = [fliplr(Tback) T(:, 6) T(:, 1) Tskip];
% tack on a zeros column to the end of T for disallowed transitions 
% index = (:, 2*lookback + 3)
T = [T repmat(0, size(T, 1), 1)];

%% figure out indexing for step types allowed into each column (non-new-state)
allowed_steps = cell(1, lookback);
for cS = 1:size(allowed_steps, 2)
    allowed_steps{cS} = (cS - lookback):1:(cS);
    allowed_steps{cS} = fliplr(allowed_steps{cS});
    allowed_steps{cS} = allowed_steps{cS} + lookback + 1;
end

% additional figuring for entering from the sKip* and Back* states
for cS = 1:size(allowed_steps, 2)
    % filling in entries for coming from sKip* state
    if allowed_steps{cS}(lookback + 1) < (lookback + 1)
        allowed_steps{cS} = [allowed_steps{cS} (allowed_steps{cS}(lookback + 1) - 1)];
    elseif allowed_steps{cS}(lookback + 1) == (lookback + 1)
        allowed_steps{cS} = [allowed_steps{cS} (lookback + 1)];
    end
    
    % filling in entries for coming from Back* state
    if allowed_steps{cS}(lookback + 2) < (lookback  - 1)
        allowed_steps{cS} = [allowed_steps{cS} (allowed_steps{cS}(lookback + 2) + 2)];
    elseif allowed_steps{cS}(lookback + 2) == (lookback - 1)
        allowed_steps{cS} = [allowed_steps{cS} (lookback + 2)];
    elseif allowed_steps{cS}(lookback + 2) == (lookback + 1)
        allowed_steps{cS} = [allowed_steps{cS} (lookback + 1)];
    end
end

%% figure out indexing for step types allowed into new-state columns
into_step_allowed = [repmat(lookback + 2, 1, lookback + 2), (2 * lookback) + 3];
into_skip_allowed = [repmat(lookback + 3, 1, lookback + 3)];
into_back_allowed = [repmat((2 * lookback) + 3, 1, lookback + 1), lookback, (2 * lookback) + 3];

%% fill in alignment matrix
% initialize alignment matrix
M = nan(N_levels, lookback + 3);
% initialize traceback matrix
B = nan(N_levels, lookback + 3);

% loop through alignment matrix, filling in row by row
for cR = 1:size(M, 1)
    for cC = 1:size(M, 2)
        % if first row, we are not "coming from" anywhere, and we can only
        % "be in" the self state (entry 1, lookback + 1)
        if cR == 1
            if cC ~= (lookback + 1)
                M(cR, cC) = -inf;
            elseif cC == (lookback + 1)
                M(cR, cC) = ...
                    S(cR, cC) ...
                    + PSA(cR);
            end
            continue
        end
        
        % if we are in the right-most 3 columns (new-level columns Step*,
        % sKip*, and Back*) we must account for the self alignment penalty
        if cC == (lookback + 1)
            % column (lookback + 1) is "new via step"
            % where came from
            [origin_score, origin] = max(M(cR - 1, :) + log(T(cR, into_step_allowed)));
            B(cR, cC) = origin;
            % new matrix element value
            % 1: best entry + path from previous row
            % 2: normalize for weight of ways in
            % 3: self alignment (new state) penalty
            % 4: level alignment score
            M(cR, cC) = ...
                origin_score ...
                + PSA(cR) ...
                + S(cR, cC);
            % - log(sum((M(cR - 1, :) > -inf) .* T(cR, into_step_allowed))) ...
            continue
        elseif cC == (lookback + 2)
            % column (lookback + 2) is "new via skip"
            % where came from
            [origin_score, origin] = max(M(cR - 1, :) + log(T(cR, into_skip_allowed)));
            B(cR, cC) = origin;
            % new matrix element value
            % 1: best entry + path from previous row
            % 2: normalize for weight of ways in
            % 3: self alignment (new state) penalty
            % 4: level alignment score
            M(cR, cC) = ...
                origin_score ...
                + PSA(cR) ...
                + S(cR, cC - 1);
            % - log(sum((M(cR - 1, :) > -inf) .* T(cR, into_skip_allowed))) ...
            continue
        elseif cC == (lookback + 3)
            % column (lookback + 3) is "new via back"
            % where came from
            [origin_score, origin] = max(M(cR - 1, :) + log(T(cR, into_back_allowed)));
            B(cR, cC) = origin;
            % new matrix element value
            % 1: best entry + path from previous row
            % 2: normalize for weight of ways in
            % 3: self alignment (new state) penalty
            % 4: level alignment score
            M(cR, cC) = ...
                origin_score ...
                + PSA(cR) ...
                + S(cR, cC - 2);
            % - log(sum((M(cR - 1, :) > -inf) .* T(cR, into_back_allowed))) ...
            continue
        end
        
        % otherwise, we do standard alignment stuff to fill where came from
        % where came from
        [origin_score, origin] = max(M(cR - 1, :) + log(T(cR, allowed_steps{cC})));
        B(cR, cC) = origin;
        % new matrix element value
        % 1: best entry + path from previous row
        % 2: normalize for weight of ways in
        % 3: level alignment score
        M(cR, cC) = ...
            origin_score ...
            + S(cR, cC);
        % - log(sum((M(cR - 1, :) > -inf) .* T(cR, allowed_steps{cC}))) ...
    end
end

%% perform traceback to construct alignment
traceback = nan(1, N_levels);

[~, traceback(N_levels)] = max(M(N_levels, :));
for cL = (N_levels - 1):-1:1
    traceback(cL) = B(cL + 1, traceback(cL + 1));
end

%% convert traceback within self alignment into an f_to_tf alignment
alignment = 1:N_levels;
traceback = traceback - lookback - 1;
markbacks = traceback == 2;
markbacks = cumsum(markbacks);
alignment = alignment + markbacks;
traceback(traceback == 2) = -2;
alignment = alignment + traceback;

%% make sure there are no unused indices in the alignment
for cL = max(alignment):-1:1
    if ~any(alignment == cL)
        alignment(alignment > cL) = alignment(alignment > cL) - 1;
    end
end

end