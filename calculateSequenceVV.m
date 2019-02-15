function [sequence, total_score, kmer_calls, kmer_confidence, kmer_half_steps, level_was_good, source_levels] = calculateSequenceVV(x, K, p, p_bad, varargin)
%
% outputs:
%   sequence
%       the maximum likelihood sequence
%   total_score
%       the total accumulated alignment score of the sequencing run
%   kmer_calls
%       called map state matched by each level
%   kmer_confidence
%       confidence in choice of each map kmer we match to
%   kmer_half_steps
%       the called step size of each level
%   level_was_good
%       boolean whether or not the level was used (good level) or discarded
%       (bad level)
%   source_levels
%       level most responsible for an individul base call
%
% inputs
%   x
%       the d x n array of feature vectors
%   K
%       the 1 x n cell array of the d x d stiffness matrices on each
%       feature vector
%   p
%       an n x (b + 1) x 3 array whose 3rd dim elements are [p(step),
%       p(skip), p(skip extend)]. for example, p(12, 3, 1) is the step
%       probability formt he 12-3 = 9th level to the 12th level. the 2nd
%       dimension size is set by the maximum number of sequential bad
%       levels allowed to be found by the sequencer; it is 1 less than the
%       dimension (equal to b in the above expression)
%   p_bad
%       an n x 1 array of bad level probabilities for each level
%
% varargin
%   backsteps
%       an n x 1 logical array determining whether each level was observed
%       to have backstepped
%   map
%       the kmer map to be used in decoding
%   transitioninfo
%       the structure of transition information (which map states can
%       transition to which other map states, and with which associated
%       penalties)
%   pindback
%       probabiliity that a level is an ATP-independent step given that we
%       observed it to backstep
%   pback
%       the probability of a backstep
%   scorecutoff
%       how much worse a match score can be than the best scoring match in
%       a row for it to still be considered. this is a negative number.
%       making its magnitude larger results in a more comprehensive search,
%       but takes longer

% set default parameters
p_ind_given_backstep = 0.975;
p_backstep = 0.025;
score_cutoff = -10;
verbose = false;

% handle varargin
for ca = 1:length(varargin)
    if ischar(varargin{ca})
        switch upper(varargin{ca})
            case 'MAP'
                map = varargin{ca+1};
            case 'TRANSITIONINFO'
                transition_info = varargin{ca+1};
            case 'BACKSTEPS'
                did_backstep = varargin{ca+1};
            case 'PINDBACK'
                p_ind_given_backstep = varargin{ca+1};
            case 'PBACK'
                p_backstep = varargin{ca+1};
            case 'SCORECUTOFF'
                score_cutoff = varargin{ca+1};
            case 'VERBOSE'
                verbose = true;
            case 'QUIET'
                verbose = false;
        end
    end
end

% load critical structures if not already specified
if ~exist('transition_info', 'var')
    transition_info = load('transition_info_hel308_6mer.mat');
    transition_info = transition_info.transition_info;
end

if ~exist('map', 'var')
    map = load('pore_model_6mer_variable_voltage.mat');
    map = map.model;
end

% find the backstep priors
% if backstep variables not passed in...
if ~exist('did_backstep', 'var')
    % record that there are no known backsteps and that we have
    % uninformative priors on this variable
    did_backstep = zeros(size(x, 2), 1);
    prior_given_backstep = 0.5 .* ones(1, size(map.mean, 2));
    prior_given_nobackstep = prior_given_backstep;
    
else
    % the probability that a level is an ATP-independent step given that we
    % did not see it backstep:
    p_ind_given_nobackstep = 0.5 - (p_ind_given_backstep * p_backstep);
    
    % an 8192 x 1 array of prior probabilities for each kmer given that a
    % level backstepped
    prior_given_backstep = repmat([p_ind_given_backstep ; 1 - p_ind_given_backstep]', 1, size(map.mean, 2) / 2);
    
    % the same thing but given that the level did not backstep
    prior_given_nobackstep = repmat([p_ind_given_nobackstep ; 1 - p_ind_given_nobackstep]', 1, size(map.mean, 2) / 2);
    
end

% name some constants
% the max number of steps we can jump before completely leaving the kmer
n = size(transition_info.p_combos, 2);

% the largest number of consecutive bad levels
max_sequential_bad = size(p, 2) - 1;

% calculate the step type transition values for sequencing
% divergent workflow depending on nature of the input p
% v1: p was n_levels x (b + 1) x 3, with the 3rd dimension representing
% step / skip / skip extend probabilities
% v2: p was n_levels x (b + 1) x 12, as generated by smartStepCounts.m,
% where the 3rd dimension now represents step1 / step2 / ... / step12
if size(p, 3) == 3
    % v1...
    % calculate the combined extension penalties for 0 through n - 1 skip
    % extensions
    p_relative_skip = (p(:,:,3).^reshape(0:n-2, 1, 1, n-1))./(sum(p(:,:,3).^reshape(0:n-2, 1, 1, n-1),3));
    
    % build the pvals matrix out of the known penalties and the known
    % structure of the transition matrix. this is a list of all unique
    % possible transition penalties, computed for each (row, col) of p.
    p_vals = zeros(size(x,2),max_sequential_bad+1,size(p_relative_skip,3)+1);
    p_vals(:,:,1) = p(:,:,1);
    p_vals(:,:,2:size(p_relative_skip,3)+1) = p(:,:,2).*p_relative_skip;
    
elseif size(p, 3) == 12
    % v2...
    % p translates directly to the p_vals we want
    p_vals = p;
    % check normalization
    for cB = 1:size(p_vals, 2)
        total_probability = sum(p_vals(:, cB, :), 3);
        total_probability = repmat(total_probability, 1, 1, size(p_vals, 3));
        p_vals(:, cB, :) = p_vals(:, cB, :) ./ total_probability;
    end
else
    error('Dimension of P is incorrect.\n');
end

% change p_vals to implement the p_combos from transition_info
p_vals = log(permute(  sum(p_vals.*permute(transition_info.p_combos,[3,4,2,1]), 3)  ,[1,2,4,3]));

%find the penalty to be applied for each bad level.
p_good = log(1 - p_bad);
p_bad = log(p_bad);

%calculate the score matrix, including priors based on backstep info
alignment_matrix = smatrixACpar(x,K,map.mean,map.stiffness) + prior_given_backstep.*(did_backstep) + prior_given_nobackstep.*(~did_backstep);

%prune the score matrix and keep track of which matches we want to consider
candidate_matches = alignment_matrix > max(alignment_matrix,[],2) + score_cutoff;
num_candidate_matches = sum(candidate_matches,2);
alignment_matrix(~candidate_matches) = -inf;

%initialize the traceback matrices- one for which row we came from, one for
%which column.
traceback_matrix_rows = zeros(size(alignment_matrix));
traceback_matrix_cols = zeros(size(alignment_matrix));

if verbose
    disp('performing traceback')
end

%loop over measured levels
for cl = 2:size(alignment_matrix,1)
    
    %track which level we are on
    if mod(cl, 10) == 0 && verbose
        disp(cl)
    end
    
    %get the transition penalties for this level
    %this is a little tricky- it turns out that there are only 48 unique
    %numbers in the 8192x8192 transition matrix, and these are different
    %sums of powers of the step penalties.
    
    %loop over each row we are considering jumping from
    this_row_best_scores = -inf(1,num_candidate_matches(cl));
    this_row_best_row_from = zeros(1,num_candidate_matches(cl));
    this_row_best_col_from = zeros(1,num_candidate_matches(cl));
    
    for cfrom = 0:min(max_sequential_bad, cl-2)
        %calculate the transition matrix for each "from" row, including the
        %bad level penalties accrued by skipping the rows between this
        %"from" row and the present "to" row.
        transition_submatrix = reshape(p_vals(cl,cfrom+1,transition_info.perms_matrix(candidate_matches(cl-cfrom-1,:), candidate_matches(cl,:))) + sum(p_bad((cl - cfrom):cl-1)), num_candidate_matches(cl-cfrom-1), num_candidate_matches(cl)) + p_good(cl - cfrom - 1);
        
        %find the best "from" option for each "to"
        [from_score, where_from] = max(transition_submatrix + alignment_matrix(cl-cfrom-1,candidate_matches(cl-cfrom-1,:))', [], 1);
        
        %find the "to" elements that are better than the best one so far
        better_scores = from_score > this_row_best_scores;
        from_indices = find(candidate_matches(cl-cfrom-1,:));
        
        %update the best-one trackers with this information
        this_row_best_col_from(better_scores) = from_indices(where_from(better_scores));
        this_row_best_row_from(better_scores) = cl - cfrom - 1;
        this_row_best_scores(better_scores) = from_score(better_scores);
    end

    %update the alignment and traceback matrices with the best options
    alignment_matrix(cl,candidate_matches(cl,:)) = this_row_best_scores + alignment_matrix(cl,candidate_matches(cl,:));
    traceback_matrix_rows(cl,candidate_matches(cl,:)) = this_row_best_row_from;
    traceback_matrix_cols(cl,candidate_matches(cl,:)) = this_row_best_col_from;
    
end

%find the starting point for the traceback
[total_score,kmer] = max(alignment_matrix(end,:));
alignment_matrix = alignment_matrix - jlog(alignment_matrix);

%make sure we start at the last level
cl = size(alignment_matrix,1);
kmer_level_source = cl;

%initialize the step size array, and 
num_half_steps = nan(1, size(alignment_matrix,1));

%the first step is always a single step by convention
num_half_steps(1) = 1;

%traceback until you get to the first level
local_score = [];
level_was_good = false(1,size(x,2));

%construct the sequence by walking backwards
while cl > 1
    
    %remember the last kmer we were on
    last_kmer = kmer(1);
    
    level_was_good(cl) = true;
    %the most likely step size to get from the last kmer to this kmer
    num_half_steps(cl) = transition_info.step_size_matrix(traceback_matrix_cols(cl,last_kmer), last_kmer);
    
    %append the previous kmer with appropriate padding in 3' to 5' order
    kmer = [traceback_matrix_cols(cl,last_kmer), nan(1, num_half_steps(cl)-1), kmer];
    kmer_level_source = [cl, nan(1, num_half_steps(cl)-1), kmer_level_source];

    %track the latest score
    local_score = [alignment_matrix(cl,last_kmer), local_score];

    %step back a level
    cl = traceback_matrix_rows(cl,last_kmer);
    
end

level_was_good(1) = true;
local_score = [local_score traceback_matrix_cols(cl,kmer(1))];
kmer_calls = zeros(1,size(x,2));
kmer_calls(level_was_good) = kmer(~isnan(kmer));
kmer_confidence = p_bad;
kmer_confidence(level_was_good) = local_score;
kmer_half_steps = num_half_steps;

%because hel308 reads 3' to 5'
kmer = fliplr(kmer);
kmer_level_source = fliplr(kmer_level_source);

%find the first and last good levels and crop the sequence
first_non_nan = find(~isnan(kmer), 1, 'first');
last_non_nan = find(~isnan(kmer), 1, 'last');
kmer = kmer(first_non_nan:last_non_nan);
kmer_level_source = kmer_level_source(first_non_nan:last_non_nan);

%jumps of 3+ nucleotides could be either jumps of 7 or 8 bases, and these
%can only be distinguished by the necessity of preserving the fact that we
%are going "even-odd-even-odd-even-odd-etc" in the kmer sequence. Using
%this constraint we find all places with 7 NaNs and add a NaN if we need
%one to preserve parity.

%find the places where the long gaps begin
long_nan_starts = strfind(isnan(kmer), ones(1,n-1));

%loop over these locations
for ii = 1:length(long_nan_starts)
    
    %the current one
    lns = long_nan_starts(ii);
    
    %calculate the parity between the last good level and the next good
    %level to determine whether to insert a nan
    parity = mod(kmer(lns-1) + kmer(lns + n-1), 2);
    
    %add the nan if necessary
    kmer = [kmer(1:lns-1) nan(1,parity) kmer(lns:end)];
    kmer_level_source = [kmer_level_source(1:lns-1) ...
                    repmat(kmer_level_source(lns),1,parity) kmer_level_source(lns:end)];
                
    %change the indexing of the long nan starts if necessary
    long_nan_starts = long_nan_starts + parity;
end

%find the first odd (pre) step
first_odd_one = 2-mod(find(mod(kmer, 2)==1, 1, 'first'), 2);

%go two by two, making it so we only have even kmers
for ii = first_odd_one:2:length(kmer)
    if ~isnan(kmer(ii))
        kmer(ii:min(ii+1, length(kmer))) = kmer(ii)+1;
    elseif ii < length(kmer)
        kmer(ii:min(ii+1, length(kmer))) = kmer(ii+1);
    end 
end

%count each determined kmer once
kmer = kmer(1:2:end);
kmer_level_source = reshape([kmer_level_source(1:end) nan(1, mod(length(kmer_level_source),2))], 2, []);

% call bases and track level -> base mapping
sequence_assembly = repmat(' ', numel(kmer), numel(kmer)*6);
sequence_assembly(1,1:6) = char(map.name(kmer(1),:));
write_position = [1,1];
for ii = 2:length(kmer)
        write_position = write_position + [1,1];
    if ~isnan(kmer(ii))
        sequence_assembly(write_position(1), write_position(2):write_position(2)+5) = char(map.name(kmer(ii),:));
    end
end

sequence = char(max(sequence_assembly,[],1));
sequence(sequence == ' ') = [];
sequence_assembly = sequence_assembly(:,1:length(sequence));
contributing_levels = sequence_assembly == sequence;
source_levels = cell(1,length(sequence));

for ii = 1:length(sequence)
    ix = find(contributing_levels(:,ii));
    for tt = 1:length(ix)
        jj = ix(tt);
        source_levels{ii} = reshape(unique([source_levels{ii} kmer_level_source(:,jj)']), 1, []);
    end
    source_levels{ii}(isnan(source_levels{ii})) = [];
end

end