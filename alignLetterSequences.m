function [aligned, alignment, accuracy, misMatchInDelVec, confusionMat, accuracyDensity] = alignLetterSequences(refseq,calledseq,isglobal)

%a global alignment requires endpoints to match up. a non-global alignment
%can align the called sequence to a subsequence of the reference
if ~exist('isglobal', 'var')
    isglobal = true;
end

%make sure everything is the same case
refseq = upper(refseq);
calledseq = upper(calledseq);

%set penalties
pIndel = 1;
pMismatch = 1;
pMatch = -1;

%name lengths of sequences
L1 = length(refseq);
L2 = length(calledseq);

%the alignment matrix is 1 larger than the lengths of the sequences, to
%enforce boundary conditions.
H = zeros(L1+1,L2+1);

%if there is a global alignment, penalize skips at the beginning
if isglobal
    H(:,1) = pIndel*(0:L1)';
end

%always penalize skips anywhere in the called sequence
H(1,:) = pIndel*(0:L2)';

%the traceback matrix is the same size as the alignment matrix
T = zeros(L1+1,L2+1);

%the top row feeds back to the top left corner
T(1,:) = 1;

%if there's a global alignment, the right column feeds up to the top left
%corner. otherwise we terminate when we reach the left column.
if isglobal
    T(:,1) = 3;
else
    T(:,1) = -1;
end

%no matter what terminate at the top left corner
T(1,1) = -1;

for ii = 2:L1+1
    for jj = 2:L2+1
        if refseq(ii-1)==calledseq(jj-1)
            matchScore = pMatch;
        else
            matchScore = pMismatch;
        end
        options = [pIndel+H(ii,jj-1), matchScore+H(ii-1,jj-1), pIndel+H(ii-1,jj)];
        [H(ii,jj), T(ii,jj)] = min(options);
    end
end

%ii is already at the end of the alignment matrix. if we are doing a local
%alignment, we need to instead find the maximum value on the entire bottom
%row instead of just tracing back from the bottom right corner.
if ~isglobal
    [~,ii] = min(H(:,end));
end

%initialize the outputs
aligned = [];
alignment = [];

while true
    
    %terminate if we've hit a terminating boundary
    if T(ii,jj) == -1
        break
        
    %1, 2, and 3 correspond to stepping left, diagonal, and up respectively
    elseif T(ii,jj) == 1
        aligned(:,end+1) = [' '; ' '];
        aligned(2,end) = calledseq(jj-1);
        alignment(:,end+1) = [ii-1;jj-1];
        jj = jj-1;
    elseif T(ii,jj) == 2
        aligned(:,end+1) = [refseq(ii-1),calledseq(jj-1)];
        alignment(:,end+1) = [ii-1;jj-1];
        jj = jj-1;
        ii = ii-1;
    elseif T(ii,jj) == 3
        aligned(:,end+1) = [' '; ' '];
        alignment(:,end+1) = [ii-1;jj-1];
        aligned(1,end) = refseq(ii-1);
        ii = ii-1;
    end
    
end

%lose the 
clearvars H T
aligned = fliplr(aligned);
% aligned = [char(aligned); goodbad((aligned(1,:) == aligned(2,:))+1)];
aligned = char(aligned);
accuracy = sum(aligned(1,:) == aligned(2,:))/size(aligned,2);

% %% calculate mismatch and indel statistics
% inserts = aligned(1,:)~=aligned(2,:) & aligned(1,:)== ' ';
% insertfreq = sum(inserts)/size(aligned,2);
% deletes = aligned(1,:)~=aligned(2,:) & aligned(2,:)== ' ';
% deletefreq = sum(deletes)/size(aligned,2);
% mismatches = aligned(1,:)~=aligned(2,:) & aligned(1,:)~= ' ' & aligned(2,:)~= ' ';
% mismatchfreq = sum(mismatches)/size(aligned,2);
% misMatchInDelVec = [mismatchfreq,insertfreq,deletefreq];
% 
% aligned(end+1,:) = '+';
% aligned(end,inserts) = '.';
% aligned(end,deletes) = '-';
% aligned(end,mismatches) = 'x';
% 
% 
% alignedNum = zeros(2,length((aligned)));
% letts = ' ACGTX';
% for ii=1:length(letts)
%     alignedNum(1,(aligned(1,:)==letts(ii)))=ii;
%     alignedNum(2,(aligned(2,:)==letts(ii)))=ii;
% end
% confusionMat = zeros(length(letts));
% for ii=1:length(alignedNum)
%    confusionMat(alignedNum(1,ii), alignedNum(2,ii))= confusionMat(alignedNum(1,ii), alignedNum(2,ii))+1;
% end
% 
% confusionMat(end,:) = [];
% confusionMat(:,end) = [];
% 
% accuracyDensity = diff(alignedNum);
% accuracyDensity(accuracyDensity==0) =10; 
% accuracyDensity(accuracyDensity~=10)=0;
% accuracyDensity = accuracyDensity/10;

aligned_call = repmat(' ', 1, size(aligned, 2));
for cA = 1:size(aligned, 2);
    if aligned(1, cA) == aligned(2, cA)
        aligned_call(cA) = '+';
    elseif aligned(1, cA) == ' '
        aligned_call(cA) = '.';
    elseif aligned(2, cA) == ' '
        aligned_call(cA) = '-';
    elseif aligned(1, cA) ~= aligned(2, cA)
        aligned_call(cA) = 'x';
    end
end

aligned = [aligned ; aligned_call];

end





