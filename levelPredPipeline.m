function [x, k, total_weight] = levelPredPipeline(sequence, map, gobig)

if ~exist('gobig', 'var')
    gobig = 0;
end

bases = 'ACGT';

baseMult = length(bases) .^ (0:(size(map.name, 2) - 1));
numseq = nan(1, length(sequence));
for cB = 1:length(bases)
    numseq(sequence == bases(cB)) = cB - 1;
end

hexIDs = conv(numseq, baseMult, 'valid') + 1;
hexIDs = reshape([hexIDs * 2 - 1; hexIDs * 2], 1, []);

if gobig == 0
    x = nan(size(map.mean_3, 1), length(hexIDs));
    k = cell(1, length(hexIDs));
    total_weight = nan(1, length(hexIDs));
elseif gobig == 1
    x = nan(size(map.mean_101, 1), length(hexIDs));
    k = cell(1, length(hexIDs));
    total_weight = nan(1, length(hexIDs));
end

for cQ = 1:length(hexIDs)
    if ~isnan(hexIDs(cQ))
        if gobig == 0
            x(:, cQ) = map.mean_3(:, hexIDs(cQ));
            k{cQ} = map.stiffness_3{hexIDs(cQ)};
            if isfield(map, 'total_weight')
                total_weight(cQ) = map.total_weight(hexIDs(cQ));
            end
        elseif gobig == 1
            x(:, cQ) = map.mean_101(:, hexIDs(cQ));
            k{cQ} = map.stiffness_101{hexIDs(cQ)};
            if isfield(map, 'total_weight')
                total_weight(cQ) = map.total_weight(hexIDs(cQ));
            end
        end
    end
end

end