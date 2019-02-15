function [x, k, e] = levelPredPipelineDC(sequence, map)

bases = 'ACGT';

baseMult = length(bases) .^ (0:(size(map.name, 2) - 1));
numseq = nan(1, length(sequence));
for cB = 1:length(bases)
    numseq(sequence == bases(cB)) = cB - 1;
end

hexIDs = conv(numseq, baseMult, 'valid') + 1;
hexIDs = reshape([hexIDs * 2 - 1; hexIDs * 2], 1, []);

x = nan(size(map.mean, 1), length(hexIDs));
e = nan(size(map.mean, 1), length(hexIDs));
k = cell(1, length(hexIDs));

for cQ = 1:length(hexIDs)
    if ~isnan(hexIDs(cQ))
        x(:, cQ) = map.mean(:, hexIDs(cQ));
        e(:, cQ) = map.error(:, hexIDs(cQ));
        k{cQ} = map.stiffness{hexIDs(cQ)};
    end
end

end