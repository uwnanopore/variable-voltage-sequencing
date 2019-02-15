function antisensesequence = antisense(sequence)

sequence = fliplr(sequence);

as = find(sequence == 'A');
cs = find(sequence == 'C');
gs = find(sequence == 'G');
ts = find(sequence == 'T');

antisensesequence(as) = 'T';
antisensesequence(cs) = 'G';
antisensesequence(gs) = 'C';
antisensesequence(ts) = 'A';

end