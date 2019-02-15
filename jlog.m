function D = jlog(d)
%computes the jacobian logarithm of a series {d_1, d_2, ..., d_n}
%D = log( sum_i ( exp (d_i) ) )
%using an approximation to avoid overflow or roundoff error

D = d(:, 1);
for ii = 2:size(d,2)
    D = max(D, d(:,ii)) + correctionFunction(D, d(:, ii));
end


end

function cf = correctionFunction(x,y)

fixes = (x == -inf & y == -inf);

cf = log(1 + exp(-abs(x - y)));
cf(fixes) = 0;


end