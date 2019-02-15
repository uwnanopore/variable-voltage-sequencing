function ld = logdet(A)
%computes the log determinant of a symmetric, positive definite matrix A

if ~issymmetric(A)
    error('A is not symmetric');
end

[R,p] = chol(A);

if p
    error('A is not positive definite');
end

ld = 2*sum(log(diag(R)));

end