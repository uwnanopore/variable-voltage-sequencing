function [mu, K_mu] = multivarSampleStats(x, K_x, nosparse)
%x is a DxN set of N measurements of a D dimensional random vector
%K_x is a cell array of the inverse covariances of each measurement
%mu is the weighted sample mean of these measurements
%K_mu is the weighted sample inverse covariance of these measurements

if ~exist('nosparse', 'var')
    nosparse = false;
end

if issparse(K_x{1})
    for cK = 1:length(K_x)
        K_x{cK} = full(K_x{cK});
    end
end

K_x = cell2mat(reshape(K_x,1,1,[]));
K_mu = sum(K_x,3)/size(x,2);

try
mu = zeros(size(x,1),1);
count = zeros(size(x,1),1);
for ii = 1:size(x,2)
    goods = ~isnan(x(:,ii));
    n_good = sum(goods);
    mu(goods) = mu(goods) + K_x(goods,goods,ii)*x(goods,ii)*n_good/size(x,1);
    count(goods) = count(goods) + n_good/size(x,1);
end
mu = K_mu\(mu./count);
catch
    disp('ok')
end

if ~nosparse
    K_mu = sparse(K_mu);
elseif nosparse
    K_mu = full(K_mu);
end

end