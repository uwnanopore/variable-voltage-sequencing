function [x_ds,K_ds] = downsampleFeatures(x,K,n)
%downsample an array of features (x) and a cell array of stiffness matrices
%(K) to only n features, weighted by errors



startix = ceil(linspace(1, size(x,1), n+1));
endix = startix(2:end);
endix(1:end-1) = endix(1:end-1) - 1;
startix = startix(1:end-1);

x_ds = zeros(n, size(x,2));
K_ds = cell(1, size(x,2));

for cl = 1:size(x,2)
    
    K_ds{cl} = zeros(n);
    
    for ii = 1:length(startix)
        ix1 = startix(ii):endix(ii);
        for jj = 1:length(startix)
            ix2 = startix(jj):endix(jj);
            K_ds{cl}(ii,jj) = sum(sum(K{cl}(ix1,ix2)));
        end
        
        x_ds(ii,cl) = mean(x(ix1,cl));
    end
    
end

end