function s = smatrixACpar(x1, K1, x2, K2)
% calculates score matrix between vector data sets x1 and x2 with
% stiffnesses K1 and K2 in parallel

% pre-allocate the constant part
s = -0.5*size(x1, 1)*log(2*pi)*ones(size(x1, 2), size(x2, 2));
pifactor = -0.5*size(x1, 1)*log(2*pi);

n_1 = size(x1, 2);
n_2 = size(x2, 2);
n_elements = n_1 * n_2;

parfor cE = 1:n_elements
    c1 = mod(cE, n_1); 
    if c1 == 0
        c1 = n_1; 
    end
    c2 = ceil(cE / n_1);
    
    this_K1 = K1{c1};
    this_x1 = x1(:, c1);
    this_K1x1 = this_K1 * this_x1;
    this_x1K1x1 = this_x1' * this_K1x1;
    
    this_K2 = K2{c2};
    this_x2 = x2(:, c2);
    this_K2x2 = this_K2 * this_x2;
    this_x2K2x2 = this_x2' * this_K2x2;
    
    this_Kx = this_K1x1 + this_K2x2;
    this_K1K2 = this_K1 + this_K2;
    
    s(cE) = pifactor + 0.5*(logdet(this_K1) + logdet(this_K2) - logdet(this_K1K2) ...
        - this_x1K1x1 - this_x2K2x2 + this_Kx'/this_K1K2*this_Kx);
    
end

end
    
    