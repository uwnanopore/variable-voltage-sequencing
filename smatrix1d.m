function S = smatrix1d(x1,K1,x2,K2)

v1 = 1./cell2mat(K1);
v2 = 1./cell2mat(K2);
V = v1' + v2;
X = (x1'-x2).^2;

S = -0.5*(log(2*pi*V) + X./V);


end