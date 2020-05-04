function [ b, yy ] = inversefit( x,y, xx,n )
N = numel(x);
NN = numel(xx);


xinv = 1./x.^n;
xxinv = 1./xx.^n;
X = [ones(N,1), xinv'];
XX = [ones(NN,1), xxinv'];
y = y';
b = (X'*X)^-1*X'*y;
yy = XX*b;

end

