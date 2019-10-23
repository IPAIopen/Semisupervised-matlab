function[d] = sinkHorn(X,Y,lambda,r)
%[d] = sinkHorn(M,lambda,r,C)
%


M = distMat(X,Y);
C = ones(size(M));
K = exp(-lambda*M);
U = ones(length(r), length(r) )/length(r);
K  = diag(1./r)*K;
for i=1:20
    U = 1./(K*(C./(K'*U))); 
end 
V = C./(K'*U);
d = sum(sum(U .* ((K .* M )*V)));

end

function A = distMat(X,Y)
    A = ((sum(X.^2,2) + sum(Y.^2,2)') - 2*X*Y');
end
