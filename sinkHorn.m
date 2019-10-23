function[d] = sinkHorn(X,Y,lambda,r)
%[d] = sinkHorn(M,lambda,r,C)
%

if nargin == 0. % Example of how to use it
    nd = 100;
    ep =0.5;
    Xo = randn(nd,2);  
    X = randn(nd,2);
    Y = 3*Xo./sqrt(sum(Xo.^2,2)+ep);

    t = linspace(2.0,4.0,201);
    figure(1)
    for i=1:length(t)
        Z = t(i)*X./sqrt(sum(X.^2,2)+ep); % + (t(i)-3);
        plot(Y(:,1),Y(:,2),'.r','markersize',30);
        hold on
        plot(Z(:,1),Z(:,2),'.b','markersize',30);
        axis([-5,5,-5,5])
        hold off
        drawnow
        domt(i)  = sinkHorn(Y,Z,100,ones(nd,1));
    end
    figure(2)
    plot(t,domt,'r','linewidth',3);
    [~,ii] = min(domt);
    hold on;
    plot(t(ii),domt(ii),'xb','markersize',20)
    title(['min = ',num2str(t(ii))]) 
    hold off
    return
end

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
