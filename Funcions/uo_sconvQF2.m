function [f,g,h,xo] = uo_sconvQF2(n,seed)
    rng(seed);
    Q=randi(10,n);
    [V,D] = eig(triu(Q)+triu(Q)');
    Q=V*diag(diag(max(D,1)))*V';
    b=rand(n,1); xo = Q^-1*b;
    f = @(x) x'*Q*x/2-b'*x; g = @(x) Q*x-b; h = @(x) Q;
end