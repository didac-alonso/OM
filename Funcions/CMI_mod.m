%{
H: H(x)
B: Hessiana modificada
%}
function [B,tau, mod] = CMI_mod(H)
    tau = 0;    
    laUB = norm(H, 'fro');
    k = 0;
    n = size(H,1);
    [R, err] = chol(H);
    mod = 0;
    l = 0;
    while(err > 0)
        l = l+1;
        tau = (1.01-1/2^l)*laUB;
        [R, err] = chol(H+tau*eye(n));
        mod = 1;
    end
    B = R'*R;
end