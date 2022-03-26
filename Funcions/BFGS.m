function [xk, dk, ak, Hk, tauk, iWk, it] = BFGS(x, f, df, amin, amax, rho, c1, c2, iW, tol, itmax, Q)
    I = eye(size(x,1));
    H = I;
    Hk = [H];
    it = 1;
    xk = [x];
    dk = [];
    ak = [];
    iWk = [];
    dfx = df(x);
    while norm(dfx) > tol & it <= itmax
        d = -H*dfx;
        dk = [dk d];
        [a, iWout] = BLS(x, f, df, d, amin, amax, rho, c1, c2, iW, Q); 
        ak = [ak a];
        iWk = [iWk iWout];
        x = x+a*d;
        dfx1 = dfx;
        dfx = df(x);
        s = x-xk(:,end); y = dfx-dfx1; 
        p = 1/(y'*s);
        xk = [xk x];
        H = (I-p*s*y')*H*(I-p*y*s')+p*(s*s');
        Hk = cat(3, Hk, H); % concatenacio 3D
        it = it + 1;
    end
end