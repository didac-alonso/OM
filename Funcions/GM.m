%{
Busca zeros a la derivada d'una funció.

Input:
    - xk: Punt inicial
    - f, df: funció i derivada
    - amin, amax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer, si iW= 0, usem ELS, iW = 1 WC, iW =
    2 SWC
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme

Output:
    - x: Punt que fa 0 df.
    - it: nombre d'iteracions usades
%}
function [xk, dk, ak, iWk, it] = GM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax, Q)
    it = 1;
    xk = [x];
    dfx = df(x);
    dk = [-dfx];
    d = -dfx;
    ak = [];
    iWk = [];
    while norm(dfx) > tol & it <= itmax
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW, Q);
        iWk = [iWk iWout];
        ak = [ak a];

        x = x + a*d;
        xk = [xk x];
        it = it + 1;
        dfx = df(x);
        d = - dfx;
        dk = [dk d];
        %end
    end
end