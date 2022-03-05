%{
Busca zeros a la derivada d'una funció.

SOL VÀLIDA PER 1 VARIABLE

Input:
    - xk: Punt inicial
    - f, df: funció i derivada
    - amin, amax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer:
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme

Output:
    - x: Punt que fa 0 df.
    - it: nombre d'iteracions usades
%}
function [xk, dk, ak, iWk, it] = GM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax)
    it = 1;
    pp = p;
    xk = [x];
    dk = [];
    dfxk = df(x);
    ak = [];
    iWk = [];
    while norm(dfxk) > tol & it <= itmax
        d = - dfxk;
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW);
        iWk = [iWk iWout];
        ak = [ak a];
        % Si no troba alpha, cal fer més exhaustiva la cerca
        if iWout < 2
            p = 1.1*p;
        else
            x = x + a*d;
            xk = [xk x];
            it = it + 1;
            p = pp;
            dfxk = df(x);
            dk = [dk d];
        end
    end
end