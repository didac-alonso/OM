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
function [xk, it] = GM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax)
    it = 0;
    pp = p;
    xk = [x];
    dfk = df(x);
    while norm(dfk) > tol & it < itmax
        d = - df(x);
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW);
        % Si no troba alpha, cal fer més exhaustiva la cerca
        if iWout < 2
            p = 1.1*p;
        else
            x = x + a*d;
            xk = [xk x];
            it = it + 1;
            p = pp;
            dfk = df(x);
        end
    end
end