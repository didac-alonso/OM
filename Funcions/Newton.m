%{
Busca zeros a la derivada d'una funció.

SOL VÀLIDA PER 1 VARIABLE

Input:
    - xk: Punt inicial
    - f, df, d2f: funció, derivada i segona derivada
    - amin, amax: rang per les alphas, s'usa en el BLS
    - c1, c2: constants per les Wolfe Conditions
    - iW: Condicions Wolfe a satisfer:
    - tol: precisió buscada
    - itmax: màxim d'iteracions que pot dur a terme l'algoritme

Output:
    - xk: Punt que fa 0 df.
    - it: nombre d'iteracions usades
%}
function [x,it]=Newton(xk, f, df, d2f, amin, amax, p, c1, c2, iW, tol,itmax)
    it=1;
    x = [xk];
    dfk = df(xk);
    while norm(dfk)>tol & it<=itmax
        d = -d2f(xk)\dfk; % Equivalent a d2f(xk)^-1*dfk, però és més eficient així
        a = 1;
        xk=xk+a*d;
        x = [x xk];
        dfk = df(xk);
        it=it+1;
    end    
end