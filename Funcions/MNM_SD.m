%{
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
function [xk, dk, alk, Hk, iWk, it]= MNM_SD(x, f, df, d2f, amin, amax, p, c1, c2, iW, tol,itmax, delta)
    it=1;
    xk = [x];
    dfk = df(x);
    Hk = [];
    dk = [];
    alk = [];
    iWk = [];
    while norm(dfk)>tol & it<=itmax
        [B,B1, mod] = SD_mod(d2f(x),delta);
        Hk = cat(3, Hk, B);
        d = -B1*dfk;
        dk = [dk d];
        %if mod
        [a, iWout] = BLS(x, f, df, d, amin, amax, p, c1, c2, iW, []);
        %else
       %     a = 1; iWout = 4;
        %end
        alk = [alk a];
        iWk = [iWk iWout];
        x=x+a*d;
        xk = [xk x];
        dfk = df(x);
        it=it+1;
    end    
end