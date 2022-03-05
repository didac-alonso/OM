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
function [xk, dk, ak, Bk, iWk, it] = CGM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax, isd, icd, irc, nu)
    it = 0;
    pp = p;
    xk = [x];
    Bk = [];
    dfxk = df(x);
    dk = [-dfxk];
    d = -dfxk;
    ak = [];
    iWk = [];
    while norm(dfxk) > tol & it < itmax
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW);
        iWk = [iWk iWout];
        ak = [ak a];
        x = x + a*d;
        xk = [xk x]; 

        if isd == 2
            if icd == 1
               B  
            end
            B = 
        end
        dfxk = -df(x);
        if isc == 1 | restart == nu
            d = -dfxk;
        else
            B = CGM(dfk, )
        end
        %{
        
        % Si no troba alpha, cal fer més exhaustiva la cerca
        if iWout < 2
            p = 1.1*p;
        %}
        it = it + 1;
        dfxk = df(x);
    end

end