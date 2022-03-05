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
    - icd:  variant de CGM icd = 1 FR, icd = 2 PR+
    - irc: Restart per la CGM, irc = no, irc = 1 si
    - nu: nombre de iteracions entre restarts en cas de que irc = 1

Output:
    - x: Punt que fa 0 df.
    - it: nombre d'iteracions usades
%}
function [xk, dk, ak, Bk, iWk, it] = CGM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax, icd, irc, nu)
    it = 0;
    xk = [x];
    Bk = [];
    dfx = df(x);
    dk = [-dfx];
    d = -dfx;
    ak = [];
    iWk = [];
    restart = 0;
    if irc == 0
        nu = 3;
    end
    while norm(dfx) > tol & it < itmax
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW);
        iWk = [iWk iWout];
        ak = [ak a];
        x = x + a*d;
        xk = [xk x];
        dfx1 = dfx;
        dfx = df(x);
        % Si irc == 0 no augmentarà
        restart = restart + irc;
        if restart == nu
            restart = 0;
            B = 0;
        elseif icd == 1
            B = dfx'*dfx/(norm(dfx1))^2;
        else
            B = max(0, (dfx'*(dfx-dfx1))/(norm(dfx1))^2);   
        end
        Bk = [Bk B];
        d = -dfx + B*dk(end);
        dk = [dk d];
        it = it + 1;
    end
end