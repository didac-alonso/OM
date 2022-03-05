%{
Busca els mínims d'una funció quadràtica a partir de la fòrmula exacta per les alphas.
La funció ha de ser de la forma f = 1/2x'Qx+b'x+k

Input:
    - x: Punt de partida
    - df: Derivada de la funció a minimitzar
    - Q: Segona derivada de f
    - b: Paràmetre de la quadràtica

Output:
    - x: Punt on hi ha un mínim
    - it: Nombre d'iteracions usades
%}

function [xk, dk, ak, iWk, it] = GM_ELS(x, df, Q, tol, itmax)
    it = 1;  
    dfx = df(x);
    d = -dfx;
    dk = [d];
    xk = [x];
    ak = [];
    while norm(dfx) > tol & it <= itmax
        a = -(dfx'*d)/(d'*Q*d);
        ak = [ak a];
        x = x + a*d;
        xk = [xk x];
        dfx = df(x);
        d = -dfx;
        dk = [d dk];
        it = it + 1;
    end
    iWk = zeros(it, 1);
end