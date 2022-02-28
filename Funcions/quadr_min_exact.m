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

function [x, it] = quadr_min_exact(x, df, Q, b, tol, itmax)
    it = 0;    
    while norm(df(x)) > tol & it < itmax
        d = -df(x);
        a = -((Q*x-b)'*d)/(d'*Q*d);
        x = x + a*d;
        it = it + 1;
    end
end