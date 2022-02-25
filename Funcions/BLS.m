%{
Input:
    - xk: punt actual
    - dk: direcció de descens
    - c1, c2: per Wolfe
    - p: rho, com fem descendir alpha
    - a: alpha, donada per interval amax, amin
    - iW: metode Wolfe a usar
Output:
    - ak: alpha òptima
    -iWout: indicador de les condicions que compleix ak
%}
function [ak, iWout] = BLS(xk, f, df, dk, amax, amin, p, c1, c2, iW)
    ak = amax;
    iWout = 0;
    [b, iWout] = WolfeC(xk, ak, f, df, dk, c1, c2, iW);
    while ak > amin & ~b
        ak = p*ak;
        [b, iWout] = WolfeC(xk, ak, f, df, dk, c1, c2, iW);
    end
    if iWout < 2
        disp('No solution found');
    end
end