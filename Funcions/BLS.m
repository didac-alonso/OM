%{
Input:
    - x: punt actual
    - d: direcció de descens
    - c1, c2: per Wolfe
    - p: rho, com fem descendir alpha
    - a: alpha, donada per interval amax, amin
    - iW: metode Wolfe a usar
        - 1: WC
        - 2: SWC
Output:
    - ak: alpha òptima
    -iWout: indicador de les condicions que compleix ak
%}
function [a, iWout] = BLS(x, f, df, d, amin, amax, p, c1, c2, iW, Q)
    if iW == 0
       a = -(df(x)'*d)/(d'*Q*d);
       iWout = 0;
    else
        a = amax;
        iWout = 0;
        [b, iWout] = WolfeC(x, a, f, df, d, c1, c2, iW);
        while a >= amin & ~b
            a = p*a;
            [b, iWout] = WolfeC(x, a, f, df, d, c1, c2, iW);
        end
        if iWout < 2
            fprintf('No solution found, iWout = %d\n', iWout);
        end
    end
end