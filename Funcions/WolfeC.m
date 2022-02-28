%{
Input:
    - x: punt actual
    - a: alpha actual
    - d: direcció de descens
    - c1, c2: constants per WC
    - iW: metode Wolfe a usar
Output:
    - b: cert si cumpleix totes les WC especificades per iW
    - iWout: indicador de les condicions que compleix ak
%}
function [b, iWout] = WolfeC(x, a, f, df, d, c1, c2, iW)
    iWout = 0;
    b = false;
    if f(x+a.*d) <= f(x) + c1.*df(x)'*d*a
        iWout = 1;
        if iW == 1
            if df(x+a.*d)'*d >= c2.*df(x)'*d 
                iWout = 2;
                b = true;
            end
        else
            if abs(df(x+a.*d)'*d) <= c2*abs(df(x)'*d)
                iWout = 3;
                b = true;
            end
        end
    end
end