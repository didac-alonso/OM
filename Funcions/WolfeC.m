%{
Input:
    - xk: punt actual
    - ak: alpha actual
    - dk: direcció de descens
    - c1, c2: constants per WC
    - iW: metode Wolfe a usar
Output:
    - b: cert si cumpleix totes les WC especificades per iW
    - iWout: indicador de les condicions que compleix ak
%}
function [b, iWout] = WolfeC(xk, ak, f, df, dk, c1, c2, iW)
    i = f(xk+ak.*dk) <= f(xk) + c1.*df(xk)'*dk*ak;
    iWout = 0;
    if i 
        iWout = 1;
    end
    if iW == 1
        b = i & df(xk+ak.*dk)'*dk >= c2.*df(xk)'*dk;
        if b 
            iWout = 2;
        end
    else
        b = i & abs(df(xk+ak.*dk)'*dk) <= c2.*abs(df(xk)'*dk);
        if b
            iWout = 3;
        end
    end
end