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