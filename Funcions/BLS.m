function [ak, iWout] = BLS(xk, f, df, dk, amax, amin, p, c1, c2, iW)
    ak = amax;
    iWout = 0;
    [b, iWout] = WolfeC(xk, ak, f, df, dk, c1, c2, iW);
    while ak > amin & ~b
        ak = p*ak;
        [b, iWout] = WolfeC(xk, ak, f, df, dk, c1, c2, iW);
    end
end