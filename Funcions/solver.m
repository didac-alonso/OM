% isd = 1 GM, ISD = 2 CGM, ISD = 3 BFGS
% icd:  variant de CGM icd = 1 FR, icd = 2 PR+
% irc: Restart per la CGM, irc = no, irc = 1 RC1, irc = 2 RC2
% nu: nombre de iteracions entre restarts en cas de que irc = 1
function [xk, dk, ak, iWk, betak, it] = solver(x, f, df, d2f,tol, itmax, amax, amin, p, c1, c2, iW, isd, icg, irc, nu)
    betak = [];
    if isd == 1
       [xk, dk, ak, iWk, it] = GM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax, d2f(0));
    elseif isd == 2
        [xk, dk, ak, betak, iWk, it] = CGM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax, icg, irc, nu, d2f(0));
    end
end