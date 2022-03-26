% isd = 1 GM, ISD = 2 CGM, ISD = 3 BFGS, 4 NM, 5 MNM-SD, 6 MNM-CMI
% icg:  variant de CGM icd = 1 FR, icd = 2 PR+
% irc: Restart per la CGM, irc = no, irc = 1 RC1, irc = 2 RC2
% nu: nombre de iteracions entre restarts en cas de que irc = 1
function [xk, dk, alk, iWk, betak, Hk, tauk, it] = solver(x, f, g, h,epsG, kmax, almax, almin, rho, c1, c2, iW, isd, icg, irc, nu, delta)
    betak = [];
    Hk = [];
    tauk = [];
    if isd == 1
        [xk, dk, alk, iWk, it] = GM(x, f, g, almin, almax, rho, c1, c2, iW, epsG, kmax, h(zeros(size(x,1))));
    elseif isd == 2
        [xk, dk, alk, betak, iWk, it] = CGM(x, f, g, almin, almax, rho, c1, c2, iW, epsG, kmax, icg, irc, nu, h(zeros(size(x,1))));
    elseif isd == 3
        [xk, dk, alk, Hk, iWk, it] = BFGS(x, f, g, almin, almax, rho, c1, c2, iW, epsG, kmax, h(zeros(size(x,1))));
    elseif isd == 4
        [xk, dk, alk, Hk, iWk, it] = Newton(x, f, g, h, almin, almax, rho, c1, c2, iW, epsG,kmax);
    elseif isd == 5
        [xk, dk, alk, Hk, iWk, it] = MNM_SD(x, f, g, h, almin, almax, rho, c1, c2, iW, epsG,kmax, delta);
    elseif isd == 6
        [xk, dk, alk, Hk, iWk, tauk, it]= MNM_CMI(x, f, g, h, almin, almax, rho, c1, c2, iW, epsG,kmax);
    end
end