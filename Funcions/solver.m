% isc = 1 GM, ISD = 2 CGM, ISD = 3 BFGS
% icd:  variant de CGM icd = 1 FR, icd = 2 PR+
% irc: Restart per la CGM, irc = no, irc = 1 si
% nu: nombre de iteracions entre restarts en cas de que irc = 1
function [xk, dk, ak, iWk, betak, it] = solver(x, f, df, d2f, amin, amax, p, c1, c2, iW, tol, itmax, isd, icd, irc, nu)
    betak = [];
    if isd == 1
       [xk, dk, ak, iWk, it] = GM(x, f, df, amin, amax, p, c1, c2, iW, tol, itmax);
    else
        
    end
    
    while norm(dfxk) > tol & it < itmax
        [a, iWout] = BLS(x,f,df,d,amin,amax, p, c1, c2, iW);
        iWk = [iWk iWout];
        ak = [ak a];
        x = x + a*d;
        xk = [xk x]; 

        if isd == 2
            if icd == 1
               B  
            end
            B = 
        end
        dfxk = -df(x);
        if isc == 1 | restart == nu
            d = -dfxk;
        else
            B = CGM(dfk, )
        end
        %{
        
        % Si no troba alpha, cal fer més exhaustiva la cerca
        if iWout < 2
            p = 1.1*p;
        %}
        it = it + 1;
        dfxk = df(x);
    end

end