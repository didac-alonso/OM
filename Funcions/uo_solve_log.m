function [gk,la1k,kappak,rk,Mk] = uo_solve_log(x1,f,g,h,epsG,kmax,almax,almin,rho,c1,c2,iW,isd,icg,irc,nu,delta,xk,dk,alk,iWk,betak,Hk,tauk,xo,xylim,logfreq)
diary ('uo_solve_log.out'); diary on;
niter = size(xk,2); n= size(x1,1); nxo = size(xo,1);
if nxo==0 xo = xk(:,niter); end
fk  = []; gk = []; rk  = []; gdk = []; rk=[]; Mk=[];
la1k=zeros(1,niter);     % if NM or MNM: lowest vap of h(xk:,k)); if QNM: lowest vap of Hk(:,:,k).
kappak = zeros(1,niter); % if NM : cond. number of h(xk(:,k))), if +def; if QNM or MNM: cond. number of Hk or Bk resp.
for k = 1:niter
    x = xk(:,k); fk = [fk,f(x)]; gk = [gk,g(x)]; 
    if isd >3 la1k(k) = min(eig(h(x))); end
    if isd == 4 & all(eig(h(x))>0) kappak(k)= cond(h(xk(:,k))); end
    if k < niter
        if isd < 4
            rk = [rk, (f(xk(:,k+1))-f(xo))/(f(xk(:,k))-f(xo))   ];
            Mk = [Mk, (f(xk(:,k+1))-f(xo))/(f(xk(:,k))-f(xo))^2 ];
        else
            rk = [ rk, norm(g(xk(:,k+1))) / norm(g(xk(:,k)))   ];
            Mk = [ Mk, norm(g(xk(:,k+1))) / norm(g(xk(:,k)))^2 ];
        end
        gdk = [gdk,gk(:,k)'*dk(:,k)];
        if     isd == 3 la1k(k) = min(eig(Hk(:,:,k))); end
        if  isd >= 3 & isd ~=4  kappak(k) = cond(Hk(:,:,k)); end
    end
end
if niter > 1
    if (isd  < 4 & nxo==0) rk(niter-1)=rk(niter-2); Mk(niter-1)= Mk(niter-2); end
    if isd <= 4  tauk(1:niter-1) = 0; end
    if isd == 5  tauk(1:niter-1) = delta; end
end
if isd == 3  la1k(niter)     = 0; end

fprintf('   [uo_solve_log]\n');
fprintf('   f= %s\n', func2str(f));
fprintf('   epsG= %3.1e, kmax= %4d\n', epsG,kmax);
if isd ~= 4 fprintf('   almax= %2d, almin= %3.1e, rho= %4.2f, c1= %3.2f, c2= %3.2f, iW= %1d\n',almax,almin,rho,c1,c2,iW); end
fprintf('   isd= %1d\n',isd);
if isd == 2 fprintf('   icg= %1d, irc= %1d, nu= %3.2f\n',icg,irc,nu); end
if isd == 5 fprintf('   delta= %3.1d\n',delta); end
if isd < 4
    fprintf('   ek = fk-f*\n');
else
    fprintf('   ek = ||gk||\n');
end
    if n==2 fprintf('   x1 = [ %+3.1e , %+3.1e ]. It k: x(k+1) <- xk+alk*dk\n', x1(1), x1(2)); end
if isd ==1
    fprintf('      k     g''*d        al iW    ||g||        f         r        M\n');
elseif isd ==2
    fprintf('      k     g''*d        al iW     beta    ||g||        f         r        M\n');
elseif isd == 3
    fprintf('      k     g''*d        al iW    la(1)    kappa    ||g||        f         r        M\n');   
else
    fprintf('      k     g''*d        al iW    la(1) del./tau    kappa    ||g||        f         r        M\n');   
end
if niter == 1
    krange=[];
elseif niter == 2
        krange =[1];
else
        krange=[1:logfreq:max(2,niter-11),max(3,niter-10):niter-1];
end
for k = krange
%for k = [1:logfreq:max(2,niter-11),max(3,niter-10):niter-1]
    if isd == 1 
        fprintf(' %6d %+3.1e %+3.2e  %1d %+3.1e %+3.1e %+3.2e %+3.1e\n', k, gdk(k), alk(k), iWk(k), norm(gk(:,k)), fk(k), rk(k),Mk(k));
    elseif isd == 2
        fprintf(' %6d %+3.1e %+3.2e  %1d %+3.1e %+3.1e %+3.1e %+3.2e %+3.1e\n', k, gdk(k), alk(k), iWk(k),  betak(k), norm(gk(:,k)), fk(k), rk(k),Mk(k));
    elseif isd == 3
        fprintf(' %6d %+3.1e %+3.2e  %1d %+3.1e %+3.1e %+3.1e %+3.1e %+3.2e %+3.1e\n', k, gdk(k), alk(k), iWk(k), la1k(k), kappak(k), norm(gk(:,k)), fk(k), rk(k),Mk(k));      
    else        
        fprintf(' %6d %+3.1e %+3.2e  %1d %+3.1e %+3.1e %+3.1e %+3.1e %+3.1e %+3.2e %+3.1e\n', k, gdk(k), alk(k), iWk(k), la1k(k), tauk(k), kappak(k), norm(gk(:,k)), fk(k), rk(k),Mk(k));      
    end
end
if isd == 1
%    fprintf('      *                      %+3.1e %+3.1e\n', norm(gk(:,niter)), fk(niter));
    fprintf(' %6d                       %+3.1e %+3.1e\n', niter, norm(gk(:,niter)), fk(niter));
    fprintf('      k     g''*d        al iW    ||g||        f         r        M\n');
elseif isd ==2
%    fprintf('      *                               %+3.1e %+3.1e\n', norm(gk(:,niter)), fk(niter));
    fprintf(' %6d                                %+3.1e %+3.1e\n', niter, norm(gk(:,niter)), fk(niter));
    fprintf('      k     g''*d        al iW     beta    ||g||        f         r        M\n');
elseif isd == 3
%    fprintf('      *                                        %+3.1e %+3.1e\n', norm(gk(:,niter)), fk(niter));
    fprintf(' %6d                                         %+3.1e %+3.1e\n', niter, norm(gk(:,niter)), fk(niter));
    fprintf('      k     g''*d        al iW    la(1)    kappa    ||g||        f         r        M\n');
else
%    fprintf('      *                      %+3.1e                   %+3.1e %+3.1e\n', la1k(niter), norm(gk(:,niter)), fk(niter));
    fprintf(' %6d                       %+3.1e                   %+3.1e %+3.1e\n', niter, la1k(niter), norm(gk(:,niter)), fk(niter));
    fprintf('      k     g''*d        al iW    la(1) del./tau    kappa    ||g||        f         r        M\n');
end

if n==2
    if isd==1 & size(h(xk(:,niter)),1)>0
        fprintf('   x* = [ %+3.1e , %+3.1e ]; rUB = %+3.2e \n', xk(1,niter), xk(2,niter),(max(eig(h(xk(:,niter))))-min(eig(h(xk(:,niter)))))^2/(max(eig(h(xk(:,niter))))+min(eig(h(xk(:,niter)))))^2);
    else fprintf('   x* = [ %+3.1e , %+3.1e ]\n', xk(1,niter), xk(2,niter));
    end
end
fprintf('   [uo_solve_log]\n');
fs = 0;
iplot = 1;
if iplot == 1 & niter > 1
    if n == 2
        if size(xylim) == [0 0] xylim = [0 0 0 0]; end
        subplot(2,2,1); uo_solve_plot(f, xk, gk, xylim, 1, fs); subplot(2,2,2); uo_solve_plot(f, xk, gk, xylim, 2,fs);
        subplot(2,2,3);plot(rk(1:niter-1),'-o');xlabel('Iterations k');title('r^k'); subplot(2,2,4); plot(Mk(1:niter-1),'-x');xlabel('Iterations k'); title('M^k');
    else
        subplot(1,2,1);plot(rk(1:niter-1),'-o');xlabel('Iterations k');title('r^k'); subplot(1,2,2); plot(Mk(1:niter-1),'-x');xlabel('Iterations k'); title('M^k');
    end
end
diary off;
end

