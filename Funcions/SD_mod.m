%{
B: Hessiana modificada, no té cap utilitat en l'algoritme, pero a cada pas
la guardem a la Hk
B1 = B^-1
%}
function [B,B1, mod] = SD_mod(H,delta)
    [Q,A] = eig(H, 'vector');
    Am = max(A, delta);
    mod = min(A == Am); % mirem si hem canviat algun eigenvalue 
    B = Q*(diag(Am))*Q'; % H simetrica => Q^-1 = Q'
    B1 = Q*diag(1./Am)*Q';
end