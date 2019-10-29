clear all;
close all;

s = tf('s');
p = 3;
z = 2;
% G = [1/(s-p) 0;0 1/(s+3)] *[sqrt(3)/2 -1/2;1/2 sqrt(3)/2] *[(s-z)/(0.1*s+1) 0;0 (s+2)/(0.1*s+1)];
% G = (s-1)*(s-3)/(s-2)/(s+1)^2
G = 1/(0.2*s+1)/(s+1)*[1 1;1+2*s 2];


[ptot,ztot] = pzmap(G);
p = ptot(find(ptot>0));
z = ztot(find(ztot>0));
np = length(p);
nz = length(z);
G = ss(G);
[V,E]=eig(G.A);
C = G.C*V  % output pole vectors
for ii = 1:np
    Yp(:,ii) = C(:,ii)/norm(C(:,ii))   %pole direction
end
for ii = 1:nz
    [U,S,V] = svd(evalfr(G,z(ii)));
    Yz(:,ii) = U(:,end)
end

Qp = (Yp'*Yp).*(1./(diag(p')*ones(np)+ones(np)*diag(p)))
Qz = (Yz'*Yz).*(1./(diag(z)*ones(nz)+ones(nz)*diag(z')))
Qzp = (Yz'*Yp).*(1./(diag(z)*ones(nz,np)-ones(nz,np)*diag(p)))
Msmin = sqrt(1+norm(sqrtm(inv(Qz))*Qzp*sqrtm(inv(Qp)))^2)

pha = acos(Yz'*Yp)
Msmin2 = sqrt(sin(pha)^2 + abs(z+p)^2/abs(z-p)^2*cos(pha)^2)
