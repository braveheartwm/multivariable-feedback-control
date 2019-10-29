% Calculates bound on peak of S and T for MIMO systems
%    
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: STboundMIMO.m,v 1.0 2004/10/31 19:30:00 kariwala Exp $

% G: must have distinct and at least one RHP-zero and one RHP-pole

[ptot,ztot] = pzmap(G);                % poles and zeros
p = ptot(find(ptot>0)); z = ztot(find(ztot>0));  % RHP poles and zeros
np = length(p); nz = length(z);

G = ss(G);
[V,E] = eig(G.A); C = G.C*V; % output pole vectors

for i = 1:np
    Yp(:,i) = C(:,i)/norm(C(:,i)); %normalize columns to get pole directions
end
    
for i = 1:nz
    [U,S,V] = svd(evalfr(G,z(i))); Yz(:,i) = U(:,end); %zero directions
end

Qp = (Yp'*Yp).*(1./(diag(p')*ones(np) + ones(np)*diag(p)));
Qz = (Yz'*Yz).*(1./(diag(z)*ones(nz) + ones(nz)*diag(z')));
Qzp = (Yz'*Yp).*(1./(diag(z)*ones(nz,np) - ones(nz,np)*diag(p)));

Msmin = sqrt(1+norm(sqrtm(inv(Qz))*Qzp*sqrtm(inv(Qp)))^2)
