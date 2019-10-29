% Calculates bound on peak of S and T for MIMO systems
%    
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: STboundMIMO.m,v 1.0 2004/10/31 19:30:00 kariwala Exp $

% G: must have distinct and at least one RHP-zero and one RHP-pole

A = [1 -0.5;0.5 1]; 
B = [0.4 0.8;0.15 0.2];
C = [0.5 0.15;0.3 0.9];
G = ss(A,B,C,eye(2));

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

Qp = Yp'*Yp.*(1./(diag(p')*ones(np) + ones(np)*diag(p)));
Qz = Yz'*Yz.*(1./(diag(z)*ones(nz) + ones(nz)*diag(z')));
Qpz = conj(Yp'*Yz).*(1./(-diag(p)*ones(np,nz) + ones(np,nz)*diag(z)));

Msmin = sqrt(1+norm(sqrtm(inv(Qp))*Qpz*sqrtm(inv(Qz)))^2)

Wu = 0.0001*eye(2);

% systemnames = 'G Wu';
% inputvar = '[r(2); u(2)]';
% outputvar = '[r-G; Wu; r-G]';
% input_to_G = '[u]';
% input_to_Wu = '[u]';
% sysoutname = 'P';
% cleanupsysic = 'yes';
% sysic;
% 
% [K,CL,gamma] = hinfsyn(P,2,2); gamma
[K,CL,gamma] = mixsyn(G,[],Wu,[]); gamma

% Pole vectors using state space realization (just to confirm)
%[V,E] = eig(G.A);
%C1 = G.C*inv(V);

