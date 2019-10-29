% SISO system
s = tf('s');
G = ss((s-2)/(s^2 - 2*s + 5));
[p,z] = pzmap(G);

[V,E] = eig(G.A);
C1 = G.C*V;

yp1 = C1(1)/norm(C1(1));  % pole directions (ref: page 133, Skogestad and Postlethwaite)
yp2 = C1(2)/norm(C1(2));

yz = 1; % zero direction

Qz = 1/(z+z');
Qp = [yp1'*yp1/(p(1)'+p(1)) yp1'*yp2/(p(1)'+p(2));yp2'*yp1/(p(2)'+p(1)) yp2'*yp2/(p(2)'+p(2))];
Qzp = [yz'*yp1/(z-p(1)) yz'*yp2/(z-p(2))];

Msmin = sqrt(1+norm(sqrtm(inv(Qz))*Qzp*sqrtm(inv(Qp)))^2)

% MIMO system

clear all
A = [1 -0.5;0.5 1]; 
B = [0.4 0.8;0.15 0.2];
C = [0.5 0.15;0.3 0.9];
G = ss(A,B,C,eye(2));
[p,z] = pzmap(G);

[V,E] = eig(G.A);
C1 = G.C*V;

yp1 = C1(:,1)/norm(C1(:,1));  % pole directions (ref: page 133, Skogestad and Postlethwaite)
yp2 = C1(:,2)/norm(C1(:,2));

[U,S,V] = svd(evalfr(G,z(1))); yz1 = U(:,end); %zero directions
[U,S,V] = svd(evalfr(G,z(2))); yz2 = U(:,end); %zero directions

Qz = [yz1'*yz1/(z(1)+z(1)') yz1'*yz2/(z(1)+ z(2)');yz2'*yz1/(z(2)+z(1)') yz2'*yz2/(z(2)+z(2)')];
Qp = [yp1'*yp1/(p(1)'+p(1)) yp1'*yp2/(p(1)'+p(2));yp2'*yp1/(p(2)'+p(1)) yp2'*yp2/(p(2)'+p(2))];
Qzp = [yz1'*yp1/(z(1)-p(1)) yz1'*yp2/(z(1)-p(2)); yz2'*yp1/(z(2)-p(1)) yz2'*yp2/(z(2)-p(2))];

Msmin = sqrt(1+norm(sqrtm(inv(Qz))*Qzp*sqrtm(inv(Qp)))^2)

