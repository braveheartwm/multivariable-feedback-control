clear all;
close all;
clc;
% LQG controller design.

% plant
G = tf(3*[-2 1],conv([5 1],[10 1]));
[a,b,c,d] = ssdata(G);

% Model dimensions.
p = size(c,1);

[n,m] = size(b);
Znm = zeros(n,m);
Zmm = zeros(m,m);
Znn = zeros(n,n);
Zmn = zeros(m,n);

% 1)design state feedback regulator
A = [a,Znm;-c Zmm];
B = [b;-d];
Q = [Znn Znm;Zmn eye(m,m)];
R = eye(m);
Kr = lqr(A,B,Q,R);
Krp = Kr(1:m,1:n);
Kri = Kr(1:m,n+1:n+m);

% 2)design Kalman filter
Bnoise = eye(n);
W = eye(n); V = eye(m);
Estss = ss(a,[b Bnoise],c,[0 0 0]);
[Kess,Ke] = kalman(Estss,W,V);

% 3)Form overall controller
Ac = [Zmm Zmn;-b*Kri a-b*Krp-Ke*c];
Bcr = [eye(m);Znm];
Bcy = [-eye(m);Ke];
Cc = [-Kri -Krp];
Dcr = Zmm;
Dcy = Zmm;
K1qg2 = ss(Ac,[Bcr Bcy],Cc,[Dcr, Dcy]);
K1qg = ss(Ac,-Bcy,Cc,-Dcy);

% simulation
sys1 = feedback(G*K1qg,1);
step(sys1,50);
sys = feedback(G*K1qg2,1,2,1,+1);
sys2 = sys*[1;0];hold;
step(sys2,50);

% [Gm,Pm,Wgm,Wpm] = margin(G*K1qg2) 

