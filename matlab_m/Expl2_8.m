%Example 2.8, Loop-shaping design for plant (2.31)
%
% Plant (2.31), G(s)=3(-2s+1)/(5s+1)(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_8.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

% plant G
clear
s=tf('s');
G=3*(-2*s+1)/(5*s+1)/(10*s+1);

% frequency response of G
points=301;
omega=logspace(-3,1,points);

%Loopshaping controller K
kc=0.05;tauf=0.33;
K=kc*(10*s+1)*(5*s+1)/s/(2*s+1)/(0.33*s+1);
[magK,phaK]=bode(K,omega); 

% ANALYSIS OF CONTROLLERS

L=minreal(G*K);
S = 1/(1+L);
T = 1-S;
time=0:0.1:50;
u=ones(size(time));
y1=lsim(T,u,time);

%input signal
KS=K*S; u1=lsim(KS,u,time);

%Figure 2.19 Bode plot L
figure(1) 
clf
margin(L);  
h=get(gcf,'children');axes(h(2));
axis([0.01 10 0.01 20])
mS=norm(S,inf,1e-4);
mT=norm(T,inf,1e-4);
disp(sprintf('[Ms,Mt]=[%5.2f,%5.2f ]',mS(1),mT(1)));

%Figure 2.20 Step response
figure(2)
clf
plot(time,y1,time,u1,time,u,':')
axis([0 50 -0.5 2.5]);
text(40,1.2,'OUTPUT (y)');text(40,0.1,'INPUT (u)');
xlabel('Time');%ylabel('y');
title('Figure 2.17');

%Figure 2.21 Bode plot K
figure(3)
clf
loglog(omega,magK(:),'b',omega,omega./omega,':')
ylabel('Magnitude')
xlabel('Frequency')
title('Figure 2.18');
axis([.001,10,.1,100])



