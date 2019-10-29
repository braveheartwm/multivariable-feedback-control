%Example 2.9 produce Figure 2.22
%
%Plant (2.62):  G(s)=200/(10s+1)(0.005s+1)^2,  Gd(s)=100/(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_9.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

% plant
clear
s=tf('s');
Eq2_62

% Controller for command following (which inverts the plant)
wc0=10;
K=wc0*(10*s+1)*(0.1*s+1)/s/200/(0.01*s+1);
L = G*K;
S = 1/(1+L);
T=1-S;
SGd = S*Gd;
time=0:0.01:3;
u=ones(size(time));
y0=lsim(SGd,u,time);
y0r=lsim(T,u,time); % rise time 0.16s, no overshoot

figure(1);
subplot(121)
plot(time,y0r,time,u,':')
axis([0,3,-0.2,1.5]);
xlabel('Time');ylabel('y');
title('TRACKING RESPONSE');
subplot(122)
plot(time,y0,time,time*0,':')
axis([0,3,-0.2,1.5]);
xlabel('Time');ylabel('y');
title('DISTURBANCE RESPONSE');
