% Figure 2.1    Sinusoidal response 
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Fig2_1.m,v 1.2 2004/01/26 16:37:11 heidisi Exp $

clear
num=5;den=[10 1];delay=2;w=0.2;
G=tf(num,den,'OutputDelay',delay); %Input- and OutputDelay equiv. in SISO
t=0:.5:100;um=1;up=0;
u=um*sin(w*t+up);   %Input
[mag,pha]=bode(G,w);pha=pha*pi/180; %alt. deg2rad in mapping toolbox
ym=um*mag;yp=up+pha;
y=ym*sin(w*t+yp);   %Output
clf
axis([0 100 -3 3])
hold;
plot(t,u,t,y,t,t-t,':');
xlabel('TIME [s]');
text(30,1,'u(t)'); text(37,2,'y(t)');
