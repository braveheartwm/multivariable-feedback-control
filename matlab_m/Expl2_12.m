%Example 2.12, Flexible loop-shaping design.
%Produce Figure 2.27
%
%Plant (2.75): G(s)=Gd(s)=2.5s(s^2+1)/(s^2+0.5^2)(s^2+2^)
%
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_12.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

clear
omega=logspace(-2,3,410);
s=tf('s');
G=2.5*s*(s^2+1)/(s^2+0.5^2)/(s^2+2^2);
Gd=G;

%Controller (1.87): K(s)=1
K=1;
L=G*K;
S=1/(1+L);
SGd=S*Gd;
[mag,pha]=bode(L,omega);
time=0:0.1:20;
y1=step(SGd,time);
yol=step(Gd,time);

subplot(221)
loglog(omega,mag(:),'b',omega,omega./omega,':')
axis([.01,100,.01,100]);
xlabel('Frequency');ylabel('Magnitude');text(10,10, 'G=Gd')
subplot(222)
plot(time,y1,time,yol,'--',time,time*0,':')
xlabel('Time');ylabel('y(t)');
text(5,1.5,'YOL');text(13,-0.5,'YCL');
text(13,-0.85,'(solid line)')

