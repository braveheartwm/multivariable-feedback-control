%Example 2.11, Two degree of freedom controller.
%Produce Figure 2.26 

%Plant (2.62):  G(s)=200/(10s+1)(0.005s+1)^2,  Gd(s)=100/(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_9.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

clear
s=tf('s');
Eq2_62

%K3: PID control
K3=0.5*(s+2)*(0.05*s+1)/(s+1e-4)/(0.005*s+1);
L3=G*K3;
T3=minreal(L3/(1+L3));

%Kr: Prefilter for the reference:
Kr=(0.5*s+1)/(0.65*s+1)/(0.03*s+1);

T32 = T3*Kr;  % get Hinf-norm of t32 equal to 1.0
time=0:0.01:3;k=length(time);
y3r=step(T3,time);
y32r=step(T32,time);
Ov3=max([y3r;1]);[x,x,x,Tr3]=margin(y3r(2:k)/0.9,y3r(2:k),time(2:k));
Ov32=max([y32r;1]);[x,x,x,Tr32]=margin(y32r(2:k)/0.9,y32r(2:k),time(2:k));
disp(sprintf('[Ov3 tr3 Ov32 tr32]=[%5.2f%5.2f%5.2f%5.2f ]',Ov3,Tr3,Ov32,Tr32));

figure(1);clf;subplot(211);
plot(time,[y32r,y3r],[0,3],[1,1],':')
axis([0,3,-0.2,1.5]);
xlabel('Time');ylabel('y');
text(0.3,1.38,'y3');text(0.5,0.85,'y3 (two degrees-of-freedom)');
