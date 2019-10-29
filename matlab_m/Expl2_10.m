%Example 2.10, Loop-shaping design for Plant (2.62)
%Produce Figure 2.23 and Table 2.3
%

%Plant (2.62):  G(s)=200/(10s+1)(0.005s+1)^2,  Gd(s)=100/(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_10.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

% plant
clear
s=tf('s');
Eq2_62

% Controllers for disturbance rejection, designed by loop-shaping
wc0=10;
K0=wc0*(10*s+1)*(0.1*s+1)/s/200/(0.01*s+1);
K1 = 0.5;
K2=0.5*(s+2)/s;
K3=K2*(0.05*s+1)/(0.005*s+1);
omega=logspace(-2,2,41);

%Table 2.2
disp(sprintf('%6s%6s%6s%6s%6s%6s%11s%13s','','','','','','',...
    'Reference','Disturbance'));
disp(sprintf('%4s%6s%6s%6s%6s%6s%6s%6s%7s%7s','','GM','PM',...
    'Wc','Ms','Mt','tr','Ov','Ymax','Y(3)'));

% 0. Inverse -based controller.
K=K0;
sysanaly;
disp(sprintf('%4s%7.2f%6.1f%6.1f%6.2f%6.2f%6.2f%6.2f%6.2f%7.2f',...
    'K0',Gm,Pm,Wc,MS,MT,tr,Ov,Ymax,Y3));

% 1. P controller
K = K1;
sysanaly;
[mag1,pha1]=bode(L,omega);
y1=yd;
disp(sprintf('%4s%7.2f%6.1f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%7.2f',...
    'K1',Gm,Pm,Wc,MS,MT,tr,Ov,Ymax,Y3));

% 2. PI controller
K = K2;
sysanaly;
[mag2,pha2]=bode(L,omega);
y2=yd;
disp(sprintf('%4s%7.2f%6.1f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%7.3f',...
    'K2',Gm,Pm,Wc,MS,MT,tr,Ov,Ymax,Y3));

% 3. PID controller 
K = K3;
sysanaly;
[mag3,pha3]=bode(L,omega);
y3=yd;
disp(sprintf('%4s%7.1f%6.1f%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f%7.3f',...
    'K3',Gm,Pm,Wc,MS,MT,tr,Ov,Ymax,Y3));

%Figure 2.24
clf;
subplot(221)
loglog(omega,mag1(:),omega,mag2(:),omega,mag3(:),omega,omega./omega,':')
axis([.01,100,.01,10000]);
text(0.02,30,'L1');text(0.18,1000,'L2, L3');
text(50,0.15,'L3');text(15,0.05,'L1, L2');
title('Loop gains');
subplot(222)
plot(time,y3,time,y1,time,y2,time,time*0,':')
axis([0,3,-0.2,1.5]);
xlabel('Time');ylabel('y');
text(0.3,0.2,'y2');text(0.53,0.5,'y3');text(1.5,0.85,'y1');
title('Disturbance responses');


