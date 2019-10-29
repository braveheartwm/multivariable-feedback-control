% Figure 5.4: Effect of increased controller gain on |S| for system with
% RHP-zero at z = 2, L(s)= k/s * 2-s/2+s
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig5_4.m,v 1.4 2004/02/03 14:12:00 vidaral Exp $
% This example illustrates the trade-off between bandwith requirements
% and the peak of the sensitivity function S. The bandwith is 
% increased by increasing the controller gain, k, for a plant with a RHP 
% zero at z=2. It is demonstrated that this results in an increase in the 
% peak of s (infinity norm)

clear all; close all;
% set(cstprefs.tbxprefs,'MagnitudeUnits','log','MagnitudeScale','log');
% constants
k1=0.1;
k2=0.5;
k3=1;
k4=2;
w=logspace(-3,1,301);

% plant g(s) = (2-s)/(2+s)
g=tf([-1 2],[1 2]);

% minimal realization of loop transfer function with controller k(s) = k/s
L1=minreal(g*tf(k1,[1 1e-10]));
L2=minreal(g*tf(k2,[1 1e-10]));
L3=minreal(g*tf(k3,[1 1e-10]));
L4=minreal(g*tf(k4,[1 1e-10]));
% frequency response of sensitivity function: S=1/(1+L)
S1=1/(1+L1); S2=1/(1+L2); S3=1/(1+L3); S4=1/(1+L4);
bodemag(S1,'-',S2,'-',S3,'-',S4,'-',tf(1),'-.',w);
 
xlabel('Frequency [rad/s]')
ylabel('Magnitude |S|')
axis([0.0100 10.0000 0.1000 30.0000]);
text(0.12,1.2,'k = 0.1')
text(0.35,1.5,'k = 0.5')
text(0.55,2.5,'k = 1.0')
text(0.75,10,'k = 2.0')
text(0.7,7.5,'(unstable)')
