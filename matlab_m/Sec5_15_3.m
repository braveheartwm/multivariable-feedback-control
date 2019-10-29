% Section 5.15.3 Application: Neutralization process
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: NeutProc.m,v 1.3 2004/02/03 14:12:00 vidaral Exp $

clear all; close all;

%Plant (5.114): G(s)=-2kd/(1+tauh*s), Gd(s)=kd/(1+tauh*s), kd=2.5e6, tauh=1000
kd=2.5e6;tauh=1000;
G=-2*kd*tf(1,[tauh 1]);
Gd=kd*tf(1,[tauh 1]);
w=logspace(-4,5,61);

%Figure 5.22
figure(1);clf;
bodemag(G,'-',Gd,'-',tf(1),'--',w); hold on; plot([2500 2500],[0.1 1],':k'); hold off;
axis([1e-4 1e5 1e-2 1e7]);
xlabel('Frequency [rad/s]');ylabel('Magnitude');
text(0.1, 3000,'|Gd|');text(.1, 2e5,'|G|');
text(1500,0.02,'\omega_d');

% n Buffer tanks, %Figure 5.24
figure(2)
hold on;
for n=1:4
	hn=(tf([0 1],[1/n 1]))^n;
	bodemag(hn,'-',w);
end
hold off;
axis([1.e-1 1e2 1.e-5 1]);
xlabel('Frequency x \tau_h');ylabel('Magnitude');
text(14,.1,'n=1');text(16,0.02,'n=2');
text(18,0.005,'n=3');text(11,0.002,'n=4');
text(1,1e-2,'|h_n|');
