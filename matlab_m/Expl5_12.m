% Example 5.12: Inputs for perfect control
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl5_12.m,v 1.2 2004/01/19 16:09:06 araujo Exp $


clear all

w=logspace(-2,2,101);
G=tf(40,conv([5 1],[2.5 1]));
Gd=3*tf([50 1],conv([10 1],[1 1]));

[magG,phaseG]=bode(G,w);
[magGd,phaseGd]=bode(Gd,w);
z(1:length(w))=1;
l1=3e-1:1e-1:13;
one(1:length(l1))=0.38;
l2=3e-1:1e-1:1;
one1(1:length(l2))=15;

%Figure 5.15
loglog(w,squeeze(magG),w,squeeze(magGd),w,z,':k',one,l1,':k',one1,l2,':k');
axis([0.01 100 0.1 100]);
xlabel('Frequency [rad/s]');
ylabel('Magnitude');
text(0.38,3e-1,'\omega_1');
text(15,3e-1,'\omega_d');
