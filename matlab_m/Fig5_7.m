% Figure 5.7: "Ideal" sensitivity function for plants with RHP-zeros
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig5_7.m,v 1.6 2004/02/03 14:12:00 vidaral Exp $

clear all; close all;
% % set(cstprefs.tbxprefs,'MagnitudeUnits','dB','MagnitudeScale','linear',); % log-log
z=1;
T=1*tf([-1 1],[1 z]);
S=1-T;
w=logspace(-2,2,41);

% Defining the asymptotes in Figure (a) based on the "procedure" on page
% 20: s=2*s/s+1
w1=logspace(-2,0,41); w2=logspace(0,2,41);
s1=tf([2 0],1); [s1mag,s1phas]=bode(s1,w1); % for w < 1, the break frequency of the denominator of s
s2=tf(2,1); [s2mag,s2phas]=bode(s2,w2); % for w > 1

figure(1)
subplot(2,1,1);
bodemag(S,w); hold on;
plot(w,1,'--',w1,s1mag(:),'--',w2,s2mag(:),'--',[z/2 z/2],[1e-2 1],':k'); hold off;
axis([.01,100,.01,10]);
text(0.42,0.02,'z/2');
xlabel('Frequency x 1/z');ylabel('Magnitude |S|');
title('(a)');

% Plot of the ideal S for plant with complex RHP-zeros, x +- jy.
% Let R = y/x
R=0.1;  S1=(1-tf([1 -2 1+R^2],[1  2 1+R^2]));
R=1;    S2=(1-tf([1 -2 1+R^2],[1  2 1+R^2]));
R=10;   S3=(1-tf([1 -2 1+R^2],[1  2 1+R^2]));
R=50;   S4=(1-tf([1 -2 1+R^2],[1  2 1+R^2]));

subplot(2,1,2)
bodemag(S1,'-',S2,'-',S3,'-',S4,'-',tf(1),'--',w);
axis([.01,100,.01,5]);
text(30,3,'y/x =  50')
text(5,3,'y/x =  10')
text(0.3,0.5,'y/x =  1')
text(0.15,2,'y/x =  0.1')
xlabel('Frequency x 1/x');ylabel('Magnitude |S|');
title('(b)');
