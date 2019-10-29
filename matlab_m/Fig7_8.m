% Figure 7.8 
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig7_8.m,v 1.2 2004/04/14 09:37:41 vidaral Exp $
close all; clear all;
% Time delay
w=logspace(-2,2,801);
theta=1;
l_delay=1-frd(tf(1,1,'ioDelay',theta),w);
subplot(121)
bodemag(abs(l_delay),'k-');hold on
plot([pi 100],[2 2],'k--',[0.01 100],[1 1],'k:',[1 1],[0.02 1],'k:')
text(0.7,0.015,'1/\theta_{max}')
ylabel('Magnitude'), xlabel('Frequency')

% Neglected lag
omega=logspace(-2,2,81);
tau=1;
l_lag=1-tf(tf(1,[tau 1]));
subplot(122)
bodemag(l_lag,'k-');hold on
plot([0.01 100],[1 1],'k:',[1 1],[0.02 1],'k:')
text(0.7,0.015,'1/\tau_{max}')
axis([0.01 100 0.01 10])
ylabel('Magnitude'), xlabel('Frequency')
