% Example 8.9 Produces fig. 8.11
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl8_9.m,v 1.1 2004/01/21 15:21:24 vidaral Exp $
clear all; close all;
omega = logspace(-3,2,1001);

tau = 75; 
G=tf([0 1],[tau 1])*[-87.8 1.4;-108.2 -1.4]; % Plant
K=tf([tau 1],[1 1e-6])*[-0.0015 0; 0 -0.075]; % Controller
wI=0.2*tf([5 1],[0.5 1]); % Input uncertainty
wIf=frd(wI,omega);

LI=(K*G);
SI=(inv(eye(2)+LI));
TI=(SI*LI);
Tjw=sigma(TI,omega);
blk=[1 1;1 1];
[bnds,muinfo]=mussv(frd(TI,omega),blk);
bodemag(1/wIf,'r--',bnds(1,1),'-',omega)
hold on
plot(omega,Tjw(1,:))
text( 1, 12, 'sv' ); text(0.002, 8, '1/WI' ); text( 1, 0.3, 'mu' );
hold off
axis auto