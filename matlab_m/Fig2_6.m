%Figure 2.6 Step responses
%explained in Example 2.1 and Example 2.2

% Plant (2.26), G(s)=3(-2s+1)/(5s+1)(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Fig2_5.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

clear
s=tf('s');
g=3*(-2*s+1)/(5*s+1)/(10*s+1);

Kc=0.5;
gc=g*Kc;
S = 1/(1+gc);
T1 = 1-S;

Kc=1.5;
gc=g*Kc;
S = 1/(1+gc);
T2 = 1-S;

Kc=2.5;
gc=g*Kc;
S = 1/(1+gc);
T3 = 1-S;

Kc=3;
gc=g*Kc;
S = 1/(1+gc);
T4 = 1-S;

clf
t=0:0.5:50; u=ones(size(t));   %input
lsim(T1,'b',T2,'g',T3,'r',T4,'c:',u,t)
hold
plot([32.5 47],[-0.26 -0.26],'g',[32.5 33.5],[-0.26 -0.23],...
    'g',[32.5 33.5],[-0.26 -0.29],'g',[46 47],[-0.23 -0.26],...
    'g',[46 47],[-0.29 -0.26],'g')
hold
axis([0 50 -0.5 2.5]);
text(11,1.95,'Kc=2.5'); text(14,1.3,'Kc=1.5'); text(8,0.2,'Kc=0.5');
text(10,2.35,'Kc=3 (Unstable)');
text(40,-0.15,'Pu');
xlabel('Time');ylabel('y');
title('')
