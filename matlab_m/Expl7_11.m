% Example 7.11 Produces fig. 7.19 Robust performance problem
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl7_11.m,v 1.2 2004/04/14 09:37:41 vidaral Exp $
w=logspace(-2,2,51);
ru=0.5;
L1=tf(0.5,[1 0]);
L2=L1*tf([-1 1],[1 1]);
Wp=frd(0.25+tf(0.1,[1 0]),w);
Wu=frd(ru*tf([1 0],[1 1]),w);
S1=inv(1+L1);
S2=inv(1+L2);
RP=inv(abs(Wp)+abs(Wu));

% Figure 7.19
bodemag(S1,'k-',S2,'k-',RP,'k:',{10^(-2),100})
text(0.4,0.4,'|S1|');
text(0.7,3.5,'|S2|');
text(10,2,'1/(|Wp|+|Wu|)');

% Alternatively, using rcast
Du=ultidyn('DeltaU',[1 1]);
L1unc=L1+Du*Wu;
L2unc=L2+Du*Wu;
WuS1unc=Wp*inv(1+L1unc);
[PERFMARG1,WCU1,REPORT1,INFO1] = robustperf(WuS1unc);
disp('For L1:')
char(REPORT1)
WuS2unc=Wp*inv(1+L2unc);
[PERFMARG2,WCU2,REPORT2,INFO2] = robustperf(WuS2unc);
disp('For L2:')
char(REPORT2)

