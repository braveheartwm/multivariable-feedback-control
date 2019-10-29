%System analysis for loop-shaping examples
%It is assumed that the plant is G=g, disturbance TFM is Gd=gd and 
%controller is K=c.
%It produce the system matrices, L=l, S=s, T=t,
%step response y=Tr, step disturbance response, yd=SGdr with time vector time.
%It also gives:Mt,Ms,Gm,Pm,Wc,W180,tr:y(tr)=0.9,Ov=max(y),
%Ymax=max(yd) and Y3=yd(3).
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: sysanaly.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

L=G*K;
S=1/(1+L);
T=1-S;
time=0:0.01:3;
yd=step(S*Gd,time);
Ymax=max(yd);Y3=yd(length(yd));
y=step(T,time);
k=length(y);
Ov=max([y;1]);[x,x,x,tr]=margin(y(2:k)/0.9,y(2:k),time(2:k));
MT=norm(T,inf,1e-4);
MS=norm(S,inf,1e-4);
[Gm,Pm,W180,Wc]=margin(L);

