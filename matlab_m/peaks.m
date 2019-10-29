% peaks.m 
% Calculated peak values and total variation of a prototype second-order system
% T(s)=1/(tau^2*s^2+2*tau*zeta*s+1)
% Variables: tau and zeta must be defined before using the procedure.
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: peaks.m,v 1.3 2004/04/14 08:11:00 vidaral Exp $

T = tf(1,[tau*tau 2*tau*zeta 1]); 
S = 1-T; 
[A,B,C,D]=ssdata(T);
if zeta<1,
%To get accurate results for small zeta (<0.1), you have to reduce 
%the value of 750, to extend simulation time, or decrease the value
%of 0.01, to get shorter integral step. Warning: This will need large
%memory and take a long time for simulation.
    t=0:0.01:750;
else,
    t=0:100/30:100;
end
y=step(A,B,C,D,1,t);
Ov=max([y;1]); %Overshoot
Tv=sum(abs(diff([y;1]))); %Total Variation
Mt=norm(T,inf,1e-4);
Ms=norm(S,inf,1e-4);


