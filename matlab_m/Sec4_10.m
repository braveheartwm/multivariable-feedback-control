%Examples to show the calculation of system norms in Section 4.10
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Sec4_10.m,v 1.3 2004/04/02 17:35:24 vidaral Exp $
close all; clear all;
G=tf(1,[1 0.5]);
G2=G*G;
% h2:   get 1 and 1.414
norm(G,2)
norm(G2,2)

%hinf: get 2 and 4
norm(G,inf)
norm(G2,inf)

% hankel: get 1 and 2.414
hsig=hankelsv(G); hsig.stab
hsig2=hankelsv(G2); max(hsig2.stab)

epps=0.0001;
f1 =tf(1,[epps 1]);    
f2 = tf([epps 0],[1 epps 1]);
% get 70.71 (h2) and 1 (hinf) and 0.5 (hankel)
norm(f1,2)
norm(f1,inf)
hsig=hankelsv(f1); hsig.stab


% get 0.007 (h2) and 1 (hinf) and 0.5 (hankel)
norm(f2,2)
norm(f2,inf)
hsig2=hankelsv(f2); max(hsig2.stab)






