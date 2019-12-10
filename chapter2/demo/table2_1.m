clear all;
close all;
clc;
tau=1;
zeta=0.1;
t=0:0.01:100;
T = nd2sys(1,[tau*tau 2*tau*zeta 1]); 
S = msub(1,T);
[A,B,C,D]=unpck(T); 
y1 = step(A,B,C,D,1,t);
overshoot=max(y1),tv=sum(abs(diff(y1)))

Mt=hinfnorm(T,1.e-4),Ms=hinfnorm(S,1.e-4)