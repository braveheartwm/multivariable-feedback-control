clear all;
close all;
clc;

L1 = tf(0.5,[1 0]);
L2 = L1*tf([-1 1],[1 1]);
Wu = 0.5*tf([1,0],[1 1]);
Wp = 0.25+tf(0.1,[1 1e-6]);

delta = ultidyn('delta',[1 1]);
S1 = inv(1+L1+delta*Wu);
[Smarg1,Dstab1,Report1] = robustperf(Wp*S1)
S2 = inv(1+L2+delta*Wu);
[Smarg2,Dstab2,Report2] = robustperf(Wp*S2)

