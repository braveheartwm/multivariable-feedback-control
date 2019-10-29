clear all;
close all;
clc;

G = 3*tf([-2 1],conv([5 1],[10 1]));
Wi = tf([10,1/3],[10/5.25 1]);
delta = ultidyn('delta',[1 1],'Bound',1);
Gp = G*(1+Wi*delta);
K = tf([12.7 1],[12.7 0]);
L1 = Gp*1.13*K;
T1 = feedback(L1,1);
[Smarg1,Dstab1,Report1] = robuststab(T1)
L2 = Gp*1.13*K;
T2= feedback(Gp*0.31*K,1);
[Smarg2,Dstab2,Report2] = robuststab(T2)

T1n = feedback(1.13*G*K,1);
T2n = feedback(1.13*G*K*0.31,1);
bode(T1n,T2n,1/Wi)
