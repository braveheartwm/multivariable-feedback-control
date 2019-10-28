clear all;
close all;
clc;

s = tf('s');
G = 4/(s-1)/(0.02*s+1)^2;


Kc = 1.25*(1+1/1.5/s);
Gc = G*Kc/(1+G*Kc);


L = G*Kc;
S = 1/(1+L);
T = 1-S;
U = Kc*S;
figure;
step(Gc,4);hold on;
step(U,4)
figure;
bode(S,L,T);hold on;
grid on;
legend('S','L','T')
