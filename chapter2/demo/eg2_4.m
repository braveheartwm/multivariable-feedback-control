clear all;
close all;
clc;

s = tf('s');
G = 3*(-2*s+1)/(10*s + 1)/(5*s + 1);
kc = [0.5 1.5 2.5 3];
for i = 1:length(kc)
  L = G*kc(i);
  T = L/(1+L);
  Gc = T;
  step(Gc,50);hold on;
  grid on;
  ylim([-0.5,2.5]);
end

Kc = 1.136*(1+1/12.7/s);
Gc = G*Kc/(1+G*Kc);
figure;
step(Gc,80);

L = G*Kc;
S = 1/(1+L);
T = 1-S;
figure;
bode(S,L,T);hold on;
grid on;
legend('S','L','T')
