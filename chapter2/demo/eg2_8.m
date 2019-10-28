clear all;
close all;
clc;

s = tf('s');
G = 3*(-2*s+1)/(10*s+1)/(5*s+1);

k = 0.05;
L = 3*k*(-2*s+1)/s/(2*s+1)/(0.33*s+1);

Kc = L/G;

T = L/(1+L);
S = 1 - T;
U = Kc*S;
step(T,U);
figure;
bode(Kc);grid on;