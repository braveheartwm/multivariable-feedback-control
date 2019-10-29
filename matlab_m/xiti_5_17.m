clear all;
close all;

K1 = 0.24;
theta1 = 1;
K2 =38;
theta2 = 5;
T = 2;
G = tf(K2,conv([30 1],[T 1]),'iodelay',0.5);
H = K1*tf([1],[1],'ioDelay',theta1);
Gd = G*H;
G0 =tf(1,1);
bode(Gd,G,G0)
