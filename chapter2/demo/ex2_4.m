clear all;
close all;
clc;

A = 0.1;
Wb = 20;
M = 10;
s = tf('s');
Wp1 = (s/M+Wb)/(s+Wb*A);
Wp2 = (s/sqrt(M)+Wb)^2/(s+Wb*sqrt(A))^2;
bode(1/Wp1,1/Wp2);grid;