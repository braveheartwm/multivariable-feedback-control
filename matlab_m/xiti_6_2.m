clear all;
close all;
clc;
 
s = tf('s');
theta = 1;
k = 1/theta;

t = 0:0.01:30';
n = length(t);

u = ones(n,1);
r1 = [u u];

G = [1 1;exp(-theta*s) 1];
K = k/s*eye(2);
T = G*K*inv(eye(2)+G*K);
y = lsim(T,r1,t);
plot(t,y(:,1),t,y(:,2))