clear all;
close all;
clc;

s = tf('s');
G = [(s+2.5)/(s-2) -(0.1*s+1)/(s-2);
    (s+2.5)/(0.1*s+1) 1];
[h1,h2]= hankelsv(G)

ksmin = 1/min(h2)

p = 2;
gp = evalfr(G,p+0.001);
[U,S,V] = svd(gp);
up = V(:,1)
g11s = (s+2.5)/(s+2);
g12s = -(0.1*s+1)/(s+2);
g21s =(s+2.5)/(0.1*s+1);
g22s = 1;
Gs = [g11s g12s;g21s g22s];
ksmin = norm(up'*inv(evalfr(Gs,p)))
% calculate pole output direction
% [num,den] = tfdata(G)
% 
% numMatrix = [num{1,1} num{1,2};num{2,1} num{2,2}]
% 
% 
% % [A,B,C,D] = tfmss(num,den)
% 
% 
% % clear all;
% close all;
% clc;
% d = [6 11 6];
% beta0 = [6,2;-6,-3];
% beta1 = [5,3;-5,-4];
% beta2 = [1,1;-1,-1];
% beta3 = [1,0;1,1];
% k = beta3;
% h = [beta0 beta1 beta2]';
% 
% [A,B,C,D] = tfmss(d,h,k)

