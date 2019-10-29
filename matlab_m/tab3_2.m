clear all;
close all;
clc;
G= tf(1,[1,1])
systemnames='G';
inputvar='[d(1);r(1);n(1);u(1)]';
input_to_G='[u]';
outputvar='[G+d-r;r-G-d-n]';
sysoutname='P';
sysic;