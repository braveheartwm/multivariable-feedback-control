clear all;
close all;
clc;

load aeroK;
% two degree H_inf controller 
K1 = K(1:3,1:3);
K2 = K(1:3,4:6);

% control plant
G = Gs;

Tref = Mo;%reference model

systemnames='G K1 K2';
inputvar = '[r(3)]';
outputvar='[G]';
input_to_G = '[ K1 + K2 ]';
input_to_K1 = '[r]';
input_to_K2 = '[G]';
sysoutname = 'M';
cleanupsysic = 'yes';
sysic;

