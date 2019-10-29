
% Generate steady-state data for column A
% Saves data in file cola_init.mat

%% Do this by simulating 20000 min with stabilized LV-model:

clear all
[t,x]=ode15s('cola_lv',[0 20000],0.5*ones(1,82)'); 
%Xinit = x(sel(size(x),1,1),:)';
lengthx=size(x); Xinit = x(lengthx(1),:)';

% Nominal inputs
LT=2.70629;                          % Reflux
VB=3.20629;                          % Boilup
D=0.5;                               % Distillate
B=0.5;                               % Bottoms
F=1.0;                               % Feedrate
zF=0.5;                              % Feed composition
qF=1.0;                              % Feed liquid fraction
Uinit = [ LT VB D B F zF qF]';

clear x; clear t; clear LT; clear VB; clear D; clear B;
clear F; clear zF; clear qF;

save cola_init  % saves in cola_init.m


