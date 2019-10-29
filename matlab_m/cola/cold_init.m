
% Generate steady-state data for column D
% Saves data in file cold_init.mat

%% Do this by simulating 20000 min with stabilized LV-model:

clear all
xoinit=0.5*ones(1,111);
Moinit=1*ones(1,109);
xallinit = [xoinit 20 Moinit 60];
[t,x]=ode15s('cold_lv',[0 20000],xallinit); 
%Xinit = x(sel(size(x),1,1),:)';
lengthx=size(x); Xinit = x(lengthx(1),:)';

% Nominal inputs
LT=11.8615732;                          % Reflux
VB=12.476099;                          % Boilup
F=1.0;                               % Feedrate
D=0.614525139;                               % Distillate
B=F-D;                               % Bottoms
zF=0.65;                              % Feed composition
qF=1.0;                              % Feed liquid fraction
Uinit = [ LT VB D B F zF qF]';

clear x; clear t; clear LT; clear VB; clear D; clear B;
clear F; clear zF; clear qF;

save cold_init  % saves in cold_init.m


