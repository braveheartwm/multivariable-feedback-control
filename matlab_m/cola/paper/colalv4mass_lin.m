function [f,g]=colalv4mass_lin(X,u)
%
% cola_lin_lv - This function designed for use with 'cola_linearize.m' to
%        create a linear model of a LV-distillation column.
%
%     x - states (liquid composition and liquid hold up)
%     u - inputs and disturbances (reflux, boilup, feedrate and feed
%         composition)
%     g - outputs (distillate and bottoms compostition)
%
%  NOTE: Everything here is row vectors rather than columnn vectors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NT=41;

% Splitting the states
x=X(1:NT);                          % Liquid composition
M=X(NT+1:2*NT);                     % Liquid hold up

% Inputs and disturbances
Mwl=30;
yd = X(NT); Mw = Mwl*yd + 40*(1-yd);           % Mole weight 
LT = u(1)/ Mw ;                         % Molar reflux
VB=u(2);  % Boilup
F= u(3);  % Feedrate
zF=u(4);  % Feed composition
qF=1.0;   % Feed liquid fraction 
          % Use qF=u(5) if qF is to be included in linear model

% P-Controllers for control of reboiler and condenser hold up.
KcB=10;  KcD=10;         % controller gains
MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
Ds=0.5; Bs=0.5;          % Nominal flows
MB=X(NT+1);  MD=X(2*NT); % Actual reboiler and condenser holdup
D=Ds+(MD-MDs)*KcD;       % Distillate flow
B=Bs+(MB-MBs)*KcB;       % Bottoms flow     

% Store all inputs and disturbances
U(1)=LT; U(2)=VB; U(3)=D; U(4)=B; U(5)=F; U(6)=zF; U(7)=qF;

% This variable is not used
t=0;

xprime=colamod(t,X',U');

% Output
f=xprime';
g=[x(NT),x(1)];


