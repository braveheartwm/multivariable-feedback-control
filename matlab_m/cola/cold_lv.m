function xprime=cold_lv(t,X) 
% sample usage:   [t,x]=ode15s('cold_lv',[0 5000],0.5*ones(1,222));
%
% cold_lv - Subroutine for simulation with LV-configuration.
%           It calls the model colamod, and 
%           includes control of condenser and reboiler level 
%           using two P-controllers with the LV-configuration. 
%
%            Inputs are reflux (LT) and boilup (VB). Disturbances
%            are feedrate and feed composition. These are set by directly
%            altering 'cold_lv.m'. Outputs are liquid composition and
%            liquid hold up for stages 1 through NT, given in x. 

% Number of stages in the column
NT=111;

% Inputs and disturbances
LT=11.8615732;                          % Reflux
VB=12.476099;                          % Boilup
F=1.0 + 0.00;                        % Feedrate
zF=0.65;                              % Feed composition
qF=1.0;                              % Feed liquid fraction

% P-Controllers for control of reboiler and condenser hold up.
KcB=10;  KcD=10;         % controller gains
MDs=60; MBs=20;        % Nominal holdups   
Ds=0.614525139; Bs=1-Ds;          % Nominal flows
MB=X(NT+1);  MD=X(2*NT); % Actual reboiler and condenser holdup
D=Ds+(MD-MDs)*KcD;       % Distillate flow
B=Bs+(MB-MBs)*KcB;       % Bottoms flow     

% Store all inputs and disturbances
U(1)=LT; U(2)=VB; U(3)=D; U(4)=B; U(5)=F; U(6)=zF; U(7)=qF;

xprime=coldmod(t,X,U);

