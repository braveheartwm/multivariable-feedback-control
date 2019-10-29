function xprime=cola_db(t,X) 
% sample usage:   [t,x]=ode15s('cola_dv',[0 5000],0.5*ones(1,82));
%
% cola_db - Subroutine for simulation with DB-configuration.
%
%
%            Inputs are distillate (D) and bottoms (B). Disturbances
%            are feedrate and feed composition. These are set by directly
%            altering 'cola_dv.m'. Outputs are liquid composition and
%            liquid hold up for stages 1 through NT, given in x. 

% Number of stages in the column
NT=41;

% Inputs and disturbances
D=0.5;                               % Distillate
B=0.5;                               % Bottoms
F=1.0 + 0.01;                        % Feedrate
zF=0.5;                              % Feed composition
qF=1.0;                              % Feed liquid fraction

% P-Controllers for control of reboiler and condenser hold up.
KcB=10;  KcD=10;         % controller gains
MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
Ls=2.70629; VBs=3.20629; % Nominal flows
MB=X(NT+1);  MD=X(2*NT); % Actual reboiler and condenser holdup
LT=Ls+(MD-MDs)*KcD;      % Reflux
VB=VBs+(MB-MBs)*KcB;     % Boilup     

% Store all inputs and disturbances
U(1)=LT; U(2)=VB; U(3)=D; U(4)=B; U(5)=F; U(6)=zF; U(7)=qF;

xprime=colamod(t,X,U);

