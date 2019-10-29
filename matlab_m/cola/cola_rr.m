function xprime=cola_rr(t,X) 
% sample usage:   [t,x]=ode15s('cola_Rv',[0 5000],0.5*ones(1,82));
%
% cola_rr - Subroutine for simulation with L/D-V/B-configuration.
%
%           There are many ways to implement the double ratio configuration.
%           In this case we set the ratios L/D and V/B externally, and let
%           the levels be controlled using D and B (controlling the levels
%           like this may not be the best solution; but for tight level
%           control it makes no difference)
%
%            Inputs are reflux (LT) and boilup (VB). Disturbances
%            are feedrate and feed composition. These are set by directly
%            altering 'cola_dv.m'. Outputs are liquid composition and
%            liquid hold up for stages 1 through NT, given in x. 

% Number of stages in the column
NT=41;

% Inputs and disturbances
RT=2.70629/0.5;                      % L/D
RB=3.20629/0.5;                      % V/B
F=1.0 + 0.00;                        % Feedrate
zF=0.5;                              % Feed composition
qF=1.0;                              % Feed liquid fraction

% P-Controllers for control of reboiler and condenser hold up.
KcB=10;  KcD=10;         % controller gains
MDs=0.5; MBs=0.5;        % Nominal holdups - these are rather small  
Ds=0.5; Bs=0.5;          % Nominal flows
MB=X(NT+1);  MD=X(2*NT); % Actual reboiler and condenser holdup
D=Ds+(MD-MDs)*KcD;       % Reflux
B=Bs+(MB-MBs)*KcB;       % Bottoms flow     

% Resulting reflux and boilup
LT = RT*D;
VB = RB*B;

% Store all inputs and disturbances
U(1)=LT; U(2)=VB; U(3)=D; U(4)=B; U(5)=F; U(6)=zF; U(7)=qF;

xprime=colamod(t,X,U);

