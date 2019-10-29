function xprime=cola4(t,X) 
% WITH 1% CHANGE IN FEED RATE
% sample usage:   [t,x]=ode15s('cola4',[0 5000],0.5*ones(1,82));
%
% cola4 - Subroutine for simulation using MATLAB integration routine
%           It calls the distillation model colamod 
%
%            SYNTAX: [t,x]=ode15s('cola4',tspan,Xinit,options);
%
%            tspan -   [t_start,t_stop]
%            Xinit -   column vector containing initial liquid composition
%                      for stages 1-NT and initial liquid hold up for stages
%                      1-NT.
%            options - see 'help ode15s'.
%
%            Inputs are reflux (LT), boilup (VB), distillate (D) and bottoms(B)
%            Disturbances are feedrate and feed composition.
%            These are set by directly altering 'cola4.m'.

% Inputs and disturbances
LT=2.70629;                          % Reflux
VB=3.20629;                          % Boilup
D=0.5;                               % Distillate
B=0.5;                               % Bottoms
F=1.0 + 0.01;                        % Feedrate
zF=0.5;                              % Feed composition
qF=1.0;                              % Feed liquid fraction

% Store all inputs and disturbances
U(1)=LT; U(2)=VB; U(3)=D; U(4)=B; U(5)=F; U(6)=zF; U(7)=qF;

xprime=colamod(t,X,U);

