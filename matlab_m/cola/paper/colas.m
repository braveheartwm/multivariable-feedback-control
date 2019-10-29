function [sys,x0] = colas(t,x,u,flag)
%
%  Simulink interface to colamod.m     
%
% Inputs:    t    - time in [min].
%            X    - State, the first 41 states are compositions of light
%                   component A with reboiler/bottom stage as X(1) and 
%                   condenser as X(41). State X(42)is holdup in reboiler/
%                   bottom stage and X(82) is hold-up in condenser. 
%            U(1) - reflux L,
%            U(2) - boilup V,
%            U(3) - top or distillate product flow D,
%            U(4) - bottom product flow B,
%            U(5) - feed rate F,
%            U(6) - feed composition, zF.
%            U(7) - liquid feed fraction, qF.
%
% Outputs:   sys and x0 as described in the SIMULINK manual.
%            when flag is 0 sys contains sizes and x0 contains 
%            initial condition. 
%            when flag is 1, sys contains the derivatives,
%            and when flag is 3 sys contains outputs; 
%            y(1)    - top composition,
%            y(2)    - bottom composition,
%            y(3)    - condenser holdup,
%            y(4)    - reboiler holdup,
%            y(5:45) - tray composition, y(5) is reboiler y(45) is top.


NT = 41;

if abs(flag) == 1
  % Return state derivatives.
  sys = colamod(t,x,u);
elseif abs(flag) == 3
  % Return system outputs.
  sys(1,1) = x(NT);         % Top composition.
  sys(2,1) = x(1);          % Bottom composition.
  sys(3,1) = x(2*NT);       % Holdup in condenser.
  sys(4,1) = x(NT+1);       % Holdup in reboiler.
  sys(5:NT+4,1)= x(1:NT);   % Compositions.
elseif flag == 0
  % Ininitalize the system
  load cola_init
  x0 = Xinit;
  sys = [2*NT, 0, NT+4, 7, 0, 0];
else
  sys = [];  
end
  




