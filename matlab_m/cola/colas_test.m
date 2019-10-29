
% This contains 
%    1. some information about the SIMULINK distillation files
%    2. some sample commands and comments

% 1. Some information.
%
%  States: X - 
%
% The numbering of the states are the same as for colamod.m. 
% The first 41 states are compositions of light component with 
% reboiler/bottom stage as X(1) and condenser as X(41). State X(42)
% is the holdup in the reboiler/bottom stage and X(82) is hold-up in the
% condenser. 
%
%  Inputs: U -
%
% U(1) - reflux L, U(2) - boilup V, U(3) - top or distillate product low D, 
% U(4) - bottom product flow B, U(5) - feed rate F, 
% U(6) - feed composition, zF and U(7) - liquid feed fraction, qF. 
%
%  Outputs:  y - 
%
% y(1) - top composition, y(2) - bottom composition, y(3) -
% condenser holdup, y(4) - reboiler holdup, 
% y(5:45) - tray compositions, e.g. y(5) is reboiler composition and y(45)
% is the condenser composition. 
%
% With the compositions available on all the stages the temperatures can be
% calculated by multiplying the value of the composition by 13.5. 
%

%
% 2. Useful commands using SIMULINK.
%

clear all

% To bring up the non-linear open loop model of distillation column
% type

colas_nonlin

% and the open-loop simulink model appears in the window.
% Choose <Start> from the <Simulation> pulldown menu,
% or press the key kombination <ctrl> <t>.
% The simulation with an 1\% increase in the feed
% flowrate F starts.
%
plot(t,y1)  % Plot the top composition. This plot also appear as 
            % the simulation propagates.
plot(t,y2)  % Plot the bottom composition. 
plot(t,y3)  % Plot the condenser level.
plot(t,y2)  % Plot the reboiler level.
%
% The SIMULINK model is the file colas.m which calls the 
% colamod.m to get the time derivaties of the states.
% An alternative to do the simulation from the SIMULINK 
% window is to simulate directly in the MATLAB command window;

clear all
load cola_init              % load the initial conditions.
% set up the input, u(1) = L, u(2) = V, u(3) = D, u(4) = B,
% u(5) = F, u(6) = zF and u(7) = qF.
% 
Us = [   0 Uinit'; 
      0.01 [Uinit(1:4).' Uinit(5)*1.01 Uinit(6:7)']; 
       100 [Uinit(1:4).' Uinit(5)*1.01 Uinit(6:7)']];
% At time t = 0 we have the initial condition.
% At time t = 0.01 we increase the feedrate u(5) with 1%
% and keep this value throughout the simulation time of 
% 100 [min]. 
[t,x,y] = rk45('colas', 100,[],[1e-5 0.01 10],Us);
% Note that the simulation may take some time to complete.
plot(t,y(:,[1])) % Plot the top composition.  
plot(t,y(:,[2])) % Plot the bottom composition.
%
% To linearize the non-linear model.
%
clear all
load cola_init % Load initial conditions.
% Contains; Xinit and Uinit
[A,B,C,D] = linmod('colas',Xinit,Uinit);
%
% Type help linmod to get more information, 
% i.e. to specify the pertubation levels etc.
%

%
% L,V configuration.
%
% Type 
colas_lv_nonlin
% and the partially controlled SIMULINK model appears 
% in the window. The condenser holdup is controlled using 
% a proportional controller with gain -10, and the reboiler 
% holdup is also controlled using a proportional controller 
% with gain -10.
% Choose <Start> from the <Simulation> pulldown menu,
% or press the key kombination <ctrl> <t>.
% The simulation with an 1\% increase in the feed
% flowrate F starts.
plot(t,y1)  % Plot the top composition. 
plot(t,y2)  % Plot the bottom composition.
plot(t,y3)  % Plot the condenser level.
plot(t,y4)  % Plot the reboiler level.

%
%
% Linearized versions. 
% Note that the variables are deviation variables.
%
clear all
load cola_init
[A,B,C,D] = linmod('colas',Xinit,Uinit); % Get the linearized model.
% Type 
colas_lin
% and the linearized model open loop model appears. 
% Simulate by pressing <ctrl> <t> and plot the result as earlier.
plot(t,y1)  % Plot the top composition. 
plot(t,y2)  % Plot the bottom composition.
plot(t,y3)  % Plot the condenser level.
plot(t,y2)  % Plot the reboiler level.


% Written by Kjetil Havre, Feb. 1997

