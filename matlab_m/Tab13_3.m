% This is Table 12.3 in Skogestad and Postlethwaite (1996).
% 
% This takes the detailed linear model with 82 states and 4 inputs
% given in the MATLAB file G4.mat and generates models for the 
% LV- DV- and DB-configurations.
%
% Uses MATLAB Mu toolbox
%

%Load linear model of distillation column:
load G4
% G4 may be generated from the nonlinear model as shown in the file cola_G4.m
%
% G4 is the detailed scaled linear model (in deviation variables) with:
%     4 inputs, 2 disturbances, 4 outputs and 82 states.
% The two disturbances appear as inputs 5 and 6 in the model.
% Order states: composition and holdup on each stage (starting from top) 
%   The states can in this model (which is in deviation variables) 
%   be initialized to zero. 
% The following scalings have been used to obtain G4 (see p. 6 in the book)
%  Du = diag([1 1 1 1]);       % max inputs (scalings)
%  Dd = diag([ .2 .1]);        % max disturbances (scalings)
%  De = diag([0.01 0.01 1 1]);  % max output errors (scalings)
% This implies the folling in terms of the scaled model G4:
   % Units for inputs (L,V,D,B):  1 = 1 kmol/min = F (the feed rate)
   % Units for disturbance 1 (F): 1 = 0.2 kmol/min (20% change)
   % Units for disturbance 2 (z_f): 1 = 0.1 mole fraction units (20% change)
   % Units for outputs 1 and 2 (y_d and x_b): 1 = 0.01 mole fraction units
   % Units for outputs 3 and 4 (M_d and M_b): 1 = 1 kmol
%   Units for states: x1, x2, ..., x41: 1 = 1 mole fraction unit, 
%                     x42, x43, ..., x82: 1 = 1 kmol
%
% REMARK 1: The model G4 is based on the "unscaled" model G4u:
%      Si = daug(Du,Dd); So = minv(De);   % scaling matrices
%      G4 = mmult (So, G4u, Si);
%   That is, to "unscale" the model back to its original units use:
%      G4u = mmult( minv(So), G4, minv(Si)); 
%
% REMARK 2: The following nominal values are just for your information 
%           (do not affect the linear model):
%   The nominal values for inputs: 
%       L=2.71 mol/s,  V=3.21 mol/s, D=0.5 mol/s, B=0.5 mol/s
%   The nominal values for disturbances: 
%       F = 1 mol/s,  z_F = 0.5 mole fraction units
%   The nominal values for outputs: 
%       y_d=0.99 and x_b=0.01 mole fraction units,
%       M_d =0.5 mol and M-b=0.5 mol
%   The nominal values for the states: 
%       For the 41 compositions (x1, x2,... x41): see file cola.dat (see X...)  
%       For the 41 holdups (x42, x43, ....  x82): 0.5 mol   

% -----------------------------------------------------------


% Now generate the linear model for the LV Configuration...
% First put in the level controller for D and B
Kd = 10; Kb = 10; % P-controllers for levels (bandwidth = 10 rad/min)
systemnames = 'G4 Kd Kb';
inputvar = '[L(1); V(1); d(2)]'; 
outputvar = '[G4(1);G4(2)]';
input_to_G4 = '[L; V; Kd; Kb; d ]';
input_to_Kd = '[G4(3)]'; 
input_to_Kb = '[G4(4)]';
sysoutname ='Glv'; cleanupsysic='yes'; sysic;

%minfo(Glv); % Glv has 4 inputs (L,V + disturbances) and 2 outputs
%[Glvb,hsig]=sysbal(Glv); %minfo(Glvb)
% G = sel(Glvb,':',[1 2]);
% Gd= sel(Glvb,':',[3 4]);
[Glvb,hsig]=balancmr(Glv,42)
G=ss(Glvb.a,Glvb.b(:,[1 2]),Glvb.c,Glvb.d(:,[1 2]));
Gd=ss(Glvb.a,Glvb.b(:,[3 4]),Glvb.c,Glvb.d(:,[3 4]));
G0 = freqresp(G,0)  % steady-state gains  - a bit different from the book
%vrga(G0)
Gd0= freqresp(Gd,0)  % steady-state gains 

% This gives a model for the LV-configurations with 82 states, which
% for practical purposes may be be replaced by the 5 state model 
% used in the book.
% REMARK: 
% I am not quite sure how the 5-state model in the book was generated
% (and Morten Hovd cannot remember) but if you try
Glv5  = hankelmr(Glvb,5); 
% then you get something reasonably close to the 5-state model in the book
pole(Glv5)


%------------------------------------------------------------
% Now generate the linear model for the DV Configuration
% First put in the level controller for L and B
Kl = 10; Kb = 10; % P-controllers for levels (bandwidth = 10 rad/min)
systemnames = 'G4 Kl Kb';
inputvar = '[D(1); V(1); d(2)]'; outputvar = '[G4(1);G4(2)]';
input_to_G4 = '[Kl; V; D; Kb; d ]';
input_to_Kl = '[G4(3)]'; input_to_Kb = '[G4(4)]';
sysoutname ='Gdv'; cleanupsysic='yes'; sysic;

freqresp(Gdv,0)

%------------------------------------------------------------
% This is to generate the linear model for the DB-configuration
Kl = 10; Kv = 10;
systemnames = 'G4 Kl Kv';
inputvar = '[D(1); B(1); d(2)]'; outputvar = '[G4(1);G4(2)]';
input_to_G4 = '[Kl; Kv; D; B; d ]';
input_to_Kl = '[G4(3)]'; input_to_Kv = '[G4(4)]';
sysoutname ='Gdb'; cleanupsysic='yes'; sysic;


% -------------------------------------------------------------
% REMARK: We can easily generate the model for other configuration!
%

