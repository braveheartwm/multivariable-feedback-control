% File cola_G4.m
% This requires the mu toolbox from MATLAB (sel,unpck, etc.)
% (This is a short version of the file cola_test.m)
% It generates and stores the 82 state linear column model in G4.mat
% which can be retrieved by load G4
%
% After you load G4 your variables are:
% 
% A         C         G4        So        
% B         D         Si        Xinit 
% 
% where:
%
% G4 - scaled linear model with 6 scaled inputs, 4 scaled outputs, 82 states
%      (in Mutools format)
% A,B,C,D - state space matrices for G4
% Si and So - input and outputs scalings. 
%    Note: Unscaled model is G4u = mmult(minv(So), G4, minv(Si));
% Xinit - nominal (initial) values for the states (G4 gives the
%         deviation from this nominal state)

clear  all
%---------------------------------------------
% First find steady-state 
%---------------------------------------------

% Do this by simulating 5000 min with stabilized LV-model:
[t,x]=ode15s('cola_lv',[0 5000],0.5*ones(1,82)'); 
Xinit = x(sel(size(x),1,1),:)'; 
  
%--------------------------------------
% Now linearize the model to obtain G4u
%--------------------------------------

% The open-loop model (with no levels closed; 6 inputs, 4 outputs)
Ls=2.70629; Vs=3.20629; Ds=0.5; Bs=0.5; Fs=1.0; zFs=0.5;
Uinit=[Ls Vs Ds Bs Fs zFs]';
[Au,Bu,Cu,Du]=cola_linearize('cola4_lin',Xinit',Uinit');
G4u =  pck(Au,Bu,Cu,Du);

%---------------------------------------
% Obtaining the scaled model G4
% --------------------------------------


% The following max. changes are used (for scaling the model):
Du = diag([1 1 1 1]);       % max inputs (scalings)
Dd = diag([ .2 .1]);        % max disturbances (scalings)
De = diag([0.01 0.01 1 1]); % max output errors (scalings)
% This implies the folling in terms of the scaled model G4:
   % Units for inputs (L,V,D,B):  1 = 1 kmol/min = F (the feed rate)
   % Units for disturbance 1 (F): 1 = 0.2 kmol/min (20% change)
   % Units for disturbance 2 (z_f): 1 = 0.1 mole fraction units (20% change)
  % Units for outputs 1 and 2 (y_d and x_b): 1 = 0.01 mole fraction units
   % Units for outputs 3 and 4 (M_d and M_b): 1 = 1 kmol
% The scaled model is then G4:
Si = daug(Du,Dd); So = minv(De);       % introduce scaling matrices
G4 = mmult(So, G4u, Si);

[A,B,C,D]=unpck(G4);

clear t;clear x;clear Si1;clear Si2;clear So1;clear So2
clear Au; clear Bu; clear Cu; clear Du; clear G4u; clear ans
clear Ls; clear Vs; clear Ds; clear Bs; clear Fs; clear zFs;

save G4
who

