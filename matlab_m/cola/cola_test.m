% This is the file cola_test.m 
% Contains a collection of useful commands for analyzing the distillation models
% The file requires the mu toolbox from MATLAB (sel,unpck, etc.), but
% can be easily changed to avoid this

% Simulate with full nonlinear model and generate linear model

%---------------------------------------------
% First find steady-state 
%---------------------------------------------

% Do this by simulating 20000 min with stabilized LV-model:
[t,x]=ode15s('cola_lv',[0 20000],0.5*ones(1,82)'); 
Xinit = x(sel(size(x),1,1),:)' % This data has been saved in cola_init.mat
                               % and can be retrieved using the command
                               % load cola_init
  
% Perform nonlinear simulation
[t,x]=ode15s('cola4',[0 500],Xinit);   
plot(t,x);  %Nothing happens so we are at steady-state...

% -------------------------------------------------------------
% Simulate a change in feed rate with level loops open 
% -------------------------------------------------------------
% FIRST GO INTO THE FILE cola.m and change F=1.00 to F=1.01. 
% SAVE THE FILE as cola_F1.m_and simulate:

[t,x]=ode15s('cola4_F1',[0 500],Xinit); 
t0 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting

% Plot reboiler holdup (Apparent delay of about 1.26 min)
plot(t0,M1); axis([0 3 0.5 0.52]); title('MB: reboiler holdup')

plot(t0,xB); title('xB: reboiler composition')
plot(t0,yD); title('yD: distillate composition')

%---------------------------------------------
% Nonlinear Simulation with other configurations (close level loops)
%---------------------------------------------

% The following is described under "Simulations with other configurations"

% FIRST GO INTO THE FILE cola_lv.m and change F to 1.01.  Then use 
[t,x]=ode15s('cola_lv_F1',[0 500],Xinit); 
tlv=t; M1lv=x(:,42); xBlv = x(:,1); yDlv = x(:,41); %save data for plot 
% Do the same for other configurations (REMEMBER TO CHANGE F!!!)
% LB-configuration
[t,x]=ode15s('cola_lb_F1',[0 500],Xinit); 
tlb=t; M1lb=x(:,42); xBlb = x(:,1); yDlb = x(:,41);
% Double ratio
[t,x]=ode15s('cola_rr_F1',[0 500],Xinit);
trr=t; M1rr=x(:,42); xBrr = x(:,1); yDrr = x(:,41); 

% Plot them together
plot(t0,M1,'-',tlv,M1lv,':',tlb,M1lb,'-.',trr,M1rr,'--'); 
title('MB: reboiler (bottom) holdup [mol]');

plot(t0,xB,'-',tlv,xBlv,':',tlb,xBlb,'-.',trr,xBrr,'--'); 
title('xB: reboiler (bottom) composition');

plot(t0,yD,'-',tlv,yDlv,':',tlb,yDlb,'-.',trr,yDrr,'--'); 
title('yD: distillate (top) composition');

%--------------------------------------
% Now linearize the model to obtain G4u
%--------------------------------------

% The open-loop model (with no levels closed; 6 inputs, 4 outputs)
Ls=2.70629; Vs=3.20629; Ds=0.5; Bs=0.5; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola4_lin',Xinit',[Ls Vs Ds Bs Fs zFs]);
G4u =  pck(A,B,C,D);

eig(A)  % Compute eigenvalues

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

% To analyze this model further, and find the model for other configurations
% see the file Table 12.3.m

%-----------------------------------------------
% Reduced-order model
%----------------------------------------------

[syss,sysu]=sdecomp(G4); [syssb,hsig]=sysbal(syss); 
sys1=hankmr(syssb,hsig,6,'d');
G4_8=madd(sys1,sysu);                   % scaled reduced-order model
G4u_8=mmult( minv(So), G4_8, minv(Si)); % scale back if desired

% Compare reduced-order and full model (not much difference!)
[A,B,C,D]=unpck(G4u);
% Simulate step in F of magnitide 1
[Y,X,T]=step(A,B,C,D,5); plot(T,Y(:,1),T, Y(:,2))  % Simulation full model
[A8,B8,C8,D8]=unpck(G4u_8); 
[Y8,X8]=step(A8,B8,C8,D8,5,T);                     % Simulate 8-state model
plot(T,Y(:,1),T,Y8(:,1));                          % not much difference!

eig(A8)                       % Compute eigenvalues for 8-state model
[T,P] = eig(A8); YP = C8*T    %   output pole directions (the two integrators 
                              %    are associated with y3=Md and y4=Mb)
eig(A8')
[Q,P] = eig(A8'); UP = B8'*Q  %   input pole directions


%-----------------------------------------------
% Linearized model for LV- and double ratio configuration
%------------------------------------------------

%The linear model for the LV-configurtation (4 inputs, 2 outputs)
Ls=2.70629; Vs=3.20629; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola_lv_lin',Xinit',[Ls Vs Fs zFs]);
Glvu =  pck(A,B,C,D);
% simulate step in F of magnitude 1
step(A,B,C,D,3);
%  --- se file cola_commands.m for further analysis of the linear LV-model

%The linear model for the double ratio-configurtation (4 inputs, 2 outputs)
R1s = 2.70629/0.5; R2s=3.20629/0.5; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola_rr_lin',Xinit',[R1s R2s Fs zFs]);
Grru =  pck(A,B,C,D);
% simulate step in F of magnitude 1
step(A,B,C,D,3);

