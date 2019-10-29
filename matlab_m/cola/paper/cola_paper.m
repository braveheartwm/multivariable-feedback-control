
% This file is used to generate the examples and figures in:
% S. Skogestad: Dynamics and control of distillation columns -
% A tutorial introduction. Presented at Distillationa and Absorbtion
% 97, Maastricht, Netherlands, 8-10 Sept. 1997

%---------------------------------------------
% First find steady-state 
%---------------------------------------------

% Do this by simulating 20000 min with stabilized LV-model:
[t,x]=ode15s('cola_lv',[0 20000],0.5*ones(1,82)'); 
Xinit = x(sel(size(x),1,1),:)'; % This data has been saved in cola_init.mat
                               % and can be retrieved using the command
                               % load cola_init
  
% Perform nonlinear simulation
[t,x]=ode15s('cola4',[0 500],Xinit);   
plot(t,x);  %Nothing happens so we are at steady-state...
%save Xinit Xinit
    
% ---------------------------------------------
% Figure 2: Composition profiles
% -------------------------------------------

xprofile0 = Xinit(1:41);
stages = [1:1:41]';

% Now after change in external flows (increase L by 0.02 with V constant)
[t,x]=ode15s('colalv1large',[0 20000],Xinit); 
xdum = x(sel(size(x),1,1),:)';
xprofile1 = xdum(1:41);

% Now after change in internal flows (increase both L and V by 1)
[t,x]=ode15s('colalv2large',[0 20000],Xinit); 
xdum = x(sel(size(x),1,1),:)';
xprofile2 = xdum(1:41);

plot(xprofile0,stages); hold
plot(xprofile1,stages,'-.'); 
plot(xprofile2,stages,'--'); hold off;
axis([0 1 1 41])
ylabel('Stage no.')
xlabel('Composition')
hl = legend('Nom','Ext','Int')
set(hl,'Visible','off','Position', [0.15 0.7 0.2 0.2] )
%print -deps profile
%!mv profile.eps ~skoge/latex/work/fig/profile.eps

% Logarithmic compositions
Xprofile0 = log(xprofile0(:)./(1-xprofile0(:)) );
Xprofile1 = log(xprofile1(:)./(1-xprofile1(:)) );
Xprofile2 = log(xprofile2(:)./(1-xprofile2(:)) );

plot(Xprofile0,stages); hold;
plot(Xprofile1,stages,'-.'); 
plot(Xprofile2,stages,'--'); hold off
axis([-6 6 1 41 ])
ylabel('Stage no.')
xlabel('Logaritmic composition')
hl = legend('Nom','Ext','Int')
set(hl,'Visible','off','Position', [0.15 0.7 0.2 0.2] )
%print -deps profilelog
%!mv profilelog.eps ~skoge/latex/work/fig/profilelog.eps
hold off

% ------------------------------------------------
% Figure 4: Respons to bottom liquid flow L_B = L_2
% ------------------------------------------------
clear all
load Xinit
tsim=7;
[t2,x]=ode23s('colalv2',[0:0.01:tsim],Xinit,[1e-5]); 
[t2n,xn]=ode15s('colalv2n',[0:0.01:tsim],Xinit); 
    taul=0.063;     	% time constant for liquid dynamics (min)

LB = x(:,43)/taul; % This is L - L0i + M0i/taul; assumes lambda=0
LT=2.70629*1.1;                          % Reflux

plot(t2,LB)

t2 = [-1; t2];
L = [LB(1);LB];
LTv = [ 0;0;LT-LT/1.1;LT-LT/1.1];
tv  = [-1;0; 0;tsim];

plot(t2,L-L(1),tv,LTv,'--');
xlabel('Time'); ylabel('Change in L');
hl = legend('LB','LT')
set(hl,'Visible','off','Position', [0.7 0.475 0.11 0.2] )
set(gca,'FontSize',14)
axis([-1 7 -0.05 0.3])
%print -deps figlv2a
%!mv figlv2a.eps ~skoge/latex/work/fig/figlv2a.eps
%!mv figlv2a.eps ../latex/fig

%-----------------------------------------------------------------
% (NO FIGURE): K2 effect (parameter lambda). Incerease V by 0.1
% ---------------------------------------------------------------a
tsim=7;
!cp colamodk2_0.m colamodk2.m
[t0,x]=ode15s('cola4k2',[0 tsim],Xinit); % lambda=0
M10=x(:,42); xb0 = x(:,1);

!cp colamodk2_05.m colamodk2.m
[t1,x]=ode15s('cola4k2',[0 tsim],Xinit); % lambda=0.5
M11=x(:,42); xb1=  x(:,1);

!cp colamodk2_1.m colamodk2.m
[t2,x]=ode15s('cola4k2',[0 tsim],Xinit); % lambda=1
M12=x(:,42); xb2 =  x(:,1);

!cp colamodk2_15.m colamodk2.m
[t3,x]=ode15s('cola4k2',[0 tsim],Xinit); % lambda=1.5
M13=x(:,42); xb3 =  x(:,1);

!cp colamodk2_m5.m colamodk2.m
[t4,x]=ode15s('cola4k2',[0 tsim],Xinit); % lambda=-0.5
M14=x(:,42); xb4 =  x(:,1);

plot(t0,M10,t1,M11,t2,M12,t3,M13,t4,M14); % Inverse response for lambda>1
plot(t0,xb0,t1,xb1,t2,xb2,t3,xb3,t4,xb4); % nothing strange for lambda=0.5

t2 = [-1; t2];
L = [LB(1);LB];
LTv = [ 0;0;LT-LT/1.1;LT-LT/1.1];
tv  = [-1;0; 0;tsim];


% -------------------------------------------------------------
% Figure 5: LV-configuration Simulate a change in reflux rate with V constant 
% -------------------------------------------------------------
% FIRST GO INTO THE FILE colalv.m and increase L by 0.1%
% SAVE THE FILE as colalv1.m simulate:
clear all
load Xinit
tsim=500;
[t,x]=ode15s('colalv1',[0 tsim],Xinit); 
t1 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD = yD(:) - 0.99; dxB = xB(:)-0.01;
% 1n. same without flow dynamics
[t,x]=ode15s('colalv1n',[0 tsim],Xinit);
t1n = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD1n = yD(:) - 0.99; dxB1n = xB(:)-0.01;
% Now do the same change in boilup
[t,x]=ode15s('colalv1v',[0 tsim],Xinit); 
t1v = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyDv = yD(:) - 0.99; dxBv = xB(:)-0.01;
% 1n. same without flow dynamics
[t,x]=ode15s('colalv1vn',[0 tsim],Xinit);
t1vn = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD1vn = yD(:) - 0.99; dxB1vn = xB(:)-0.01;

plot(t1,dxB,'--',t1,dyD,t1n,dxB1n,':',t1n,dyD1n,':');
hold on
plot(t1v,dxBv,'--',t1v,dyDv,t1vn,dxB1vn,':',t1vn,dyD1vn,':'); 
%title('0.1% increase in external flows'); xlabel('Time (min)')
set(gca,'Fontsize',14)
xlabel('Time')
ylabel('Composition')
hl = legend('xB','yD','NoH')
set(hl,'Visible','off','Position', [0.15 0.7 0.2 0.2] )
text(200,1e-3,'Increase in L with V constant');
text(200,-1e-3,'Increase in V with L constant');
%print -deps figlv1
%!mv figlv1.eps ~skoge/latex/work/fig/figlv1.eps
hold off


% -------------------------------------------------------------
% Figure 6: LV-configuration Simulate a change in  internal flows
% -------------------------------------------------------------
% Increse L by 10% and increase V by same amount som D is constant
% SAVE THE FILE as colalv2.m simulate:
clear all
load Xinit
tsim=500;
[t,x]=ode15s('colalv2',[0 tsim],Xinit); 
t2 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD = yD(:) - 0.99; dxB = xB(:)-0.01;


%plot(t2,dxB,t2,dyD); title('10% increase in internal flows') 
% 2n. Do the same without flow dynamics
[t,x]=ode15s('colalv2n',[0 tsim],Xinit); 
t2n = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD2n = yD(:) - 0.99; dxB2n = xB(:)-0.01;

plot(t2,dxB,'--',t2,dyD,t2n,dxB2n,':',t2n,dyD2n,':'); 
%title('10% increase in internal flows'); xlabel('Time (min)') 
set(gca,'Fontsize',14)
xlabel('Time')
ylabel('Composition')
hl = legend('xB','yD','NoH')
set(hl,'Visible','off','Position', [0.15 0.475 0.2 0.2] )
%print -deps figlv2
%!mv figlv2.eps ~skoge/latex/work/fig/figlv2.eps

%-----------------------------------------------
% Linearized model for LV- configuration
%------------------------------------------------

clear all
load Xinit
%The linear model for the LV-configurtation (4 inputs, 2 outputs)
Ls=2.70629; Vs=3.20629; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola_lv_lin',Xinit',[Ls Vs Fs zFs]);
Glvu =  pck(A,B,C,D);
eig(A)
% compare with model with no flow dynamics
[A,B,C,D]=cola_linearize('cola_lvn_lin',Xinit',[Ls Vs Fs zFs]);
eig(A)

%-----------------------------------------------
% Figure 7: Compare linear and nonlinear simulations
%------------------------------------------------
tsim=30;

% Nonlinear of magnitude 0.1%
[t,x]=ode15s('colalv1',[0 tsim],Xinit);
t3 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD3 = (yD(:) - 0.99)/(2.70629*0.001);
YD = log(yD./(1-yD)); dYD3 = (YD(:) - YD(1))/(2.70629*0.001);

% simulate linear step in L of magnitude 1
t3l = linspace(0,tsim,10*tsim); zero = zeros(10*tsim,1);
[Y,X] = step(A,B,C,D,1,t3l);
dyDl = Y(:,1);

% Nonlinear of magnitude 1%
[t,x]=ode15s('colalv1_1',[0 tsim],Xinit);
t3_1 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD3_1 = (yD(:) - 0.99)/(2.70629*0.01);
YD = log(yD./(1-yD)); dYD3_1 = (YD(:) - YD(1))/(2.70629*0.01);

% Nonlinear of magnitude 10%
[t,x]=ode15s('colalv1_10',[0 tsim],Xinit);
t3_10 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD3_10 = (yD(:) - 0.99)/(2.70629*0.1);
YD = log(yD./(1-yD)); dYD3_10 = (YD(:) - YD(1))/(2.70629*0.1);

% Nonlinear of magnitude 50%
[t,x]=ode15s('colalv1_50',[0 tsim],Xinit);
t3_50 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD3_50 = (yD(:) - 0.99)/(2.70629*0.5);
YD = log(yD./(1-yD)); dYD3_50 = (YD(:) - YD(1))/(2.70629*0.5);

% Plotting

plot(t3,dyD3); hold on
plot(t3l,dyDl,'--'); 
plot(t3l,zero,':');
plot(t3_1,dyD3_1); 
plot(t3_10,dyD3_10); 
plot(t3_50,dyD3_50); 
hold off
text(25,0.112,'lin','HorizontalAlignment','left',...
                    'VerticalAlignment','bottom');
text(25,0.104,'+0.1','HorizontalAlignment','left','VerticalAlignment','top');
text(25,0.09,'+1','HorizontalAlignment','left','VerticalAlignment','top');
text(25,0.034,'+10','HorizontalAlignment','left',...
                    'VerticalAlignment','bottom');
text(25,0.007,'+50','HorizontalAlignment','left',...
                    'VerticalAlignment','bottom');

%title('Nonlinearity xD for increase in L');  xlabel('Time (min)')
set(gca,'FontSize',14,'Position',[0.13 0.11 0.775 0.8]); 
set(gcf,'Position',[296 420 560 420],'PaperPosition',[0.25 2.5 8 6])
xlabel('Time'); ylabel('Distillate composition')
%print -deps figlv3
%!mv figlv3.eps ~skoge/latex/work/fig/figlv3.eps

plot(t3,dYD3); hold;
plot(t3_1,dYD3_1);
plot(t3_10,dYD3_10);
plot(t3_50,dYD3_50);
text(25,10.9,'+0.1','HorizontalAlignment','left',...
                     'VerticalAlignment','bottom');
text(25,10.5,'+1','HorizontalAlignment','left','VerticalAlignment','top');
text(25,9.1,'+10','HorizontalAlignment','left',...
                    'VerticalAlignment','top');
text(25,6.8,'+50','HorizontalAlignment','left',...
                    'VerticalAlignment','bottom');
%title('Nonlinearity with log compositions');  xlabel('Time (min)')
set(gca,'FontSize',14,'Position',[0.13 0.11 0.775 0.8]); 
set(gcf,'Position',[296 420 560 420],'PaperPosition',[0.25 2.5 8 6])
xlabel('Time'); ylabel('Distillate composition')
%print -deps figlv3log
%!mv figlv3log.eps  ~skoge/latex/work/fig/figlv3log.eps
hold;

%---------------------------------------------
% Figure 8: Effect of mass flows
%---------------------------------------------
clear all
load Xinit
tsim = 1000;

% Decrease zf by 1 % (from 0.50 to 0.495)
[t,x]=ode15s('colalv4',[0:1:tsim],Xinit);
t4 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD4 = yD(:) - 0.99; dxB4 = xB(:)-0.01;
plot(t4,dyD4,t4,dxB4); 
hold

% Same with MASS reflux constant
[t,x]=ode15s('colalv4mass',[0:1:tsim],Xinit);
%[t,x]=ode15s('dum',[0:1:tsim],Xinit);
%[t,x]=ode23s('colalv4mass',[0:1:tsim],Xinit);
%[t,x]=ode45('colalv4mass',0, 1000, Xinit);
t4mass = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD4mass = yD(:) - 0.99; dxB4mass = xB(:)-0.01;

plot(t4mass,dyD4mass,'--',t4mass,dxB4mass,'--'); 
set(gca,'FontSize',14); xlabel('Time'); ylabel('Composition')
hl = legend('--','MR','-','MF')
set(hl,'Visible','off','Position', [0.15 0.12 0.15 0.18] )
text(400,-0.003,'xB'); text(400,-0.02,'yD');
hold off;
%title('Effect of mass reflux');  xlabel('Time (min)')
axis([0 1000 -0.05 0] )

%print -deps figlv4
%!mv figlv4.eps  ~skoge/latex/work/fig/figlv4.eps

% CHECK:
% check stability by computing eigenvalues along the trajectory
Mwl=30;
Mw0 = Mwl*0.99 + 40*0.01; LTW = 2.70629*Mw0;   % Mass reflux constant
Vs=3.20629; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('colalv4mass_lin',Xinit',[LTW Vs Fs zFs]);
P = eig(A); Ps = sort(P);
[A2,B2,C2,D2]=cola_linearize('colalv4mass_lin',Xinit',[LTW Vs Fs zFs*0.98]);
P2 = eig(A2); P2s = sort(P2);

UstPo = [];
EUstPo = [];
PT = [];
for i=1:size(x,1)
  disp(['Point number: ', int2str(i)] );
  [Ah,Bh,Ch,Dh]=cola_linearize('colalv4mass_lin',...
           x(i,:),[LTW Vs Fs zFs*0.98]);
  Ph = eig(Ah); mPh = max(Ph);
  if mPh >= 0
    disp(['Point ', int2str(i),' is unstable: ' ])
    UstPo  = [UstPo; i];
    EUstPo = [EUstPo; Ph.']
  end
  Phs = sort(Ph);
  PT = [PT; Phs(1:10).'];
end

% save EigsOnTra UstPo EUstPo PT 
%
% NOTE: DOES INDEED BECOME UNSTABLE ALONG THER PATH FOR ML=30.
% END CHECK

%---------------------------------------------
% Figure 9: Nonlinear Simulation with other configurations (close level loops)
%---------------------------------------------

% -------------------------------------------------------------
% Simulate a change in feed rate with level loops open 
% -------------------------------------------------------------
% FIRST GO INTO THE FILE cola.m and change F=1.00 to F=1.01. 
% SAVE THE FILE as cola_F1.m_and simulate:

clear all
load Xinit

[t,x]=ode15s('cola4_F1',[0 500],Xinit); 
t0 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting

% Plot reboiler holdup (Apparent delay of about 1.26 min)
plot(t0,M1); axis([0 3 0.5 0.52]); title('MB: reboiler holdup')
plot(t0,xB); title('xB: reboiler composition')
plot(t0,yD); title('yD: distillate composition')

% LV-configuration
%FIRST GO INTO THE FILE cola_lv.m and change F to 1.01.  Then use 
[t,x]=ode15s('cola_lv_F1',[0 500],Xinit); 
tlv=t; M1lv=x(:,42); xBlv = x(:,1); yDlv = x(:,41); %save data for plot 
% Do the same for other configurations (REMEMBER TO CHANGE F!!!)
% LB-configuration
[t,x]=ode15s('cola_lb_F1',[0 500],Xinit); 
tlb=t; M1lb=x(:,42); xBlb = x(:,1); yDlb = x(:,41);
% Double ratio
[t,x]=ode15s('cola_rr_F1',[0 500],Xinit);
trr=t; M1rr=x(:,42); xBrr = x(:,1); yDrr = x(:,41); 
% DB
[t,x]=ode15s('cola_db_F1',[0 500],Xinit);
tdb=t; M1db=x(:,42); xBdb = x(:,1); yDdb = x(:,41); 
taul=0.063; LB = x(:,43)/taul;
plot(tdb,LB);

% Plot them together
plot(t0,M1,'-',tlv,M1lv,':',tlb,M1lb,'-.',trr,M1rr,'--'); 
title('MB: reboiler (bottom) holdup [mol]');
set(gca,'FontSize',14); xlabel('Time')
%hl = legend('OL','LV','LB','RR')
%set(hl,'Visible','off','Position', [0.15 0.68 0.15 0.25] )

plot(t0,xB,'-',tlv,xBlv,':',tlb,xBlb,'--',trr,xBrr,'--',tdb,xBdb,'--'); 
set(gca,'FontSize',14); xlabel('Time'); %ylabel('xB')
text(350,0.016,'OL','VerticalAlignment','top');
text(225,0.0151,'LV','VerticalAlignment','bottom',...
                    'HorizontalAlignment','right');
text(350,0.01,'RR','VerticalAlignment','bottom');
text(350,0.00695,'LB','VerticalAlignment','bottom');
text(350,0.0032,'DB','VerticalAlignment','bottom');
%print -deps fig5xb
%!mv fig5xb.eps  ~skoge/latex/work/fig/fig5xb.eps
%!mv fig5xb.eps ../latex/fig

plot(t0,yD,'-',tlv,yDlv,':',tlb,yDlb,'--',trr,yDrr,'--',tdb,yDdb,'--'); 
set(gca,'FontSize',14); xlabel('Time'); %ylabel('yD') 
text(350,0.99665,'DB','VerticalAlignment','top',...
                     'HorizontalAlignment','left');
text(350,0.9924,'OL','VerticalAlignment','bottom');
text(150,0.9918,'LV','VerticalAlignment','top');
text(350,0.9900,'RR','VerticalAlignment','top');
text(350,0.9849,'LB','VerticalAlignment','bottom');

%print -deps fig5yd
%!mv fig5yd.eps  ~skoge/latex/work/fig/fig5yd.eps
%!mv fig5yd.eps ../latex/fig

%----------------------------------------------------------------
% Figure 12: Comparisons with the perfect operator. One-point control
%----------------------------------------------------------------

% 1. Feedback control one point. No meas. delay
clear all
L0 = 2.70629;
V0 = 3.20629;
TL = [0    L0;
      1000 L0];
TV = [0    V0;
      9.9  V0;
      10   V0*1.01
      1000 V0*1.01];
colas_PItop;
%----------------------------------------------------------
% IMPORTANT: YOU MUST NOW Start simulation in simulink
%----------------------------------------------------------
t1cr=t; M1cr = y4; xBcr = y2-0.01; yDcr = y1-0.99; 

% 2. Perfect operator at steady-state, use the value of L from 
% feedback control save as cola_lv_op.m
load Xinit
tsim = 500;
[t0,x]=ode15s('cola_lv_op',[0 tsim],Xinit); 
t0op = t0; M1op=x(:,42); xBop = x(:,1)-0.01; yDop = x(:,41)-0.99; 

close all

%PLOT TOP COMPOSITION
plot(t1cr,xBcr,'--',t0op,xBop,'-',t1cr,yDcr,'--',t0op,yDop,'-')
set(gca,'FontSize',14); xlabel('Time'); %ylabel('Compsitions')
hl = legend('--','LC','-','OP')
set(hl,'Visible','off','Position', [0.2 0.45 0.15 0.17] )
text(400,-0.79e-3,'xB')
text(400,1.2e-4,'yD')
%print -deps fig6c
%!mv fig6c.eps ~skoge/latex/work/fig/fig6.eps
%!mv fig6c.eps fig6c.eps ../latex/fig

% PLOT REFLUX:  u1 is L
DL = TL(1,2)+u1(size(u1,1),1)
plot(t1cr,u1+TL(1,2),'--',[0 9.9 10 tsim],[TL(1,2) TL(1,2) DL DL],'-')
set(gca,'FontSize',14); xlabel('Time'); %ylabel('Control value')
hl = legend('--','LC','-','OP')
set(hl,'Visible','off','Position', [0.2 0.45 0.15 0.17] )
%axis([0 tsim 0 0.05])

%print -deps fig6u
%!mv fig6u.eps  ~skoge/latex/work/fig/fig6.eps
%!mv fig6u.eps  ../latex/fig


%-------------------------------------------------------------
% Figure 13:  Same with two-point control 
%-------------------------------------------------------------

% 1. First two-point control
clear all 
close all
L0 = 2.70629;
V0 = 3.20629;
TL = [0    L0;
      1000 L0];
TV = [0    V0;
      1000 V0];
rxB = [0    0.01;
       1000 0.01];
ryD = [0    0.99;
       9.9  0.99;
       10   0.995;
       1000 0.995];
colas_PItu1_nodelay   % with no time delay in measurements
% Simulation start, the simulation takes abot 2-4 min.
% Then:
t1cr=t; M1cr = y4; xBcr = y2-0.01; yDcr = y1-0.99; 

% 2. Now do the "open-loop change (Perfect operator)
load Xinit
tsim = 500;
[t,x]=ode15s('cola_lv_op2',[0 tsim],Xinit); 
t0op = t; M1op=x(:,42); xBop = x(:,1)-0.01; yDop = x(:,41)-0.99; 

close all
plot(t1cr,xBcr,'--',t1cr,yDcr,'--',t0op,yDop,'-',...
       t0op,xBop,'-',[0 tsim],[5e-3 5e-3],':');
set(gca,'FontSize',14); xlabel('Time'); %ylabel('Compsitions')
hl = legend('--','LC','-','OP')
set(hl,'Visible','off','Position', [0.3 0.53 0.1 0.17] )
text(400,-0.7e-3,'xB')
text(400,4e-3,'yD')
%print -deps fig7c
%!mv fig7c.eps  ~skoge/latex/work/fig/fig7c.eps

% To plot the inputs (L and V),

DL = TL(1,2)+u1(size(u1,1),1)
DL = TL(1,2)+u1(size(u1,1),1)
DV = TV(1,2)+u2(size(u1,1),1)
plot(t1cr,u1+TL(1,2),'--',t1cr,u2+TV(1,2),'--',...
     [0 9.9 10 tsim],[TL(1,2) TL(1,2) DL DL],'-',...
     [0 9.9 10 tsim],[TV(1,2) TV(1,2) DV DV],'-')
set(gca,'FontSize',14); xlabel('Time'); %ylabel('Control value')
hl = legend('--','LC','-','OP')
set(hl,'Visible','off','Position', [0.2 0.53 0.1 0.17] )
%axis([0 tsim 0 0.05])
text(400,DL+0.04,'L');
text(400,DV-0.02,'V','VerticalAlignment','top');

%print -deps fig7u
%!mv fig7u.eps  ~skoge/latex/work/fig/fig7u.eps


%------------------------------------------------------------
% Simulate same with time delay (not much effect; see also Fig. 18)
%------------------------------------------------------------

clear all 
colas_PI_dist
L0 = 2.70629;
V0 = 3.20629;
TL = [0    L0;
      1000 L0];
TV = [0    V0;
      1000 V0];
TF = [0    1;
      1000 1];
TzF= [0     0.5;
      1000  0.5];
rxB = [0    0.01;
%       99.9, 0.01;
%       100  0.005;
%       1000 0.005];
       1000 0.01];
ryD = [0    0.99;
       9.9  0.99;
       10   0.995;
       1000 0.995];
% Simulation start, the simulation takes abot 2-4 min.
% Then:
t1=t; M1 = y4; xB1 = y2-0.01; yD1 = y1-0.99; 
L1 = u1; V1 = u2;
close all
tsim=200;
figure(1)
plot(t1,xB1,'--',t1,yD1,'-',[0 tsim],[5e-3 5e-3],':',...
                           [0 tsim],[0 0],':');
text(150,0.00001,'xB','VerticalAlignment','bottom')
text(150,0.005-0.00001,'yD','VerticalAlignment','top');
set(gca,'FontSize',14); xlabel('Time'); ylabel('Compsitions')
axis([0 tsim -0.006 0.006])
%print -deps fig9tu1
%!mv fig9tu1.eps  ~skoge/latex/work/fig/fig9tu1.eps

%------------------------------------------------------------
% Figure 18: Simulation with disturbances and time delay
%------------------------------------------------------------
clear all 
L0 = 2.70629;
V0 = 3.20629;
TL = [0    L0;
      1000 L0];
TV = [0    V0;
      1000 V0];
rxB = [0    0.01;
       1000 0.01];
ryD = [0    0.99;
       199.9  0.99;
       200   0.995;
       1000  0.995];
TF = [0    1;
      9.9  1;
      10   1.2
      1000 1.2];

TzF= [0     0.5;
      99.9  0.5
      100   0.6
      1000  0.6];

colas_PI_dist
% Simulation start, the simulation takes abot 2-4 min.
% Then:
t1=t; M1 = y4; xB = y2-0.01; yD = y1-0.99; 
L = u1; V = u2;

subplot(2,1,1)
plot(t1,yD,t1,xB,'--')
set(gca,'FontSize',14); xlabel('Time'); ylabel('Compsitions')
text(30,2.2e-3,'xB','VerticalAlignment','bottom')
text(30,-3.2e-3,'xD','VerticalAlignment','top')
axis([0 300 -0.006 0.006])
%print -deps dist5yn
%!mv dist5yn.eps  ~skoge/latex/work/fig/dist5yn.eps



% -----------------------------------------------------
% Figure 10: Effect of slow level tuning for DV-configuration 
% -----------------------------------------------------

clear all 
load Xinit

tsim = 400;
% Increase V by 1 % (from 0.50 to 0.495)
% First for condenser level gain K=10
[t,x]=ode15s('coladv8_10',[0 tsim],Xinit);
t8 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD8 = yD(:) - 0.99; dxB8 = xB(:)-0.01;
plot(t8,dyD8,t8,dxB8,'--');
hold on

% Condenser level gain K=1
[t,x]=ode15s('coladv8_1',[0 tsim],Xinit);
t8 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD8 = yD(:) - 0.99; dxB8 = xB(:)-0.01;
plot(t8,dyD8,t8,dxB8,'--');

% Condenser level gain K=0.1
[t,x]=ode15s('coladv8_01',[0 tsim],Xinit);
t8 = t; M1=x(:,42); xB = x(:,1); yD = x(:,41); % Save the data for plotting
dyD8 = yD(:) - 0.99; dxB8 = xB(:)-0.01;
plot(t8,dyD8,t8,dxB8,'--');
hold off

%title('Effect of level tuning for DV');  
xlabel('Time')
axis([0 tsim -20e-4 6e-4])
text(50,3.3e-4,'K=10')
text(50,-5.4e-4,'K=10')
text(150,2.7e-4,'K=1','VerticalAlignment','top')
text(150,-6.1e-4,'K=1','VerticalAlignment','top')
text(250,-0.8e-4,'K=0.1','VerticalAlignment','top')
text(250,-10e-4,'K=0.1','VerticalAlignment','top')
hl = legend('yD','xB')
set(hl,'Visible','off','Position', [0.75 0.1 0.12 0.2] )

%print -deps figdv8
%!mv figdv8.eps  ~skoge/latex/work/fig/figdv8.eps
%!mv figdv8.eps ../latex/fig

%------------------------------------------------------------
%  Linear analysis (Freq.dep. RGA etc.)
% ---------------------------------------------------------

clear all
load Xinit

% LV-configurations
Ls=2.706; Vs=3.206; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola_lv_lin',Xinit',[Ls Vs Fs zFs]);
Glvu =  pck(A,B,C,D);  % Model has 4 inputs and 2 outputs
% Scale the model:
Si = diag([1 1 .2 .1]); So = diag([1/0.01 1/0.01]); %scalings
Glv = mmult (So, Glvu, Si);           

%For the double ratio (L/D-V/B) configuration (with 4 inputs and 2 outputs):
R1s = 2.70629/0.5; R2s=3.20629/0.5; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola_rr_lin',Xinit',[R1s R2s Fs zFs]);
Grru =  pck(A,B,C,D);
Si = diag([1 1 .2 .1]); So = diag([1/0.01 1/0.01]); %scalings
Grr = mmult (So, Grru, Si);           

% We could alternatively start directly from the linearized 4x4 model
% I'll do that for the rest
Ls=2.706; Vs=3.206; Ds=0.5; Bs=0.5; Fs=1.0; zFs=0.5;
[A,B,C,D]=cola_linearize('cola4_lin',Xinit',[Ls Vs Ds Bs Fs zFs]);
G4u =  pck(A,B,C,D); 
Si = diag([1 1 1 1 .2 .1]); So = diag([1/0.01 1/0.01 1 1]); %scalings
G4 = mmult (So, G4u, Si);               % scaled 82 state model

% Use sysic to generate DV-configuration
Kl = 10; Kb = 10; % P-controllers for levels (bandwidth = 10 rad/min)
systemnames = 'G4 Kl Kb';
inputvar = '[D(1); V(1); d(2)]'; outputvar = '[G4(1);G4(2)]';
input_to_G4 = '[Kl; V; D; Kb; d ]';
input_to_Kl = '[G4(3)]'; input_to_Kb = '[G4(4)]';
sysoutname ='Gdv'; cleanupsysic='yes'; sysic;

% DB-configuration
Kl = 10; Kv = 10;
systemnames = 'G4 Kl Kv';
inputvar = '[D(1); B(1); d(2)]'; outputvar = '[G4(1);G4(2)]';
input_to_G4 = '[Kl; Kv; D; B; d ]';
input_to_Kl = '[G4(3)]'; input_to_Kv = '[G4(4)]';
sysoutname ='Gdb'; cleanupsysic='yes'; sysic;

    % Now do the same as a function of freq
w = sort([logspace(-4,0,161) 1.1:0.1:10 2.02:0.02:3]);     

% For each configuration do the following
G = sel(Glv,':',[1 2]); Gd= sel(Glv,':',[3 4]);
Glvf = frsp(G,w); Gdlvf=frsp(Gd,w);  % Frequency repsonse of G and Gd
G = sel(Gdv,':',[1 2]); Gd= sel(Gdv,':',[3 4]);
Gdvf = frsp(G,w); Gddvf=frsp(Gd,w);  % Frequency repsonse of G and Gd
G = sel(Gdb,':',[1 2]); Gd= sel(Gdb,':',[3 4]);
Gdbf = frsp(G,w); Gddbf=frsp(Gd,w);  % Frequency repsonse of G and Gd
G = sel(Grr,':',[1 2]); Gd= sel(Grr,':',[3 4]);
Grrf = frsp(G,w); Gdrrf=frsp(Gd,w);  % Frequency repsonse of G and Gd

%-------------------------
% Figure 11: Diagonal RGA-elements
%-------------------------
rlv = sel(vrga(Glvf),1,1);
rdv = sel(vrga(Gdvf),1,1);
rdb = sel(vrga(Gdbf),1,1);
rrr = sel(vrga(Grrf),1,1);

vplot('liv,lm',rlv,rdv,'--',rdb,'--',rrr,'--',1,':');  
axis([0.0001,10,.1,100])
text(.00015,35,'LV','VerticalAlignment','bottom')
text(.00015,17.2,'DB','VerticalAlignment','bottom')
text(.00015,0.44,'DV','VerticalAlignment','top')
text(.00015,3.25,'L/D V/B','VerticalAlignment','top')
xlabel('Frequency [rad/min]')
ylabel('Magnitude of diagonal RGA-element')
%print -deps figrga
%!mv figrga.eps  ~skoge/latex/work/fig/figrga.eps
%!mv figrga.eps ../latex/fig

%--------------------------------
% Figure 14: Open-loop effect of disturbances
%--------------------------------


% Consider effect of feed rate on top composition
dlv = sel(Gdlvf,1,[1 ]);
ddv = sel(Gddvf,1,[1 ]);
ddb = sel(Gddbf,1,[1 ]);
drr = sel(Gdrrf,1,[1 ]);

vplot('liv,lm',dlv,ddv,'-.',ddb,'--',drr,1,':');
axis([0.001,1,.01,100])
text(.002,65,'DB')
text(.002,9,'LV, DV')
text(.004,0.016,'L/D V/B')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
%print -deps figdf0
%!mv figdf0.eps  ~skoge/latex/work/fig/figdf0.eps
%!mv figdf0.eps  ../latex/fig

% Consider effect of feed composition on top composition
dlv = sel(Gdlvf,1,[2 ]);
ddv = sel(Gddvf,1,[2 ]);
ddb = sel(Gddbf,1,[2 ]);
drr = sel(Gdrrf,1,[2 ]);

vplot('liv,lm',dlv,ddv,'-.',ddb,'--',drr,1,':');
axis([0.001,1,.01,100])
text(.0015,8,'LV, DV, DB, L/D V/B (all configurations)',...
      'VerticalAlignment', 'bottom')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
%print -deps figdz0
%!mv figdz0.eps  ~skoge/latex/work/fig/figdz0.eps
%!mv figdz0.eps  ../latex/fig

%------------------------------------------------
% Figure 15: Close bottom composition loop (partial control)
%------------------------------------------------
gd1 = sel(Gdlvf,1,[1 2]); gd2 = sel(Gdlvf,2,[1 2]); 
g12 = sel(Glvf,1,2);  g22 = sel(Glvf,2,2);
Pdlv = msub(gd1,mmult(g12,minv(g22),gd2));
gd1 = sel(Gddvf,1,[1 2]); gd2 = sel(Gddvf,2,[1 2]); 
g12 = sel(Gdvf,1,2);  g22 = sel(Gdvf,2,2);
Pddv = msub(gd1,mmult(g12,minv(g22),gd2));
gd1 = sel(Gddbf,1,[1 2]); gd2 = sel(Gddbf,2,[1 2]); 
g12 = sel(Gdbf,1,2);  g22 = sel(Gdbf,2,2);
Pddb = msub(gd1,mmult(g12,minv(g22),gd2));
gd1 = sel(Gdrrf,1,[1 2]); gd2 = sel(Gdrrf,2,[1 2]); 
g12 = sel(Grrf,1,2);  g22 = sel(Grrf,2,2);
Pdrr = msub(gd1,mmult(g12,minv(g22),gd2));

% Feed rate disturbance
plv = sel(Pdlv,1,1);
pdv = sel(Pddv,1,1);
pdb = sel(Pddb,1,1);
prr = sel(Pdrr,1,1);

vplot('liv,lm',plv,pdv,'-.',pdb,'--',prr,1,':');
axis([0.001,1,.01,100])
text(.0015,18,'DV, DB','VerticalAlignment','bottom')
text(.0015,1.3,'LV','VerticalAlignment','bottom')
text(.004,0.016,'L/D V/B','VerticalAlignment','bottom')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
%print -deps figdf1
%!mv figdf1.eps  ~skoge/latex/work/fig/figdf1.eps
%!mv figdf1.eps  ../latex/fig

% Feed composition disturbance
plv = sel(Pdlv,1,2);
pdv = sel(Pddv,1,2);
pdb = sel(Pddb,1,2);
prr = sel(Pdrr,1,2);

vplot('liv,lm',plv,pdv,'-.',pdb,'--',prr,1,':');
axis([0.001,1,.01,100])
text(.0015,0.0325,'LV','VerticalAlignment','bottom')
text(.0015,19,'DV, DB','VerticalAlignment','bottom')
text(.0015,1.3,'L/D V/B','VerticalAlignment','bottom')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
%print -deps figdz1
%!mv figdz1.eps  ~skoge/latex/work/fig/figdz1.eps
%!mv figdz1.eps  ../latex/fig

%-------------------------------------------
% Figure 16: Two-point control (CLDG)
%-------------------------------------------

gdiag=vdiag(vdiag(Glvf)); prga = mmult(gdiag,minv(Glvf)); 
cldglv=mmult(prga,Gdlvf);
gdiag=vdiag(vdiag(Gdvf)); prga = mmult(gdiag,minv(Gdvf)); 
cldgdv=mmult(prga,Gddvf);
gdiag=vdiag(vdiag(Gdbf)); prga = mmult(gdiag,minv(Gdbf)); 
cldgdb=mmult(prga,Gddbf);
gdiag=vdiag(vdiag(Grrf)); prga = mmult(gdiag,minv(Grrf)); 
cldgrr=mmult(prga,Gdrrf);

% Feed rate disturbance on xD
clv = sel(cldglv,1,1);
cdv = sel(cldgdv,1,1);
cdb = sel(cldgdb,1,1);
crr = sel(cldgrr,1,1);

vplot('liv,lm',clv,cdv,'-.',cdb,'--',crr,1,':');
axis([0.001,1,.01,100])
text(.0015,45,'LV','VerticalAlignment','bottom')
text(.004,15,'DB','VerticalAlignment','bottom')
text(.0015,8.5,'DV','VerticalAlignment','bottom')
text(.0015,0.0185,'L/D V/B','VerticalAlignment','bottom')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
%print -deps figdf2
%!mv figdf2.eps  ~skoge/latex/work/fig/figdf2.eps
%!mv figdf2.eps  ../latex/fig

% Feed composition disturbance on xD
clv = sel(cldglv,1,2);
cdv = sel(cldgdv,1,2);
cdb = sel(cldgdb,1,2);
crr = sel(cldgrr,1,2);

vplot('liv,lm',clv,cdv,'-.',cdb,'--',crr,1,':');
axis([0.001,1,.01,100])
text(.0015,1.15,'LV','VerticalAlignment','bottom')
text(.0015,35.5,'DB','VerticalAlignment','bottom')
text(.0015,8.5,'DV','VerticalAlignment','bottom')
text(.0015,4,'L/D V/B','VerticalAlignment','top')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
%print -deps figdz2
%!mv figdz2.eps  ~skoge/latex/work/fig/figdz2.eps
%!mv figdz2.eps  ../latex/fig/figdz2.eps

% Now evaluate the CLDGs for the BOTTOM:
% Feed rate disturbance on xB
clv = sel(cldglv,2,1);
cdv = sel(cldgdv,2,1);
cdb = sel(cldgdb,2,1);
crr = sel(cldgrr,2,1);

vplot('liv,lm',clv,cdv,'-.',cdb,'--',crr,1,':');
axis([0.001,1,.01,100])
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
text(.0065,43,'LV','VerticalAlignment','bottom')
text(.0065,10.5,'DB','VerticalAlignment','bottom')
text(.0065,1.15,'DV','VerticalAlignment','bottom')
text(.002,0.046,'L/D V/B','VerticalAlignment','bottom')
%print -deps figdf2b
%!mv figdf2b.eps  ~skoge/latex/work/fig/figdf2b.eps
%!mv figdf2b.eps  ../latex/fig/figdf2b.eps

% Feed composition disturbance on xB
clv = sel(cldglv,2,2);
cdv = sel(cldgdv,2,2);
cdb = sel(cldgdb,2,2);
crr = sel(cldgrr,2,2);

vplot('liv,lm',clv,cdv,'-.',cdb,'--',crr,1,':');
axis([0.001,1,.01,100])
text(.0015,9.0,'LV','VerticalAlignment','bottom')
text(.0015,36,'DB','VerticalAlignment','bottom')
text(.0015,0.125,'DV','VerticalAlignment','bottom')
text(.0015,6.2,'L/D V/B','VerticalAlignment','top')
xlabel('Frequency [rad/min]')
ylabel('Magnitude of xB (scaled)')
%print -deps figdz2b
%!mv figdz2b.eps  ~skoge/latex/work/fig/figdz2b.eps
%!mv figdz2b.eps  ../latex/fig/figdz2b.eps


%------------------------------------------------------------
% Figure 17: Two point control LV configuration, with given controllers
%------------------------------------------------------------

k1 = nd2sys([3.76 1],[3.76 1e-8],26.1/100);
k2 = nd2sys([3.31 1],[3.31 1e-8],-37.5/100);

Kdiag = daug(k1,k2);

Kdiagjw = frsp(Kdiag,w);
Llvjw = mmult(Glvf,Kdiagjw);
gdiaglv=vdiag(vdiag(Glvf)); prgalv = mmult(gdiaglv,minv(Glvf)); 
cldglv=mmult(prgalv,Gdlvf);

figure(1)
vplot('liv,lm', sel(Llvjw,1,1),'--',sel(prgalv,1,1),sel(prgalv,1,2),...
                sel(cldglv,1,1),sel(cldglv,1,2),1,':')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
axis([1e-3 10 1e-2 1e3])
text(1e-2,3e2,'l1','VerticalAlignment','bottom');
text(2.5,3.4,'p11','VerticalAlignment','bottom');
text(3.5,0.28,'p12','VerticalAlignment','top');
text(1e-2,20,'c11','VerticalAlignment','bottom');
text(0.1,0.14,'c12','VerticalAlignment','bottom');
%print -deps fPIfl1
%!mv fPIfl1.eps  ~skoge/latex/work/fig/fPIfl1.eps


figure(2)
vplot('liv,lm', sel(Llvjw,2,2),'--',sel(prgalv,2,1),sel(prgalv,2,2),...
                sel(cldglv,2,1),sel(cldglv,2,2), 1,':')
xlabel('Frequency [rad/min]')
ylabel('Magnitude')
axis([1e-3 10 1e-2 1e3])
text(2e-2,1.5e2,'l2','VerticalAlignment','bottom');
text(0.6e-2,1.5e2,'l2d','VerticalAlignment','bottom');
text(2e-3,62,'c21','VerticalAlignment','bottom');
text(2.5,5.7,'p21','VerticalAlignment','bottom');
text(1.7,0.6,'p22','VerticalAlignment','top');
text(2e-3,9,'c22','VerticalAlignment','bottom');
%print -deps fPIfl2
%!mv fPIfl2.eps  ~skoge/latex/work/fig/fPIfl2.eps


% fjerne figurtekst: close all
% Etterpaa figurer: ps2frag *.eps


