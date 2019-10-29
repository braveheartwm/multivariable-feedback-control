%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3.17 Simple 2x2 plant with RHP zero. (Figure 3.9 & 3.10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl3_8.m,v 1.4 2004/02/09 12:01:24 vidaral Exp $
clear all; close all;
G=tf(1,[0.2 1])*tf(1,[1 1])*[1 1;tf([2 1],[0 1]) 2];
G=ss(G);
poles=pole(G) 	% 4 poles, minimum realization
z=zero(G)	% zero at z=1.5

% Analyze zero direction
% G05 = frsp(G,-0.5i)
G05=evalfr(G,z);
[u,s,v]=svd(G05);  % look at last column in u and v

%Open-loop simulation
figure(1)                   % Figure 3.11
u1=[1*ones(101,1) 0*ones(101,1)];
u2=[0*ones(101,1) 1*ones(101,1)];
u3=[1*ones(101,1) -1*ones(101,1)];
t=[0:0.1:10];
y1=lsim(G,u1,t);
y2=lsim(G,u2,t);
y3=lsim(G,u3,t);
subplot(231);plot(t,y1); text(3,1.2,'y_2'); text(3,0.8,'y_1');
subplot(232);plot(t,y2); text(3,0.8,'y_1'); text(3,1.7,'y_2');
subplot(233);plot(t,y3); text(3,0.1,'y_1'); text(1.5,0.9,'y_2');

% H-INFINITY S-KS design

% Design 1: Require  wb, slope 1 and peak less than M
Mi=1.5; wbi=z/2; Ai=1e-4;
wpi = tf([1/Mi wbi], [1 wbi*Ai]); % Same weight both outputs
Wp=[wpi 0;0 wpi]; Wu=1.0*eye(2);
% systemnames = 'G Wp Wu';
% inputvar = '[ r(2); u(2)]';
% outputvar = '[Wp; Wu; r-G]';
% input_to_G = '[u]';
% input_to_Wp = '[r-G]';
% input_to_Wu = '[u]';
% sysoutname = 'P';
% cleanupsysic = 'yes';
% sysic;
% nmeas=2; nu=2; gmn=0.667; gmx=20; tol=0.001;
% [khinf1,cl1,hinf1] = hinfsyn(P,nmeas,nu,'GMAX',gmx,'GMIN',gmn,'TOLGAM',tol); 

[khinf1, hinf1, gopt1] = mixsyn(G,Wp,Wu,[]);

K1=khinf1;
disp(sprintf('H_inf norm design 1: %0.5g',gopt1)); % Optimal hinf norm

S1=inv(eye(2)+G*K1);
T1=eye(2)-S1;
omega=logspace(-3,2,121);
% Simulation of bad (y1) and good (y2) direction
yd1=lsim(T1,u3,t);

% Design 2: New performance weight
Mi=1.5; wb1=z/2; wb2=50*z; wp1 = tf([1/Mi wb1], [1 wb1*Ai]); wp2=tf([1/Mi wb2],[1 wb2*Ai]); 
Wp=[wp1 0;0 wp2]; 
% systemnames = 'G Wp Wu';
% inputvar = '[ r(2); u(2)]';
% outputvar = '[Wp; Wu; r-G]';
% input_to_G = '[u]';
% input_to_Wp = '[r-G]';
% input_to_Wu = '[u]';
% sysoutname = 'P';
% cleanupsysic = 'yes';
% sysic;
% nmeas=2; nu=2; gmn=0.667; gmx=20; tol=0.001;
% [khinf2,cl2,hinf2] = hinfsyn(P,nmeas,nu,'GMAX',gmx,'GMIN',gmn,'TOLGAM',tol);

[khinf2, ghinf2, gopt2] = mixsyn(G,Wp,Wu,[]);

K2=khinf2;
disp(sprintf('H_inf norm design 2: %0.5g',gopt2));

S2=inv(eye(2)+G*K2);
T2=eye(2)-S2;
yd2=lsim(T2,u3,t);

% Plot Figure 2
figure(2);            %Figure 3.12
set(cstprefs.tbxprefs,'MagnitudeUnits','abs','MagnitudeScale','log');
subplot(121);sigma(S1,'b-',S2,'g--',tf(1),':k',omega);
text(1e0,3e-2,'\sigma_{min}(S)'); text(1e0,6e-1,'\sigma_{min}(S)');
text(5e-2,9e-1,'\sigma_{max}(S)'); text(0.9e1,2.3e0,'\sigma_{max}(S)');
axis([0.01,100,0.001,10]);
subplot(122);
plot(t,yd1,'b-',t,yd2,'g--',t,ones(1,length(t)),':')
xlabel('Time');
text(4,-1.2,'y_2');text(4,0.8,'y_1');
axis([0,5,-2,2.5]);
