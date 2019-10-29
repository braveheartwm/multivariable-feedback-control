%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark 3 in Section 3.7.2, S-KS DESIGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Rem3_372.m,v 1.1 2004/01/21 10:55:10 standber Exp $

G0 = [87.8 -86.4; 108.2 -109.6];
dyn = tf(1,[75 1]); 
G=dyn*G0;

A= 1.e-4; M=2; wb=0.05;
wp1 = tf([1/M wb], [1 wb*A]); wp2=wp1;
Wp= [wp1 0;0 wp2];
Wu = 1.0*eye(2);

% systemnames = 'G Wp Wu'
% inputvar = '[ r(2); u(2)]';
% outputvar = '[Wp; Wu; r-G]';   % NEGATIVE feedback
% input_to_G = '[u]';
% input_to_Wp = '[r-G]';
% input_to_Wu = '[u]';
% sysoutname = 'P';
% cleanupsysic = 'yes';
% sysic;
% nmeas=2; nu=2; gmn=1; gmx=100; tol=0.01;
% [khinf,ghinf,gopt,info] = hinfsyn(P,nmeas,nu,'GMAX',gmx,'GMIN',gmn,'TOLGAM',tol,'DISPLAY','On');

[khinf, hinf, gopt] = mixsyn(G,Wp,Wu,[]);

K=khinf; 

% TIME  simulation
% Nominal
I2=eye(2);
cls=loopsens(G*K,I2);
Kr=tf(1,[5 1]); % 5 min filter on reference change
Tr = cls.Ti*Kr;
u1=[1*ones(1001,1) 0*ones(1001,1)];
t=[0:0.1:100];
y=lsim(Tr,u1,t);
u=lsim(K*cls.Si*Kr,u1,t);

% With 20% uncertainty
Unc = [1.2 0; 0 0.8];
GKu=G*Unc*K;
clsu=loopsens(GKu,I2);
Tru=clsu.Ti*Kr;
yu=lsim(Tru,u1,t);
uu=lsim(K*clsu.Si*Kr,u1,t);

subplot(211);plot(t,y,t,yu,'--');title('OUTPUTS')
axis([0 100 0 4]);
text(60,3.5,'Nominal plant:   solid line');
text(60,3.0,'Perturbed plant: dashed line');
subplot(212);plot(t,u,t,uu,'--');title('INPUTS');
xlabel('TIME (min)');
% Simulation shows that the controller is very sensitive
% to uncertainty.  Large and Long peak in y1 and y2.
