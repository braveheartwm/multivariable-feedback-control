% Example 5.9: H_inf design for plant with RHP-pole and RHP-zero.
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl5_9.m,v 1.3 2004/02/03 14:12:00 vidaral Exp $

clear all; close all;

z=4; % zero=4
p=1; % pole=1 
G=tf([1 -z],conv([1 -p],[0.1 1])); G=ss(G);

Wu=1;
M=2; wb=1;
Wp=tf([1/M wb],[1 1e-8]); Wp=ss(Wp);

% systemnames  = 'G Wp Wu';
% inputvar     = '[ r(1); u(1)]';
% outputvar    = '[Wp; Wu; r-G]';
% input_to_G   = '[u]';
% input_to_Wp  = '[r-G]';
% input_to_Wu  = '[u]';
% sysoutname   = 'P';
% cleanupsysic = 'yes';
% sysic;
% 
% nmeas=1; nu=1; gmn=0.5; gmx=20; tol=0.001;

[K,clp,gmfin] = mixsyn(G,Wp,Wu,[]);
% [K,clp,gmfin] = hinfsyn(P,nmeas,nu);

K=tf(K);
L = G*K; 
S = 1/(1+L);  T = L/(1+L);
KS = K*S;
if all(eig(T)<0)
   T = balreal(T);
end
if all(eig(KS)<0),
   KS = balreal(KS);
end

norm(S,inf,1e-7)
norm(T,inf,1e-7)
w=logspace(-2, 2, 401);

t=0:0.001:5;
y=step(T,t);
u=step(KS,t);

%Figure 5.13
figure(1)
subplot(2,2,1);
bodemag(S,'-',T,'-',tf(1),':',w); hold on; plot([p p],[2e-2 1],'k--',[z z],[2e-2 1],'k--'); hold off;
axis([.01,100,.01,4]);
text(.05,.17,'|S|');
text(.05,1.5,'|T|');
text(p,1.8e-2,'p');
text(z,1.8e-2,'z');
xlabel('Frequency [rad/s]');
ylabel('Magnitude')

subplot(2,2,2);
plot(t,y,'-',t,u,'--',t,ones(1,length(t)),':k');
axis([0,5,-2.1,2.1]);
xlabel('Time [sec]');
text(4,1.5,'y(t)');
text(4,0.6,'u(t)');
