% Example 5.13: Unstable plant and input constraints
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl5_13.m,v 1.3 2004/02/03 14:12:00 vidaral Exp $
clear all; close all;

%PLANT (5.88)
w=logspace(-2,2,501);
p=1.0;  
g1=tf([1],[1 -p]); dyn=tf(5,[10 1]); G=g1*dyn;
kd=0.5; Gd=tf(kd,[0.2 1.2 1]);

% PID
Kc=0.4;
c1=tf([10 1],[10 1.e-5]); c1=Kc*c1; c2=tf([10 1],[.1 1]); 
c3 = tf(1,[.1 1]); c2=c2*c3;
K=c1*c2;

L = G*K; 

subplot(2,2,1);
bodemag(G,'-',Gd,'-',tf(1),'--',w); hold on; plot([p p],[2e-2 1],':k'); hold off;
axis([1e-2 20 1e-2 1e1])
xlabel('Frequency [rad/s]');
ylabel('Magnitude')
text(p,2e-2,'p');
text(.35,2,'|G|');
text(3,.2,'|Gd|');

[Ac,Bc,Cc,Dc]=ssdata(K);
[Ag,Bg,Cg,Dg]=ssdata(G);
[Ad,Bd,Cd,Dd]=ssdata(Gd);

Expl5_13s;    %display block diagram
[t,x,y]=sim('Expl5_13s', 10, [],[]);

% Data from simulink
time=tout; u1=uout; u2=u2out; y1=yout; y2=y2out;  % WITH CONSTRAINTS
z1(1:length(time))=0;

subplot(2,2,2);
plot(time,y2,'--',time,y1,time,u1,time,u2,'--',time,z1,':');
axis([0,10,-1.5,1.5]);
xlabel('Time [sec]');
text(5,-0.7,'y(t)');
text(6,1.15,'u(t)');
legend('Unconstrained','Constrained');
legend('boxoff');
