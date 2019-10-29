% Example 5.6 and Example 5.7: Limitations at low and high frequencies (Effect of
% RHP-zeros. Try different sign of gain at high frequency)
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl5_6_7.m,v 1.4 2004/02/03 14:08:37 vidaral Exp $
clear all; close all;

w=logspace(-2,2,141);
G=tf([-1 1],[1 0]);
K=tf(1,[0.05,1]);

Kc=0.2;
GK=G*K*Kc;
S1=1/(1+GK);
T1=GK/(1+GK);
[y1,t1]=step(T1);
z(1:length(w))=1;
z1(1:length(t1))=0;
z2(1:length(t1))=1;

Kc=0.5;
GK=G*K*Kc;
S2=1/(1+GK);
T2=GK/(1+GK);
[y2,t2]=step(T2);

Kc=0.8;
GK=G*K*Kc;
S3=1/(GK+1);
T3=GK/(1+GK);
[y3,t3]=step(T3);

figure(1) %Figure 5.7
subplot(1,2,1); bodemag(S1,'-',S2,'-',S3,'-',tf(1),':',w);
axis([0.01 100 0.01 10]);
text(2.5,8,'Kc=0.8'); text(2.5,3,'Kc=0.5'); text(.03,0.5,'Kc=0.2');
xlabel('Frequency [rad/s]');ylabel('Magnitude |S|');
title('SENSITIVITY FUNCTION');

% Plotting time responses
subplot(1,2,2); plot(t1,y1,t2,y2,t3,y3,t1,z1,':k',t1,z2,':k');
axis([0 5 -2 2]);
text(0.25,1.7,'Kc=0.8'); text(0.8,0.7,'Kc=0.5'); text(2,0.55,'Kc=0.2');
text(3.5,1.3,'Setpoint');
title('STEP RESPONSE');
xlabel('Time [sec]');ylabel('y(t)');

% Effect of rhp-zeros. Try different sign of gain at high frequency
%Plant is -s+1/ s+1. Try controller c = K s with different gains

G=tf([-1 1],[1 1]);

%Controller
K=tf([1 0],conv([0.02 1],[0.05 1]));

% Transient step change
t=0:0.001:0.2; u(1:101)=1; u(102:length(t))=0;

Kc2=-0.1;
GK=G*K*Kc2;
S1=1/(1+GK);
T1=GK/(1+GK);
[y1,t1]=lsim(T1,u,t);
z(1:length(w))=1;
z1(1:length(t1))=0;
z2(1:length(t1))=1;

Kc2=-0.5;
GK=G*K*Kc2;
S2=1/(1+GK);
T2=GK/(1+GK);
[y2,t2]=lsim(T2,u,t);

Kc2=-0.9;
GK=G*K*Kc2;
S3=1/(1+GK);
T3=GK/(1+GK);
[y3,t3]=lsim(T3,u,t);

figure(2),     %Figure 5.8
subplot(1,2,1);
bodemag(S1,'-',S2,'-',S3,'-',tf(1),':',w);
axis([.01,100,.01,10]);
text(1.5,5,'Kc=0.9'); text(2,2,'Kc=0.5'); text(15,.3,'Kc=0.1');
xlabel('Frequency [rad/s]');ylabel('Magnitude |S|');
title('SENSIVITY FUNCTION');

% Plotting time response
subplot(1,2,2),plot(t1,y1,t2,y2,t3,y3,t,u,':k');
axis([0 0.14 -0.3 1.2]);
text(0.06,0.95,'Kc=0.9'); text(.06,0.65,'Kc=0.5'); text(.06,.3,'Kc=0.1');
text(0.06,1.1,'Setpoint');
title('STEP RESPONSE');
xlabel('Time [sec]');ylabel('y(t)');


% Exercício 5.8

U1=K*S1;
[u1,t1]=lsim(U1,u,t);

U2=K*S2;
[u2,t2]=lsim(U2,u,t);

U3=K*S3;
[u3,t3]=lsim(U3,u,t);

figure(3)
plot(t1,u1,'g:',t2,u2,'r--',t3,u3,'b')
text(0.06,.8,'K_c=0.9'); text(.06,2.1,'K_c=0.5'); text(.06,5,'K_c=0.1');
xlabel('Time');ylabel('u(t)');
%axis([0 0.14 -0.3 1.2]);


