%Example 2.4 , PI-control

% Plant (2.26), G(s)=3(-2s+1)/(5s+1)(10s+1)
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_4.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

clear
s=tf('s');
G=3*(-2*s+1)/(5*s+1)/(10*s+1);

% frequency response of G
points=301;
w=logspace(-3,1,points);

[Ku,pm,Wu,wcp]=margin(G);    % Ku=2.5 and Wu=15.2
Pu=2*pi/Wu;Kc=Ku/2.2;ti=Pu/1.2; %1.137,12.7, Ziegler-Nichols' tuning
K=Kc*(1+1/ti/s);  % PI-Controller
L = G*K;
S = inv(1+L);
T = 1-S;

%Figure 2.7
% Step response
figure(1);clf;
time=0:0.5:50;
step(T,'b',K*S,'g--',time)
axis([0 50 -0.5 2.5]);
text(40,1.3,'OUTPUT (y)');text(40,0.0,'INPUT (u)');
xlabel('Time');

%Figure 2.8
% Frequency responses of L, S, T
w=logspace(-2,1,41);
figure(2);clf;
[mag1,pha1]=bode(L,w);
[mag2,pha2]=bode(S,w);
[mag3,pha3]=bode(T,w);
subplot(2,1,1)
loglog(w,mag1(:),w,mag2(:),w,mag3(:))
ylabel('Magnitude')
hold %on
plot(w,w./w,':')
hold %off
subplot(2,1,2)
semilogx(w,pha1(:),w,pha2(:),w,pha3(:),w,w*0,':')
axis([.01,10,-300,120])
xlabel('Frequency [rad/sec]'),ylabel('Phase')

% Find GM, PM and the peak values of $S$ and $T$ 
[Gm,Pm,W180,Wc]=margin(L)
[norm_S,wS]=norm(S,inf,1e-4)  % find peak is 3.92 at w=0.27
[norm_T,wT]=norm(T,inf,1e-4)  % find peak is 3.35 at w=0.26
