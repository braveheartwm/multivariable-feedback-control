% Section 5.15.2 Application: Room heating
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: RoomHeat.m,v 1.4 2004/02/05 12:23:14 vidaral Exp $

clear all; clear all;
set(cstprefs.tbxprefs,'MagnitudeUnits','log','MagnitudeScale','log');

% Plant (5.111): G(s)=20/(1000s+1), Gd(s)=10/(1000s+1)
G=tf(20,[1000 1]);
Gd=tf(10,[1000 1]);
w=logspace(-5,0,61);

%Figure 5.19
figure(1)
bodemag(G,'-',Gd,'-',tf(1),'--',w); hold on; plot([1e-2 1e-2],[1.5e-2 1],':k'); hold off;
axis([1.e-5 3.e-1 1.e-1 40])
xlabel('Frequency [rad/s]');ylabel('Magnitude');
text(1.e-4,25,'|G|');text(1.e-4,7,'|Gd|'); text(1e-2,2e-1,'\omega_d');

%PID controller: K(s)=0.4(200s+1)(60s+1)/200s(0.1*60s+1)
kc=0.4;taui=200;taud=60;

%measurement delay :
thetam=100;

%Simulation using SIMULINK:
r=0;d=1; %Disturbance step response.

rhmodel;    %display block diagram
[td,xd,yd]=sim('rhmodel',1000, [],[]);
stepydpid=yd(:,1);
stepudpid=yd(:,2);

%Figure 5.20
figure(2);subplot(2,2,1);
plot(td,stepydpid,td,stepudpid);
axis([0 1000 -0.8 1.25]);
text(400,0.4,'y(t)');text(400,-0.4,'u(t)');

%Reference step response
r=1;d=0;

[ts,xs,ys]=sim('rhmodel',1000,[],[]);
stepyrpid=ys(:,1);
stepurpid=ys(:,2);

%Figure 4.21
subplot(2,2,2);
z1(1:length(ts))=3;
plot(ts,stepyrpid,ts,stepurpid,ts,z1,':k');
axis([0 1000 -0.1 3.5]);
text(200,2,'y(t)');text(400,0.5,'u(t)');text(500,2.8,'r(t)');
