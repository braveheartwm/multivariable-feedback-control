%Section 5.15.2 Application: Room heating
%
clear all
% Plant (5.111): G(s)=20/(1000s+1), Gd(s)=10/(1000s+1)
G=tf(20,[1000 1]);
Gd=tf(10,[1000 1]);
w=logspace(-5,0,61);
G_f=freqresp(G,w);
Gd_f=freqresp(Gd,w);

for i=1:length(G_f)
    G_f2(i)=norm(G_f(1,1,i));
    Gd_f2(i)=norm(Gd_f(1,1,i));
end

%Figure 5.19
figure(1)
loglog(w,G_f2,w,Gd_f2)
axis([1.e-5,1,1.e-1,40])
xlabel('Frequency [rad/s]');ylabel('Magnitude');
text(1.e-4,25,'|G|');text(1.e-4,7,'|Gd|');

%PID controller: K(s)=0.4(200s+1)(60s+1)/200s(0.1*60s+1)
kc=0.4;taui=200;taud=60;
%measurement delay :
thetam=100;

%Simulation using SIMULINK:
r=0;d=1; %Disturbance step response.

rhmodel;    %display block diagram
[t,x,y]=sim('rhmodel',1000);
% stepydpid=vpck(y(:,1),t);
% stepudpid=vpck(y(:,2),t);

%Figure 5.20
figure(2);subplot(221);
plot(t,y,t,3);
text(400,0.4,'y');text(400,-0.4,'u');
%Reference step response
r=1;d=0;

[t,x,y]=sim('rhmodel',1000);
% stepyrpid=vpck(y(:,1),t);
% stepurpid=vpck(y(:,2),t);

%Figure 4.21
subplot(222);
plot(t,y);
text(400,3,'y');text(400,0.5,'u');

