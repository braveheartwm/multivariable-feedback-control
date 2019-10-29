%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 5.15.3 Application: Neutralization process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plant (5.114): G(s)=-2kd/(1+tauh*s), Gd(s)=kd/(1+tauh*s), kd=2.5e6, tauh=1000
kd=2.5e6;tauh=1000;
% G=nd2sys(1,[tauh 1],-2*kd);
% Gd=nd2sys(1,[tauh 1],kd);
G=tf(-2*kd,[tauh 1]);
Gd=tf(kd,[tauh 1]);
w=logspace(-4,4,61);
% G_f=frsp(G,w);
% Gd_f=frsp(Gd,w);
G_f=freqresp(G,w);
Gd_f=freqresp(Gd,w);

for i=1:length(G_f)
    G_f2(i)=G_f(1,1,i);
    Gd_f2(i)=Gd_f(1,1,i);
end

%Figure 5.22
figure(3);clf
wd = [1 2500; 0.05 2500; 2 Inf]
wd2 = [0.02 2500; 0.01 2500; 2 Inf]
% vplot('liv,lm',Gd_f,G_f,1,'w:',wd,'w:',wd2,'w:');
loglog(w,G_f2,w,Gd_f2,2500,logspace(-2,0,61),w,1)
xlabel('FREQUENCY [RAD/S]');ylabel('MAGNITUDE');
text(0.1, 3000,'Gd');text(.1, 2e5,'G');
text(1500,0.03,'2500');

% n Buffer tanks
w=logspace(-2,4,121);
tank1=tf(1, [1 1]);
tank1_f = freqresp(tank1,w);
for i=1:length(tank1_f)
    tank1_f2(i)=norm(tank1_f(1,1,i));
end
%Figure 5.24
figure(4);clf;
% vplot('liv,lm',tank1_f);
loglog(w,tank1_f2)
axis([1.e-1,1e2,1.e-3,1]);
hold;

n=4;
for i=2:n 
tank0=tf(1, [1/i 1]); tank=tank0;
for j=2:i tank=tank*tank0;end
tank_f = freqresp(tank,w);
for k=1:length(tank_f)
    tank_f2(k)=norm(tank_f(1,1,k));
end

loglog(w,tank_f2);end
xlabel('Frequency*tauh');ylabel('Magnitude');
text(14,.1,'n=1');text(16,0.02,'n=2');
text(18,0.005,'n=3');text(11,0.002,'n=4');

