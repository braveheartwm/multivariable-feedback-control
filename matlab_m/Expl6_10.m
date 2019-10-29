%Example 6.10 Feedback control of distillation process, DV-model
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl6_8.m,v 1.1 2004/01/20 10:47:35 aske Exp $

clear all 
close all

%Model (6.92)
G0=[-87.8 1.4;-108.2 -1.4];
w=logspace(-1,2,61);
dyn = tf(1,[75 1]); 
Dyn = append(dyn,dyn); 
G = Dyn*G0;

%RGA(G)
Lambda=rga(G0)  

%Inverse-based controller
kdyn=tf(0.7,[1 1e-6]);
K=append(kdyn,kdyn)*inv(G0);

w1=tf(0.2);
w1f=frd(w1,w);
w2=tf(0.2*[5 1],[0.5 1]);
w2f=abs(frd(w2,w));
gamma=1.1096;
s=inv(1+tf(0.7, [1 0]));
sf=abs(frd(s,w));
t=1-s;
tf=abs(frd(t,w));

wt1=w1*t;
wt1f=abs(frd(wt1,w));
wt2=w2*t;
wt2f=abs(frd(wt2,w));

L1=sf*(1+wt1f*inv(1+wt1f));
L2=sf*(1+wt2f*inv(1+wt2f));
U11=gamma*sf*inv(1-w1f*tf);
U12=gamma*sf*inv(1-w2f*tf);
U21=sf*inv(1-gamma*w1f*tf);
U22=sf*inv(1-gamma*w2f*tf);

%Figure 6.3
subplot(221);
uplot('liv,lm',L1,U11,U21)
hold on
plot(w,1,'k:')
axis([0.1 100 0.6 2.5]);
ylabel('Magnitude');xlabel('Frequency |rad/min|');
text(2,0.95,'L1');text(2,1.06,'U1');text(2,1.18,'U2');
set(gca,'Ytickmode','manual');
set(gca,'Ytick',[1 2]);
subplot(222);
uplot('liv,lm',L2,U12,U22)
hold on
plot(w,1,'k:');
axis([0.1 100 0.6 2.5]);
ylabel('Magnitude');xlabel('Frequency |rad/min|');
text(2,1.3,'L1');text(1.5,2.2,'U1');text(20,1.25,'U2');
set(gca,'Ytickmode','manual');
set(gca,'Ytick',[1 2]);

