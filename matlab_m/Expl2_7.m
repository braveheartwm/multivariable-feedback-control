%Example 2.7
%Step response and Frequency response for
%Plant: T=(-s+z)/(s+z)(tau*s+1), z=0.1, tau=1, and S=1-T
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_5.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

%Figure 2.14 Step response 
clear
figure(1);clf;
tau=1;z=0.1; w=logspace(-2,1,101);
s=tf('s');
L=(-s+z)/s/(tau*s+tau*z+2);
T=(-s+0.1)/(s+0.1)/(s+1);
S=1-T;

%step respons T
step(T,50);hold;plot([0 50],[0 0],'g:');
ylabel('y');

%Work out margins, peak values and frequencies.
[Gm,Pm,W180,Wc]=margin(L);
mT=norm(T,inf,1e-4); mS=norm(S,inf,1e-4); 
omega=logspace(-2,1,301);
[x,x,x,Wb]=margin(S/0.707);
warning off
[x,x,x,WbT]=margin(T/0.707);
warning on

disp(sprintf('[Gm,Pm,Ms,Mt]=[%4.1f,%5.1f,%5.2f,%4.1f ]',Gm,Pm,mS(1),mT(1)));
disp(sprintf('[W180,Wc,Wb,Wbt]=[%5.2f,%6.3f,%6.3f,%4.1f ]',W180,Wc,Wb,WbT));

%Figure 2.15 Frequency response
figure(2);clf;
[mag1,pha1]=bode(T,w);
[mag2,pha2]=bode(S,w);
loglog(w,mag1(:),'g',w,mag2(:),'b',w,w./w,'r:',[0.01 WbT],...
       [0.707 0.707],'r:',[Wb Wb],[0.707 1],'r:',[WbT WbT],[0.707 1],'r:')
axis([.01,10,.1,4]);


