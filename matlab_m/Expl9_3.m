%Example 9.3 Glover-McFarlane H_inf loop shaping
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite

% $Id: Expl9_3.m,v 1.2 2004/04/15 08:10:13 vidaral Exp $

clear all; close all;

time=0:0.01:3;
dist=ones(size(time)); 
omega=logspace(-2,3,401);

%Plant (9.75):        
G=tf([0 200],[10 1])*tf([0 1],[0.0025 0.1 1]);
Gd=tf([0 100],[10 1]);
% Weight (9.76)
W1=tf([1 2],[1 0]);
Gs=G*W1;
[a,b,c,d]=ssdata(Gs);

%  Robustify this shaped plant with respect to coprime uncertainty
gammarel=1.1;
[Ac,Bc,Cc,Dc]=coprimeunc(a,b,c,d,gammarel);
[num,den]=ss2tf(Ac,-Bc,Cc,-Dc);
Ks=tf(num,den);
K = W1*Ks; 
Lr =G*K;

L=Gs;
S=inv(1+L);
SGd=S*Gd;
[magL,phaL]=bode(L,omega);
ys=lsim(SGd,dist,time);

l=Lr;
s=inv(1+l);
sgd=s*Gd;
[magl,phal]=bode(l,omega);
yr=lsim(sgd,dist,time);

subplot(121) 
loglog(omega,magL(:),'b--',omega,magl(:),'r',omega,omega./omega,':')
axis([.01,100,.01,10000]); 
xlabel('Frequency')
ylabel('Magnitude')
text(0.3,800,'Gs');
text(0.05,200,'GsKs');
drawnow;

subplot(122)
plot(time,yr,time,ys,'--',time,time*0,':')
axis([0,3,-0.5,1]);
xlabel('Time');
