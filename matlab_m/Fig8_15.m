% Figure 8.15 (Table 8.1, mu-analysis)
%    
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig8_15.m,v 1.3 2004/04/14 13:16:22 vidaral Exp $

clear all; close all;
G0=[87.8 -86.4; 108.2 -109.6];
G=tf([1],[75 1])*G0;
G=minreal(ss(G)); % ss

% Inversed based controller
Kinv=0.7*tf([75 1],[1 1e-5])*inv(G0);
Kinv=minreal(ss(Kinv)); % ss

% Weights
wp=0.5*tf([10 1],[10 1e-5]);
wp=minreal(ss(wp)); % ss
wi=tf([1 0.2],[0.5 1]);
wi=minreal(ss(wi)); % ss
Wp=wp*eye(2); % 2x2
Wi=wi*eye(2); % 2x2

% Generalized plant P
systemnames = 'G Wp Wi';
inputvar    = '[ydel(2); w(2) ; u(2)]';
outputvar   = '[Wi ; Wp ; -G-w]';
input_to_G  = '[u+ydel]';
input_to_Wp = '[G+w]';
input_to_Wi = '[u]';
sysoutname  = 'P';
cleanupsysic= 'yes';
sysic;
N=minreal(lft(P,Kinv));

% mu for RP
blk=[1 1; 1 1; 2 2];
omega = logspace(-3,3,61);
Nf=frd(N,omega);
[mubnds,muinfo]=mussv(Nf,blk,'c');
muRP=mubnds(:,1); muRPinf=fnorm(muRP); [muRPinf,muRPw] = norm(muRP,inf) % muRPinf=5.7726

% Worst case weigthed sensitivity
Nsys=ltisys(N.a,N.b,N.c,N.d,1);  % E=ones(dim(A))
delta=ublock(2,1);
[muup,muupfreq]=muperf(Nsys,delta)   % muup=44.92

% mu for RS
Nrs=Nf(1:2,1:2); % Picking out wITi
blk=[1 1; 1 1];
[mubnds,muinfo]=mussv(Nrs,blk,'c');
muRS=mubnds(:,1); muRSinf=fnorm(muRS); [muRSinf,muRSw]=norm(muRS,inf) % muRSinf=0.5242

% mu for NS (=max. singular value of Nrp)
Nnp=Nf(3:4,3:4); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,blk,'c');
muNS=mubnds(:,1); muNSinf=fnorm(muNS); [muNSinf,muNSw]=norm(muNS,inf) % muNSinf=0.500
bodemag(muRP,'',muRS,'--',muNS,'-.',omega)
xlabel('Frequency');ylabel('ssv');
text(0.01,0.5,'RS');
text(20,0.8,'NP');
text(0.3,4.5,'RP');
