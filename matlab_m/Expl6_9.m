%Example 6.9 Feedback control of distillation process.
% Uses the function rga
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl6_7.m,v 1.3 2004/04/14 08:13:34 vidaral Exp $

clear all; close all;
%Model (6.88):
G0 = [87.8 -86.4; 108.2 -109.6];
[U,sval,V]=svd(G0);
tau=75; dyn = tf(1,[tau 1]); 
G = dyn*G0;

% 1. Inverse-based control
dynk = tf(0.7*[tau 1],[1 1.e-6]); 
Kinv = dynk*inv(G0);

% Check of peak value: Inverse based controller (Kinv) Check with RGA-bound
w = logspace(-2,2,81);
k=0.7;  Tinv = tf(k, [1 k]);  Sinv=1-Tinv; 
Tinvf=frd(Tinv,w);  Sinvf=frd(Sinv,w);
rgamaxrow=norm(rga(G0),'inf') % 69.1
wit=abs(0.2*Tinvf); 
witrga = rgamaxrow*wit*inv(1+wit);
Smaxrga =abs(Sinvf)*(1+witrga);

% Peak value 6.81 (at w=0.7943), which is less that 14.2
[Sinv_rga, Sinv_rgaw]=norm(Smaxrga,inf)

% Actual peak value, with uncertainty Ei=diag(0.2, -0.2).
L=minreal(G*[1.2 0;0 0.8]*Kinv);
Sinv2=inv(eye(2)+L);

% Actual peak value 14.2077 (at w=0.6880)
[Sinv_act, Sinv_actw]=norm(Sinv2,inf) % 14.2

% Compare with worst-case performance computed with wcsens
% Uncertainty weight
Gunc=G*(eye(2)+[0.2*ultidyn('DeltaI1',[1 1]) 0; 0 0.2*ultidyn('DeltaI2',[1 1])]) ;
Lunc=Gunc*Kinv;
wcs=wcsens(Lunc,eye(2),'So');

% The worst case exact is 14.3606 at w=0.7908
Sinv_wc=wcs.So.UnweightedMaxGain;


% 2. Diagonal feedback controller.
k2=2.4e-2;
kdyn = tf(k2*[tau 1],[1 1e-6]); 
Kdiag=kdyn*[1 0;0 -1];

% Norminal S and T:
Sdiag=inv(eye(2)+G*Kdiag); Tdiag=eye(2)-Sdiag;
gok=1;wi=0.2;
Ssigma=sigma(frd(Sdiag,w)); Tsigma=sigma(frd(Tdiag,w));
Ssmax=max(max(Ssigma)); Tsmax=max(max(Tsigma));

% Upper bound of uncertain ||S||, using (6.82)
Subd=Ssmax*inv(1-gok*wi*Tsmax);
Sdiag_ubd=norm(Subd,inf) %1.26

% Actual peak value, with uncertainty Ei=diag(0.2, -0.2).
 L=G*[1.2 0; 0 0.8]*Kdiag;
 Sdiag2=inv(eye(2)+L);
 Sdiag_act=norm(Sdiag2,inf)

% Compare with worst-case performance computed with wcsens
Lunc=Gunc*Kdiag;
wcs=wcsens(Lunc,eye(2),'So');

% Worst case exact is  1.0492 at w=1.2753
Sdiag_wc=wcs.So.UnweightedMaxGain;

% Result:
clc
disp('1. Inverse-based controller:');
disp(sprintf('   Lower bound of ||S|| estimated using RGA:%6.2f',Sinv_rga));
disp(sprintf('   Actual ||S|| with E_I=diag(0.2,-0.2):    %6.2f',Sinv_act(1)));
disp(sprintf('   Worst-case performance using wcsens:     %6.2f',Sinv_wc));
disp('2. Diagonal controller:');
disp(sprintf('   Upper bound of ||S|| using (6.77):       %6.2f',Sdiag_ubd));
disp(sprintf('   Actual ||S|| with E_I=diag(0.2,-0.2):    %6.2f',Sdiag_act(1)));
disp(sprintf('   Worst-case performance using wcsens:     %6.2f',Sdiag_wc));