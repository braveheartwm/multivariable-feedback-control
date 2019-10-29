clear all
close all
clc;
set( cstprefs.tbxprefs, 'MagnitudeUnits', 'abs', 'MagnitudeScale', 'log' ) ;


% The shaped plant of the model and the unscaled controller 
% are saved in aeroK.mat from the design, in Sec12_33.m.
% Also the reference model is inside this file.

load aeroK

% control plant
[As,Bs,Cs,Ds] = ssdata(balreal(Gs));

% reference model 
[Ar,Br,Cr,Dr] = ssdata(Mo);

% reference model states 
[nr,nr] = size(Ar);
[lr,mr] = size(Dr);

[ns,ns] = size(As);
[ls,ms] = size(Ds);


Rs = eye(ls) + Ds*Ds'
Ss = eye(ms) + Ds'*Ds

A = (As - Bs * inv(Ss) * Ds'* Cs)
B = Cs'*sqrt(inv(Rs));
Q = Bs * inv(Ss)*Bs';
[Zs,ZAMP,G,REP] = care(A,B,Q);


rho = 2;
A = blkdiag(As,Ar);
B1 = [zeros(ns,mr) ((Bs*Ds'+Zs*Cs'))*inv(sqrt(Rs));
    Br zeros(nr,ls)];
B2 = [Bs;zeros(nr,ms)];

C1 = [zeros(ms,ns+nr);Cs zeros(ls,nr);rho*Cs -rho*rho*Cr];
C2 = [zeros(mr,ns+nr);Cs zeros(ls,nr)];

D11 = [zeros(ms,mr+ls);zeros(ls,mr) sqrt(Rs); -rho*rho*Dr rho*sqrt(Rs)];
D12 = [eye(ms);Ds;rho*Ds];
D21 = [rho*eye(mr) zeros(mr,ls);zeros(ls,mr) sqrt(Rs)];
D22 = [zeros(mr,ms);Ds];

B = [B1 B2];
C = [C1;C2];
D = [D11 D12;D21 D22];
P = ss(A,B,C,D);

[l1,m2] = size(D12);
[l2,m1] = size(D21);
nmeas = 6;ncon=m2;
gmin = 1;
gmax = 5;
gtol = 0.01;



% S = eye(size(D'*D))+D'*D;
% R = eye(size(D*D'))+D*D';
% Rinv = inv(R);
% Sinv = inv(S);
% A1 = (A-B*Sinv*D'*C);
% R1 = S;
% B1 = B;
% Q1 = C'*Rinv*C;
% [X,XAMP,G] = care(A1,B1,Q1,R1);
% A2 = A1';
% Q2 = B*Sinv*B';
% B2 = C';
% R2 = R;
% [Z,ZAMP,G] = care(A2,B2,Q2,R2);
% 
% XZ = X*Z;
% gammin = sqrt(1+max(eig(XZ)));
% 
% gam = 1.1*gammin;
% gam2 = gam*gam;
% gamconst = (1-gam2)*eye(size(XZ));
% Lc = gamconst +XZ;
% Li = inv(Lc');
% Fc = -Sinv*(D'*C+B'*X);
% Ac =A + B*Fc + gam2*Li*Z*C'*(C+D*Fc);
% Bc = gam2*Li*Z*C';
% Cc = B'*X;
% Dc = -D';
% K = ss(Ac,Bc,Cc,Dc)

[K,Gnclp,gam] = hinfsyn(P,nmeas,ncon,gmin,gmax,gtol);


Kf1=K(1:3,1:3);
Kf2=K(1:3,4:6);

systemnames  = 'Gs Kf1 Kf2';
inputvar     = '[ r( 3 ) ]';
outputvar    = '[ Gs ]';
input_to_Gs   = '[ Kf1 + Kf2 ]';
input_to_Kf1 = '[ r ]';
input_to_Kf2 = '[ Gs ]';
sysoutname   = 'Ms';
cleanupsysic = 'yes';
sysic;

Mo0 = freqresp( Mo, 0 ) ;
M0  = freqresp( Ms, 0 ) ;
Kf1 = Kf1*inv(M0)*Mo0;

systemnames  = 'Gs Kf1 Kf2';
inputvar     = '[ r( 3 ) ]';
outputvar    = '[ Gs ]';
input_to_Gs   = '[ Kf1 + Kf2 ]';
input_to_Kf1 = '[ r ]';
input_to_Kf2 = '[ Gs ]';
sysoutname   = 'Ms';
cleanupsysic = 'yes';
sysic;

t=0:0.01:2;

y1=lsim(Ms,ones(size(t'))*[1 0 0],t);
y2=lsim(Ms,ones(size(t'))*[0 1 0],t);
y3=lsim(Ms,ones(size(t'))*[0 0 1],t);

y10=lsim(Mo,ones(size(t'))*[1 0 0],t);
y20=lsim(Mo,ones(size(t'))*[0 1 0],t);
y30=lsim(Mo,ones(size(t'))*[0 0 1],t);

% Figure 13.13

figure(1); subplot(2,3,1);
plot(t,y1(:,1),'-',t,y1(:,2),'--',t,y1(:,3),'-.',t,1,':',t,y10(:,1),t,y10(:,2),t,y10(:,3));
axis([0 2 -0.1 1.1 ]);
xlabel('Time');

subplot(2,3,2);
plot(t,y2(:,1),'-',t,y2(:,2),'--',t,y2(:,3),'-.',t,1,':',t,y20(:,1),t,y20(:,2),t,y20(:,3));
axis([0 2 -0.1 1.1]);
xlabel('Time');

subplot(2,3,3);
plot(t,y3(:,1),'-',t,y3(:,2),'--',t,y3(:,3),'-.',t,1,':',t,y30(:,1),t,y30(:,2),t,y30(:,3));
axis([0 2 -0.1 1.1]);
xlabel('Time');





