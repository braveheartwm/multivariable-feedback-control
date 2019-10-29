% Section 13.3.3          Aero-engine control, Part II
%
% Design of multivariable controller for aero engine, 
% reviewed by Kjetil Havre 13/4-1995.
%
% Dependencies: data file aero1.mat
%               function align.m
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
%
% $Id: Sec13_33.m,v 1.2 2004/04/20 13:35:13 vidaral Exp $

clear all;close all;
set( cstprefs.tbxprefs, 'MagnitudeUnits', 'abs', 'MagnitudeScale', 'log' ) ;

% Load the engine model

load aero1
G=G5; 
% Singular values of plant

%Figure 13.12
figure(1); subplot(2,2,1);
sigma(G);
axis([1e-1 1e3 1e-4 1e2]);


% Build weight Wpbar=I/s (where I is a 3*3 identity matrix)

sys1=tf(1, [1 0]);
Wpbar(1,1)=sys1;
Wpbar(2,2)=sys1;
Wpbar(3,3)=sys1;

% Weight crossover equal to 1 rad/s
% Build weight W2=I

W2=eye(3);

% Build W2*G*Wpbar and find the align gain Ka (at 7 rad/sec)

W2GWpbar=W2*G*Wpbar;
[a,b,c,d]=ssdata(W2GWpbar);
ff=c*inv(j*7*eye(size(a))-a)*b+d;	
Ka=align(ff);

% Find condition number of Ka.
% Compare with condition number of G(j7).
fprintf('\nThe condition number of Ka is %0.3f\n', cond(Ka));

% Condition number of G(j7) = 3.77
fprintf('\nThe condition number of G(j7) is %0.3f\n', cond(ff));

% Build Kg, W1 and the final shaped plant W2*G*W1

Kg=diag([1,2.5,0.3]);
W1=Wpbar*Ka*Kg;
W2GW1=balreal(W2*G*W1);
Gs=W2GW1;


% Plot singular values of the shaped plant

subplot(2,2,2);
sigma(Gs);
axis([1e-1 1e3 1e-4 1e2]);
[a,b,c,d]=ssdata(W2GW1);

% Find gamma_0 by solving the two Riccati equations
S=eye(size(d'*d))+d'*d;
Sinv=S\eye(3);
R=eye(size(d*d'))+d*d';
Rinv=R\eye(3);

% Using care
A1=a-b*Sinv*d'*c;
B1=b;
Q1=c'*Rinv*c;
R1=S;

[X,Xleig,Xgain,XWP]=care(A1,B1,Q1,R1);
if (XWP==-1)
 disp( '*** WARNING! X-Riccati Equation may be ill-posed. ***' )
end

% Using care
A2=A1;
B2=c';
Q2=b*Sinv*b';
R2=R;
[Z,Zleig,Zgain,ZWP]=care(A2,B2,Q2,R2);
% Warnings displayed
if (ZWP==-1)
 disp( '*** WARNING! Z-Riccati Equation may be ill-posed. ***' )
end

gamma_0=real(sqrt(1+max(eig(X*Z))))

%  Initial gamma is: 2.3297

rootR=sqrtm(R);
H=-(b*d'+Z*c')*Rinv;

% Build the reference model Mo

sys1=tf(1,[0.018 1]);
sys2=tf(1,[0.008 1]);
sys3=tf(1,[0.2 1]);
Mo(1,1)=sys1;
Mo(2,2)=sys2;
Mo(3,3)=sys3;

Mo=ss(Mo);

% Set rho to 1 and build the generalized plant

rho=1;
A=[a,zeros(size(a,1),size(Mo.a,2));zeros(size(Mo.a, 1),size(a,2)), Mo.a];
invrootR=rootR\eye(3);

B1_12=(b*d'+Z*c')*invrootR;
[m1,n1]=size(B1_12);
[mBo,nBo]=size(Mo.b);
[mb,nb]=size(b);
B1=[zeros(m1,nBo) B1_12;
       -Mo.b zeros(mBo,n1)];
B2=[b; zeros(mBo,nb)];
[mc,nc]=size(c);
[mCo,nCo]=size(Mo.c);
[mDo,nDo]=size(Mo.d);
[mR,nR]=size(rootR);
[md,nd]=size(d);

C1=[zeros(nb,nc) zeros(nb,nCo);
    c zeros(mc,nCo);
    rho*c     rho^2*Mo.c];
C2=[zeros(mc,nc) zeros(mc,nCo);
      c zeros(mc,nCo)];
                 
D11=[zeros(nb,nDo) zeros(nb,nR);
       zeros(mR,nDo)     rootR;
         -( rho^2 )*Mo.d   rho*rootR];
D12=[eye(nb,nd);
         d;
      rho*d];
D21=[rho*eye(mc,nDo) zeros(mc,nR);
            zeros(mR, nDo) rootR];
D22=[zeros(mc,nd);d];

% Do gamma iterations, and calculate the H-inf controller
% given by. gam is the achieved suboptimal gamma.
% Use hinfsyn, since hinfopt is removed
B=[B1 B2]; C=[C1;C2]; D=[D11 D12;D21 D22];

P=ss(A,B,C,D);
[l1,m2]=size(D12); [l2,m1]=size(D21);
nmeas=l2; ncon=m2; gmin=0.2; gmax=5; gtol=0.01;
[K,Gnlcp,gam]=hinfsyn(P,nmeas,ncon,'gmin',gmin,'gmax',gmax,'tolgam',gtol,'display','on');

fprintf('\nThe final gamma value of the H_inf iteration is %0.4f\n', 1 / gam);

% scale the prefilter of the 2-DOF controller to give perfect steady-state
% tracking

% 
%  Form the unscaled controller without integral action.
%  to be used in the model reduction.
%

disp( 'Controller is formed without scaling of prefilter.' )

[acp,bcp,ccp,dcp]=ssdata(K);


dcK1=dcgain(acp,bcp(:,1:3),ccp,dcp(:,1:3));
dcK1inv=dcK1\eye(3);
dcK2=dcgain(acp, bcp(:,4:6),ccp,dcp(:,4:6));
bcp(:,1:3)=bcp(:,1:3)*dcK1inv*(-dcK2);
dcp(:,1:3)=dcp(:,1:3)*dcK1inv*(-dcK2);
K=ss(acp,bcp,ccp,dcp);


% Build the final controller W1*[K1 K2]*[I 0;0 W2]

Kf=W1*K;		% W2=I so no need for post-multiplication by [I 0;0 W2]

Kf1=Kf(1:3,1:3);
Kf2=Kf(1:3,4:6);

systemnames  = 'G Kf1 Kf2';
inputvar     = '[ r( 3 ) ]';
outputvar    = '[ G ]';
input_to_G   = '[ Kf1 + Kf2 ]';
input_to_Kf1 = '[ r ]';
input_to_Kf2 = '[ G ]';
sysoutname   = 'Ms';
cleanupsysic = 'yes';
sysic;
t=0:0.01:2;

y1=lsim(Ms,ones(size(t'))*[1 0 0],t);
y2=lsim(Ms,ones(size(t'))*[0 1 0],t);
y3=lsim(Ms,ones(size(t'))*[0 0 1],t);


% Figure 13.13

figure(2); subplot(2,3,1);
plot(t,y1(:,1),'-',t,y1(:,2),'--',t,y1(:,3),'-.',t,1,':');
axis([0 2 -0.1 1.1 ]);
xlabel('Time');

subplot(2,3,2);
plot(t,y2(:,1),'-',t,y2(:,2),'--',t,y2(:,3),'-.',t,1,':');
axis([0 2 -0.1 1.1]);
xlabel('Time');

subplot(2,3,3);
plot(t,y3(:,1),'-',t,y3(:,2),'--',t,y3(:,3),'-.',t,1,':');
axis([0 2 -0.1 1.1]);
xlabel('Time');

% Plot the sensitivity function
L =G*Kf2;
LI=Kf2*G ;
Se=inv(eye(3)-L);
SI=inv(eye(3)-LI);
T =eye(3)-Se;
TI=eye(3)-SI;

% Peak in S.

[HinfS HinfSw]=norm(Se,inf,1e-6) ;
fprintf( '\nThe peak value for S is %0.2f and occurs at %0.2f rad/s.\n', HinfS, HinfSw ) ;

% Peak in S is: 1.4396 for \w = 2.7499e+01
%
% Figure 13.14

figure(3); subplot(2,2,1);
sigma(Se);
axis([ 1e-2 1e3 1e-3 3e0]);
xlabel('Frequency');

% Output multiplicative robustness

subplot(2,2,2);
sigma(T(1,1),TI(1,1));
axis([1e-2 1e3 1e-3 3e0]);
xlabel('Frequency');
