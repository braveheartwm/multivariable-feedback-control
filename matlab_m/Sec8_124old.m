%Section 8.12.4 Example: mu-synthesis with DK-iteration
%
%DK-iteration (Table 8.2).
clear all;close all;
%Plant (8.143):
G0 = [87.8 -86.4; 108.2 -109.6];
dyn = nd2sys(1,[75 1]); Dyn = daug(dyn,dyn); G = mmult(Dyn,G0);
%  Weights
wp = nd2sys([10 1],[10 1.e-5],0.5);
Wp = daug(wp,wp); % Approximated integrator
wi = nd2sys([1 0.2],[0.5 1]); Wi = daug(wi,wi);
%  Generalized plant P
systemnames = 'G Wp Wi';
inputvar = '[di(2); w(2) ; u(2)]'; outputvar = '[Wi; Wp; -G-w]';
input_to_G = '[u+di]'; input_to_Wp = '[G+w]'; input_to_Wi = '[u]';
sysoutname = 'P'; cleanupsysic = 'yes'; sysic;
%  Initialize
omega = logspace(-3,3,61);
blk = [1 1; 1 1; 2 2]; nmeas=2; nu=2; gmin=0.9; gamma=2; tol=0.01;
d0 = 1; dsysl = daug(d0,d0,eye(2),eye(2)); dsysr=dsysl;
%  START ITERATION
color=['y';'g';'b'];
mutextx=[0.01;0.05;0.05];
mutexty=[1.15;1.05;1.0];
mutext=['Iter. 1';'Iter. 2';'Iter. 3'];
dtextx=[200;200];
dtexty=[0.11;0.035];
dtext=['Iter. 1';'Iter. 2'];
for i=1:3,

%  STEP 1: Find H-infinty optimal controller with given scalings:
  DPD = mmult(dsysl,P,minv(dsysr)); gmax=1.05*gamma;
  [K,Nsc,gamma] = hinfsyn(DPD,nmeas,nu,gmin,gmax,tol);
  gamma=gamma         % 1.1823, 1.0238, 1.0189
  Nf=frsp(Nsc,omega); % (Remark: Without scaling: N=starp(P,K);)
%  STEP 2: Compute mu using upper bound:
  [mubnds,rowd,sens,rowp,rowg] = mu(Nf,blk,'c');
  murp=pkvnorm(mubnds,inf)  % 1.1818, 1.0238, 1.0189
  muS = sel(mubnds,':',1);
  figure(1); %Figure 8.17
  vplot('liv,m',muS,color(i));text(mutextx(i),mutexty(i),mutext(i,:));
  xlabel('Freqency');ylabel('mu');title('MU FOR RP');
  axis([0.001,1000,0.5,1.2]);hold on;drawnow;
% STEP 3:  Fit resulting D-scales:
  d1S=d1data(dsysl,rowd); %Save D-scales data for plotting.
  figure(2);
  [dsysl,dsysr]=musynflp(dsysl,rowd,sens,blk,nmeas,nu); % use 4th order
%New Version:
%[dsysL,dsysR]=msf(Nf,mubnds,rowd,sens,blk); & \% order: 4, 4, 0.\cr
%              dsysl=daug(dsysL,eye(2)); dsysr=daug(dsysR,eye(2)); & \cr
% or: 
%  [dsysl,dsysr]=muftbtch(dsysl,rowd,sens,blk,nmeas,nu,[.26,.1,4,4]);
%  GOTO STEP 1 (unless satisfied with murp)
  d1fitS = sel(dsysl,1,1);
  d1fit=frsp(d1fitS,omega);
  if i<3,
    figure(3); %Figure 8.18
    vplot('liv,lm',d1S,color(i),d1fit,[color(i) ':']);
    axis([.001,1000,.01,5]);text(dtextx(i),dtexty(i),dtext(i,:));
    xlabel('Freqency');ylabel('Magnitude');title('D-SCALE (d1)');
    hold on;drawnow;
  end
end
K3=K;
%Optimial results (8.144):
j = sqrt(-1);
dz = [-1000; -0.25; -0.054]; 
dgain = 2.0e-3;
dp = [-0.013; -0.67+0.56*j; -0.67-0.56*j];
di = zp2sys(dz,dp,dgain);
dif=frsp(di,omega); %vplot('liv,lm',dif);
D = daug(di,di,eye(4));
Pmu=P;
DPD = mmult(D,Pmu,minv(D));
[Kopt,Minf,gamma] = hinfsyn(DPD,2,2,0.9,1.1,0.001);  % 0.9737 !!
M = Minf; %M=starp(Pmu,Kopt);
Mf=frsp(M,omega);
[mubnds,rowd,sens,rowp,rowg] = mu(Mf,blk,'Uc');
[pk,pkomega,vindex]=pkvnorm(mubnds,inf) % 0.9737
muopt = sel(mubnds,':',1);
d1opt=dif;
[U,sval,V]=svd(G0);
Kdiag=mmult(V',frsp(Kopt,omega),U); xtract(Kdiag,.001) % DIAGONAL!!

% PLOT OF DK-ITERATION

figure(1); %Figure 8.17
vplot('liv,m',muopt,'--');
text(30,0.95,'Optimal');drawnow

figure(3); %Figure 8.18
vplot('liv,lm',1,d1opt,'--');
text(200,1.2,'Initial');
text(6,0.035,'Optimal');drawnow

Kinff=frsp(K,omega);

%Compute mu fpr NP, RS, RP
Pmuf = frsp(Pmu,omega);
Mf=starp(Pmuf,Kinff);
RS = sel(Mf,1:2,1:2);
NP = sel(Mf,3:4,3:4);
blk = [1 1; 1 1; 2 2];

[mubnds,rowd,sens,rowp,rowg] = mu(Mf,blk,'c');
muRP = sel(mubnds,':',1); pkvnorm(muRP)
%actual worst-case performance
[delworst,mulow,muup] = wcperf(Mf,blk,1); 
delworst=delworst
skewed_mu=vunpck(vinterp(muup,[1 1;1 Inf],1))

[mubnds,rowd,sens,rowp,rowg] = mu(RS,[1 1; 1 1],'c');
muRS = sel(mubnds,':',1); pkvnorm(muRS)

[mubnds,rowd,sens,rowp,rowg] = mu(NP,[2 2],'c');
muNP = sel(mubnds,':',1); pkvnorm(muNP)
%Figure 8.19
figure(2);clf;
vplot('liv,m',muRP,muRS,muNP);
xlabel('Freqency'); ylabel('mu');
%This is for mu-optimal controller
  title(' MU - PLOTS for controller K3');
  text(.01,1.09,'RP'); text(.01,0.89,'NP'); text(.01,.28,'RS');
  axis([.001,1000,0,1.2]);drawnow

% Now compute the sensitivity for  a few cases

wpfinv = minv(frsp(wp,omega));
% Check nominal sensitivity.
Sf = minv(madd(eye(2),mmult(frsp(G,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
figure(4)
vplot('liv,lm',sel(Sfs,1,1),wpfinv,'--',1,':');drawnow
s0=sel(Sfs,1,1);
% OK! Both singular values less than bound.

% Three extreme cases of the gain uncertainty

Gunc1 = mmult(G,[1.2 0; 0 1.2]);
Gunc2 = mmult(G,[0.8 0; 0 1.2]);
Gunc3 = mmult(G,[1.2 0; 0 0.8]);
Gunc4 = mmult(G,[0.8 0; 0 0.8]);

% Check sensitivity
Sf = minv(madd(eye(2),mmult(frsp(Gunc1,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
vplot('liv,lm',sel(Sfs,1,1),wpfinv,1);drawnow
s1=sel(Sfs,1,1);

Sf = minv(madd(eye(2),mmult(frsp(Gunc2,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
vplot('liv,lm',sel(Sfs,1,1),wpfinv,'--',1,':');drawnow
s2=sel(Sfs,1,1);

Sf = minv(madd(eye(2),mmult(frsp(Gunc3,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
vplot('liv,lm',sel(Sfs,1,1),wpfinv,'--',1,':');drawnow
s3=sel(Sfs,1,1);

Sf = minv(madd(eye(2),mmult(frsp(Gunc4,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
vplot('liv,lm',sel(Sfs,1,1),wpfinv,'--',1,':');drawnow
s4=sel(Sfs,1,1);
% Mu-optimal controller: All cases OK!
%    Plant 4 is the worst at low frequency because of low gain
%    Plants 2 and 3 have the highest peak

% The uncertainty weight has magnitude 0.2 (5s+1)/(0.5s+1)
% One allowed perturbations is l(s) = 0.2 (-5s+1)/(0.5s+1)
%  which yields the input gain k1 = 1+l(s) = 1.2 (-0.417s+1)/(0.5s+1)
% Another one is l(s)= -0.2 (5s+1)/(0.5s+1)
%  which yields the input gain k2 = 1+l(s) = 0.8 (-0.633s+1)/(0.5s+1)
k1 = nd2sys([-0.417 1],[0.5 1],1.2);
k2 = nd2sys([-0.633 1],[0.5 1],0.8);
Gunc5 = mmult(G,daug(k1,k1));
Gunc6 = mmult(G,daug(k2,k1));

% Check sensitivity
Sf = minv(madd(eye(2),mmult(frsp(Gunc5,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
vplot('liv,lm',sel(Sfs,1,1),wpfinv,'--',1,':');drawnow
s5=sel(Sfs,1,1);
% Mu-opt: OK as expected, but note how it almost reaches the bound 
% at high frequency

Sf = minv(madd(eye(2),mmult(frsp(Gunc6,omega),Kinff)));
[udum,Sfs,vdum] = vsvd(Sf); 
vplot('liv,lm',sel(Sfs,1,1),wpfinv,'--',1,':');drawnow
s6=sel(Sfs,1,1);
% Mu-opt: OK as expected. This one almost touches the bound both 
% at low and high frequency

%Finally, use the worst-case uncertainty computed by wcperf for the
% mu-optimal controller
delworst = [ -1.8518         0   -1.9250         0    2.0000;
                   0   -0.1737         0   -0.5896         0;
              1.9250         0    1.0005         0         0;
                   0   -0.5896         0   -1.0005         0;
                   0         0         0         0      -Inf]
Guncwc = mmult(G, madd(eye(2),mmult(Wi,delworst)));
Sf = minv(madd(eye(2),mmult(frsp(Guncwc,omega),Kinff)));
[udum,Sfs1,vdum] = vsvd(Sf); s7 = sel(Sfs1,1,1);
vplot('liv,lm',s7,wpfinv,'--',1,':');drawnow
[Peak_of_S,iv,vindex] = pkvnorm(s7,'inf')
%Figure 8.20
figure(4);clf;
vplot('liv,lm',s0,wpfinv,'--',s1,':',s2,':',s3,':',s4,':',s5,':',s6,':',s7)
axis([0.01,10,.1,3]);
xlabel('Frequency'), ylabel('sv'); text(5,2.2,'1/WP');

% TIME  simulation
I2=eye(2);
K=K3;
% Nominal
GK = mmult(G,K);
S = minv(madd(I2,GK)); T = msub(I2,S);
kr=nd2sys(1,[5 1]); Kr=daug(kr,kr); Tr = mmult(T,Kr);
y = trsp(Tr,[1;0],100,.1);
u = trsp(mmult(K,S,Kr),[1;0],100,.1);
% With 20% uncertainty
Unc = [1.2 0; 0 0.8];
GKu = mmult(G,Unc,K);
Su = minv(madd(I2,GKu)); Tu = msub(I2,Su);
hinfnorm(S), hinfnorm(Su)  % Kinv: Peak of S  is 1 (nominally) and 14.2 (with uncertainty)
hinfnorm(T), hinfnorm(Tu)  % Kinv: Peak of T  is 1 (nominally) and 14.2 (with uncertainty)
Tru = mmult(Tu,Kr);
yu = trsp(Tru,[1;0],100,.1);
uu = trsp(mmult(K,Su,Kr),[1;0],100,.1);
%Figure 8.21
figure(5);clf;
vplot(y,yu,'--');
axis([0 100 -0.1 1.1]); text(77,0.9,'y1'), text(77,0.1,'y2');
xlabel('Time');

