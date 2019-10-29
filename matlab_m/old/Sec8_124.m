%Section 8.12.4 Example: mu-synthesis with DK-iteration
%
%DK-iteration (Table 8.2).
clear all;close all;
%Plant (8.143):
G0 = [87.8 -86.4; 108.2 -109.6];
dyn = tf(1,[75 1]); G = dyn*eye(2)*G0;
%  Weights
wp=.5*tf([10 1],[10 1.e-5]);
Wp = .5*tf([10 1],[10 1.e-5])*eye(2);
Wi = tf([1 0.2],[0.5 1])*eye(2); 
%  Generalized plant P
systemnames = 'G Wp Wi';
inputvar = '[udel(2); w(2) ; u(2)]';
outputvar = '[Wi; Wp; -G-w]';
input_to_G = '[u+udel]';
input_to_Wp = '[G+w]'; input_to_Wi = '[u]';
sysoutname = 'P'; cleanupsysic = 'yes'; sysic;
P=minreal(ss(P));
%  Initialize
omega = logspace(-3,3,61);
% blk = [1 1; 1 1; 2 2; 2 2]; 
blk = [1 1; 1 1; 2 2]; 
nmeas=2; nu=2; d0 = 1; gmin=0.9; gamma=2; tol=0.001;
d1= append(d0,d0,tf(eye(2)),tf(eye(2))); %dsysr=dsysl;
D= append(d0,d0,tf(eye(2)),tf(eye(2)));
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
DPD=d1*P/d1;
gmax=1.05*gamma;
% [K,Nsc,gamma,info] = hinfsyn(DPD,nmeas,nu,'Gmin',gmin,'Gmax',gmax,...
%     'method','lmi','Tolgam',1e-2);
  [K,Nsc,gamma] = hinfsyn(DPD,nmeas,nu,gmin,gmax,tol);
gamma=gamma
Nf = frd(lft(P,K),omega);
size(Nf)


% STEP 2: Compute mu using upper bound:

[mubnds, Info] = mussv(Nf,blk,'c');
murp = norm(mubnds(1,1),inf)
figure(1)
uplot('liv,m',mubnds(1,1),color(i));text(mutextx(i),mutexty(i),mutext(i,:));
 xlabel('Frequency');ylabel('mu');title('MU FOR RP');
  axis([0.001,1000,0.5,1.2]);hold on;drawnow;
% SETP 3: Fit resulting D-scales:
figure(2)
bodemag(mubnds(1,1),omega);
[dsysl,dsysr] = mussvunwrap(Info);
[VDelta,VSigma,VLmi] = mussvextract(Info)
size(dsysl)
dsysl = dsysl/dsysl(3,3);
size(dsysl)

d1=fitfrd(genphase(dsysl(1,1)),4);

d1fit=frd(dsysl(1,1),omega);

  if i<3,
    figure(3); %Figure 8.18
    uplot('liv,lm',dsysl(1,1),color(i),d1fit,[color(i) ':']);
    axis([.001,1000,.01,5]);text(dtextx(i),dtexty(i),dtext(i,:));
    xlabel('Frequency');ylabel('Magnitude');title('D-SCALE (d1)');
    hold on;
    drawnow;
  end
  
%   DPD = mmult(dsysl,P,minv(dsysr)); gmax=1.05*gamma;
%   [K,Nsc,gamma] = hinfsyn(DPD,nmeas,nu,gmin,gmax,tol);
%   gamma=gamma         % 1.1823, 1.0238, 1.0189
%   Nf=frsp(Nsc,omega); % (Remark: Without scaling: N=starp(P,K);)
% %  STEP 2: Compute mu using upper bound:
%   [mubnds,rowd,sens,rowp,rowg] = mu(Nf,blk,'c');
%   murp=pkvnorm(mubnds,inf)  % 1.1818, 1.0238, 1.0189
%   muS = sel(mubnds,':',1);
%   figure(1); %Figure 8.17
%   vplot('liv,m',muS,color(i));text(mutextx(i),mutexty(i),mutext(i,:));
%   xlabel('Freqency');ylabel('mu');title('MU FOR RP');
%   axis([0.001,1000,0.5,1.2]);hold on;drawnow;
% % STEP 3:  Fit resulting D-scales:
%   d1S=d1data(dsysl,rowd); %Save D-scales data for plotting.
%   figure(2);
%   [dsysl,dsysr]=musynflp(dsysl,rowd,sens,blk,nmeas,nu); % use 4th order
%New Version:
%[dsysL,dsysR]=msf(Nf,mubnds,rowd,sens,blk); & \% order: 4, 4, 0.\cr
%              dsysl=daug(dsysL,eye(2)); dsysr=daug(dsysR,eye(2)); & \cr
% or: 
%  [dsysl,dsysr]=muftbtch(dsysl,rowd,sens,blk,nmeas,nu,[.26,.1,4,4]);
%  GOTO STEP 1 (unless satisfied with murp)
%   d1fitS = sel(dsysl,1,1);
%   d1fit=frsp(d1fitS,omega);
%   if i<3,
%     figure(3); %Figure 8.18
%     vplot('liv,lm',d1S,color(i),d1fit,[color(i) ':']);
%     axis([.001,1000,.01,5]);text(dtextx(i),dtexty(i),dtext(i,:));
%     xlabel('Freqency');ylabel('Magnitude');title('D-SCALE (d1)');
%     hold on;drawnow;
%   end
end
K3=K;
%Optimal results (8.144):
j = sqrt(-1);
dz = [-1000; -0.25; -0.054]; 
dgain = 2.0e-3;
dp = [-0.013; -0.67+0.56*j; -0.67-0.56*j];
di = zpk(dz,dp,dgain);
dif=frd(di,omega); %vplot('liv,lm',dif);
D = append(di,di,eye(4));
Pmu=P;
DPD = D*Pmu/D;
[Kopt,Minf,gamma] = hinfsyn(DPD,2,2,0.9,1.1,0.001);  % 0.9737 !!
M = Minf; %M=starp(Pmu,Kopt);
Mf=frd(M,omega);
[mubnds,Info] = mussv(Mf,blk,'Uc');
% [pk,pkomega,vindex]=pkvnorm(mubnds,inf) % 0.9737
murp = norm(mubnds(1,1),inf,1e-6);

% muopt = sel(mubnds,':',1);
d1opt=dif;
[U,sval,V]=svd(G0);
Kdiag=V'*frd(Kopt,omega)*U; %xtract(Kdiag,.001) % DIAGONAL!!

% PLOT OF DK-ITERATION

figure(1); %Figure 8.17
uplot('liv,m',mubnds(1),'--');
text(30,0.95,'Optimal');drawnow

figure(3); %Figure 8.18
uplot('liv,lm',d1opt,'--');
text(200,1.2,'Initial');
text(6,0.035,'Optimal');drawnow

Kinff=frd(K,omega);

%Compute mu fpr NP, RS, RP
Pmuf = frd(Pmu,omega);
Mf=lft(Pmuf,Kinff);
RS = Mf(1:2,1:2);
NP = Mf(3:4,3:4);
blk = [1 1; 1 1; 2 2];

[mubnds,Info] = mussv(Mf,blk,'c');
muRP = mubnds(1); %norm(muRP)
%actual worst-case performance
[maxgain,delworst,info] = wcgain(Mf); 
% delworst=delworst
% skewed_mu=vunpck(vinterp(muup,[1 1;1 Inf],1))

[mubnds,Info] = mussv(RS,[1 1; 1 1],'c');
muRS = mubnds(1);% pkvnorm(muRS)

[mubnds,Info] = mussv(NP,[2 2],'c');
muNP = mubnds(1); %pkvnorm(muNP)
%Figure 8.19
figure(2);clf;
uplot('liv,m',muRP,muRS,muNP);
xlabel('Frequency'); ylabel('mu');
%This is for mu-optimal controller
  title(' MU - PLOTS for controller K3');
  text(.01,1.09,'RP'); text(.01,0.89,'NP'); text(.01,.28,'RS');
  axis([.001,1000,0,1.2]);drawnow

% Now compute the sensitivity for  a few cases

wpfinv = inv(frd(wp,omega));
% Check nominal sensitivity.
Sf = inv(eye(2)+frd(G,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
figure(4)
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
axis([.001,1000,0,1.2]);drawnow
s0=Sfs(1,1);
% OK! Both singular values less than bound.

% Three extreme cases of the gain uncertainty

Gunc1 = G*[1.2 0; 0 1.2];
Gunc2 = G*[0.8 0; 0 1.2];
Gunc3 = G*[1.2 0; 0 0.8];
Gunc4 = G*[0.8 0; 0 0.8];

% Check sensitivity
Sf = inv(eye(2)+frd(Gunc1,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
s1=Sfs(1,1);

Sf = inv(eye(2)+frd(Gunc2,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
s2=Sfs(1,1);

Sf = inv(eye(2)+frd(Gunc3,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
s3=Sfs(1,1);

Sf = inv(eye(2)+frd(Gunc4,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
s4=Sfs(1,1);
% Mu-optimal controller: All cases OK!
%    Plant 4 is the worst at low frequency because of low gain
%    Plants 2 and 3 have the highest peak

% The uncertainty weight has magnitude 0.2 (5s+1)/(0.5s+1)
% One allowed perturbations is l(s) = 0.2 (-5s+1)/(0.5s+1)
%  which yields the input gain k1 = 1+l(s) = 1.2 (-0.417s+1)/(0.5s+1)
% Another one is l(s)= -0.2 (5s+1)/(0.5s+1)
%  which yields the input gain k2 = 1+l(s) = 0.8 (-0.633s+1)/(0.5s+1)
k1 = 1.2*tf([-0.417 1],[0.5 1]);
k2 = 0.8*tf([-0.633 1],[0.5 1]);
Gunc5 =G*append(k1,k1);
Gunc6 = G*append(k2,k1);

% Check sensitivity
Sf = inv(eye(2)+frd(Gunc5,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
s5=Sfs(1,1);
% Mu-opt: OK as expected, but note how it almost reaches the bound 
% at high frequency

Sf = inv(eye(2)+frd(Gunc6,omega)*Kinff);
[udum,Sfs,vdum] = svd(Sf); 
uplot('liv,lm',Sfs(1,1),wpfinv,'--')%[omega' ones(size(omega))'],':');
s6=Sfs(1,1);
% Mu-opt: OK as expected. This one almost touches the bound both 
% at low and high frequency

%Finally, use the worst-case uncertainty computed by wcperf for the
% mu-optimal controller
% delworst = [ -1.8518         0   -1.9250         0    2.0000;
%                    0   -0.1737         0   -0.5896         0;
%               1.9250         0    1.0005         0         0;
%                    0   -0.5896         0   -1.0005         0;
%                    0         0         0         0      -Inf]

A=[-1.8518 0; 0 -0.1737];
B=[-1.9159 0; 0 -0.5896];
C=[1.9250 0; 0 -.5896];
D=[1 0; 0 -1];
delworst=ss(A,B,C,D);

Guncwc = G*(eye(2)+Wi*delworst);
Sf = inv(eye(2)+frd(Guncwc,omega)*Kinff);
[udum,Sfs1,vdum] = svd(Sf); s7 = Sfs1(1,1);
uplot('liv,lm',s7,wpfinv,'--');drawnow
% [Peak_of_S,iv,vindex] = pkvnorm(s7,'inf')
%Figure 8.20
figure(4);clf;
uplot('liv,lm',s0,wpfinv,'--',s1,':',s2,':',s3,':',s4,':',s5,':',s6,':',s7)
axis([0.01,1000,.1,3]);
xlabel('Frequency'), ylabel('sv'); text(5,2.2,'1/WP');

% TIME  simulation
I2=eye(2);
K=K3;
% Nominal
GK = G*K;
S = inv(I2+GK); T = I2-S;
kr=tf(1,[5 1]); Kr=append(kr,kr); Tr = minreal(T*Kr);
t=0:1:100;
r=[ones(length(t),1) zeros(length(t),1)];
y =lsim(Tr,r,t);
u = lsim(K*S*Kr,r,t);
% With 20% uncertainty
Unc = [1.2 0; 0 0.8];
GKu = G*Unc*K;
Su = inv(I2+GKu); Tu = I2-Su;
norm(S,inf), norm(Su,inf)  % Kinv: Peak of S  is 1 (nominally) and 14.2 (with uncertainty)
norm(T,inf), norm(Tu,inf)  % Kinv: Peak of T  is 1 (nominally) and 14.2 (with uncertainty)
Tru = Tu*Kr;
yu = lsim(Tru,r,t);
uu = lsim(K*Su*Kr,r,t);
%Figure 8.21
figure(5);clf;
plot(t,y,t,yu,'--');
axis([0 100 -0.1 1.1]); text(77,0.9,'y1'), text(77,0.1,'y2');
xlabel('Time');

