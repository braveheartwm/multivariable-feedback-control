clear all;close all;

set(cstprefs.tbxprefs,'MagnitudeUnits','abs','MagnitudeScale','log');

% System description
G0 = [87.8 -86.4; 108.2 -109.6];
dyn = tf(1,[75 1]); G=dyn*eye(2)*G0;

% Weights
Wp = 0.5*tf([10 1],[10 1.e-5])*eye(2);
Wi = tf([1 0.2],[0.5 1])*eye(2);

% Generalized Plant P
systemnames = 'G Wp Wi';
inputvar = '[di(2); w(2) ; u(2)]'; outputvar = '[Wi; Wp; -G-w]';
input_to_G = '[u+di]'; input_to_Wp = '[G+w]'; input_to_Wi = '[u]';
sysoutname = 'P'; cleanupsysic = 'yes'; sysic;
P = minreal(ss(P)); 

%% Iteration 1
nmeas=2; nu=2;  
d10 = 1; d20 = 1; D=append(d10,d20,tf(eye(2)),tf(eye(2)))
P1 = minreal(D*P*inv(D));
[K1,Nsc1,gamma(1),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);
%[K1,Nsc1,gamma(1),info]=hinfsyn(P1,nmeas,nu,'Tolgam',1e-2);
% [K1,Nsc1,gamma(1),info] = mixsyn(G,Wp,Wi,[]);  
%% Compute mu using upper bound
omega = logspace(-3,3,61);
Nf=frd(lft(P,K1),omega);
blk = [1 1; 1 1; 2 2];
[mubnds,Info] = mussv(Nf,blk,'c'); murp(1)=norm(mubnds(1,1),inf,1e-6)
bodemag(mubnds(1,1),omega); 

%Fit resulting D-scales
[dsysl1,dsysr]=mussvunwrap(Info); dsysl1 = dsysl1/dsysl1(3,3);
d11 = fitfrd(genphase(dsysl1(1,1)),4); d21 = fitfrd(genphase(dsysl1(2,2)),4);

Dtemp = append(d11,d21,tf(eye(2)));
Act_norm(1) = norm(Dtemp*lft(P,K1)*inv(Dtemp),inf);

%% Iteration 2
D=append(d11,d21,tf(eye(2)),tf(eye(2))); 
P1 = minreal(D*P*inv(D));
[K2,Nsc2,gamma(2),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);
%[K2,Nsc2,gamma(2),info]=hinfsyn(P1,nmeas,nu,'Tolgam',1e-2);
  
%% Compute mu using upper bound
Nf=frd(lft(P,K2),omega);
[mubnds,Info] = mussv(Nf,blk,'c'); murp(2)=norm(mubnds(1,1),inf,1e-6);
hold on;
bodemag(mubnds(1,1),omega); 

%Fit resulting D-scales
[dsysl2,dsysr]=mussvunwrap(Info); dsysl2 = dsysl2/dsysl2(3,3);
d12 = fitfrd(genphase(dsysl2(1,1)),4); d22 = fitfrd(genphase(dsysl2(2,2)),4);

Dtemp = append(d12,d22,tf(eye(2)));
Act_norm(2) = norm(frd(Dtemp*lft(P,K2)*inv(Dtemp),omega),inf);

%return

%% Iteration 3
D=append(d12,d22,tf(eye(2)),tf(eye(2))); 
P1 = minreal(D*P*inv(D));

[K3,Nsc3,gamma(3),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);
  
%% Compute mu using upper bound
Nf=frd(lft(P,K3),omega);
[mubnds,Info] = mussv(Nf,blk,'c'); murp(3)=norm(mubnds(1,1),inf,1e-6);
bodemag(mubnds(1,1),omega); 

%Fit resulting D-scales
[dsysl3,dsysr]=mussvunwrap(Info); dsysl3 = dsysl3/dsysl3(3,3);
d13 = fitfrd(genphase(dsysl3(1,1)),4); d23 = fitfrd(genphase(dsysl3(1,1)),4);

Dtemp = append(d13,d23,tf(eye(2)));
Act_norm(3) = norm(frd(Dtemp*lft(P,K3)*inv(Dtemp),omega),inf);

% More iterations

d1 = d13; d2 = d23;

for i = 1:5
    D=append(d1,d2,tf(eye(2)),tf(eye(2))); 
    P1 = minreal(D*P*inv(D));

    [K,Nsc,gamma1(i),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);
  
    %% Compute mu using upper bound
    Nf=frd(lft(P,K),omega);
    [mubnds,Info] = mussv(Nf,blk,'c'); murp1(i)=norm(mubnds(1,1),inf,1e-6);
    
    %Fit resulting D-scales
    [dsysl,dsysr]=mussvunwrap(Info); dsysl = dsysl/dsysl(3,3);
    d1 = fitfrd(genphase(dsysl(1,1)),4); d2 = fitfrd(genphase(dsysl(1,1)),4);

end

    