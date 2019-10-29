% File updated August 2008
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

%--------------------------------------------------------------

display(' 1. Manual mu-synthesis (DK-iteration)')

%% Iteration 1
nmeas=2; nu=2;  
d10 = 1; d20 = 1; D=append(d10,d20,tf(eye(2)),tf(eye(2))); 
P1 = minreal(D*P*inv(D));
[K1,Nsc1,gamma(1),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-4);
  
%% Compute mu using upper bound
omega = logspace(-3,3,61);
Nf=frd(lft(P,K1),omega);
blk = [1 1; 1 1; 2 2];
[mubnds1,Info] = mussv(Nf,blk,'c'); murp(1)=norm(mubnds1(1,1),inf,1e-6) % mu=1.1798

% Iteration 2 
%% Fit resulting D-scales
[dsysl1,dsysr]=mussvunwrap(Info); dsysl1 = dsysl1/dsysl1(3,3);
d11 = fitfrd(genphase(dsysl1(1,1)),4); d21 = fitfrd(genphase(dsysl1(2,2)),4);
%% Design controller
D=append(d11,d21,tf(eye(2)),tf(eye(2))); 
P1 = minreal(D*P*inv(D));
[K2,Nsc2,gamma(2),info]=hinfsyn(P1,nmeas,nu,'Tolgam',1e-3);
%% Compute mu using upper bound
Nf=frd(lft(P,K2),omega);
[mubnds2,Info] = mussv(Nf,blk,'c'); murp(2)=norm(mubnds2(1,1),inf,1e-6) % mu=1.0275

%% Iteration 3
[dsysl2,dsysr]=mussvunwrap(Info); dsysl2 = dsysl2/dsysl2(3,3);
d12 = fitfrd(genphase(dsysl2(1,1)),4); d22 = fitfrd(genphase(dsysl2(2,2)),4);
D=append(d12,d22,tf(eye(2)),tf(eye(2))); 
P1 = minreal(D*P*inv(D));
[K3,Nsc3,gamma(3),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);
Nf=frd(lft(P,K3),omega);
[mubnds3,Info] = mussv(Nf,blk,'c'); murp(3)=norm(mubnds3(1,1),inf,1e-6) % mu=1.0207

%--------------------------------------------------------------
% --- in book we stop after 3 iterations, but let us continue some more...

display('More iterations...')
for i = 4:9
    [dsysl,dsysr]=mussvunwrap(Info); dsysl = dsysl/dsysl(3,3);
    d1 = fitfrd(genphase(dsysl(1,1)),4); d2 = fitfrd(genphase(dsysl(1,1)),4);
    D=append(d1,d2,tf(eye(2)),tf(eye(2))); 
    P1 = minreal(D*P*inv(D));
    [K,Nsc,gamma1(i),info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);
    Nf=frd(lft(P,K),omega);
    [mubnds,Info] = mussv(Nf,blk,'c'); murp(i)=norm(mubnds(1,1),inf,1e-6)
end    
% Get for the first 11 iterations:
%  murp = 1.1798    1.0275    1.0207    1.0130    1.0069    1.0023    0.9998    0.9984    0.9975  0.9976 0.9977
% --- time to stop!  mu is startying to increase  ..


%--------------------------------------------------------------

display('2. Optimal controller from Lundstr�m')
s = tf('s');
dopt = 2*(0.001*s+1)*(s+0.25)*(s+0.054)/((s+0.67)^2 + 0.56^2)/(s+0.013);

D=append(dopt,dopt,tf(eye(2)),tf(eye(2))); 
P1 = minreal(D*P*inv(D));

[Kopt,Nscopt,gammaopt,info]=hinfsyn(P1,nmeas,nu,'method','lmi','Tolgam',1e-3);

Nf=frd(lft(P,Kopt),omega);'method','lmi'
[mubndsopt,Info] = mussv(Nf,blk,'c'); murpopt=norm(mubndsopt(1,1),inf,1e-6) % mu=0.9740


%--------------------------------------------------------------

display('3. Plotting the results') 
bodemag(mubnds1(1,1),mubnds2(1,1),mubnds3(1,1),mubndsopt(1,1)) %Changes in mu

figure,bodemag(d11,d12,dopt,omega) %Changes in D-scale

Snom=minreal(inv(eye(2)+G*K3)); Tnom = minreal(eye(2)-Snom);

figure,sigma(Wp*Snom,Wi*Tnom,frd(mubnds3(1,1),omega),omega) %Mu plots with K3

display('4. Worst case analysis')
Delta = [ultidyn('Delta1',[1 1]) 0;0 ultidyn('Delta2',[1 1])];
Gpert = G*(eye(2) + Wi*Delta);
Spert = inv(eye(2)+Gpert*K3);
WpSpert = Wp*Spert;
Spertf = ufrd(Spert,omega); WpSpertf = ufrd(WpSpert,omega); 
[perfmarg,perfmargunc,report,info] = robustperf(WpSpertf);
WpSwc = usubs(WpSpertf,perfmargunc); Swc = usubs(Spertf,perfmargunc); 
norm(WpSwc,inf), norm(Swc,inf)
'method','lmi'
% For the wosrt case perturbation norm of S is close to 1. This is not
% suprising, as there exists perturbations for which norm of S is close to
% 2, but these perturbations are not worst w.r.t Wp*S.

%Perturbation Analysis

E1 = diag([1.2 1.2]); E2 = diag([0.8 1.2]);
E3 = diag([1.2 0.8]); E4 = diag([0.8 0.8]);
f1 = 1 + tf([-1 0.2],[0.5 1]); f2 = 1 - tf([1 0.2],[0.5 1]);
E5 = [f1 0;0 f1]; E6 = [f2 0;0 f1];

for i = 1:6
    eval(['E=E',int2str(i),';']);
    temp = sigma(minreal(inv(eye(2)+G*E*K3)),omega);
    Ssig = temp(1,:);
    eval(['Ssig',int2str(i),'=Ssig;']);
end

temp = sigma(Snom,omega); Ssignom = temp(1,:);
temp = sigma(Swc,omega); Ssigwc = temp(1,:);
Wpf = abs(freqresp(Wp(1,1),omega));

figure, loglog(omega,Ssig1,':',omega,Ssig2,':',omega,Ssig3,':',omega,Ssig4,':',omega,Ssig5,':',omega,Ssig6,':');
hold
loglog(omega,Ssignom,omega,Ssigwc,'r');
loglog(omega,1./Wpf(:),'--');

display(' 5. Time responses')
S3=minreal(inv(eye(2)+G*E3*K3));T3 = minreal(eye(2)-S3);
R = [tf(1,[5 1]) 0;0 1];
[y,t] = step(T3*R,0:0.1:100);
yp = y(:,:,1);

[y,t] = step(Tnom*R,0:0.1:100);
ynom = y(:,:,1);

figure,plot(t,yp(:,:),'--',t,ynom(:,:))

return

% ------------------------------------------------------------------
display(' 6. Using Automatic software (dksyn)')
Punc = lft([ultidyn('DeltaI_1',[1 1]) 0;0 ultidyn('DeltaI_2',[1 1])],P);

opt=dkitopt('FrequencyVector',logspace(-3,3,61),'DisplayWhileAutoIter','on');

[K,clp,bnd,dkinfo]=dksyn(Punc,nmeas,nu,opt); % 
% Remark August 2008: 
%The mu-values for the first 4 iterations of the automated dk-iteration are:  mu = 1.212, 1.049, 1.050, 1.053
% This is better than what I found in 2005 (1.094 after 4 iterations) 
% ... but still not as good as the manual iterations which give for the first 3 iterations:     mu = 1.180, 1.027, 1.021

% Remark: It is possible to specify an initial controller in the dksyn. But
% with this controller the results are not so good. 





    
