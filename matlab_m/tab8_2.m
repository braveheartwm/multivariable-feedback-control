clear all;close all;

set(cstprefs.tbxprefs,'MagnitudeUnits','abs','MagnitudeScale','log');

% System description
G0 = [87.8 -86.4; 108.2 -109.6];
dyn = tf(1,[75 1]); G=dyn*eye(2)*G0;
nmeas=2; nu=2;  
% Weights
Wp = 0.5*tf([10 1],[10 1.e-5])*eye(2);
Wi = tf([1 0.2],[0.5 1])*eye(2);
omega = logspace(-3,3,61);
% Generalized Plant P
systemnames = 'G Wp Wi';
inputvar = '[di(2); w(2) ; u(2)]'; outputvar = '[Wi; Wp; -G-w]';
input_to_G = '[u+di]'; input_to_Wp = '[G+w]'; input_to_Wi = '[u]';
sysoutname = 'P'; cleanupsysic = 'yes'; sysic;
P = minreal(ss(P)); 

Delta = [ultidyn('di_1',[1 1]) 0;0 ultidyn('di_2',[1 1])];
Punc = lft(Delta,P);
opt = dkitopt('FrequencyVector',omega);
[K,clp,bnd,dkinfo]  = dksyn(Punc,nmeas,nu,opt)