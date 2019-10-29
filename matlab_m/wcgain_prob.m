%Worst case sensitivity for input uncertainty

G0=[87.8 -86.4; 108.2 -109.6];
G=tf([1],[75 1])*G0;
G=minreal(ss(G)); % ss

% Inversed based controller
Kinv=0.7*tf([75 1],[1 1e-3])*inv(G0); % Integrator moved to 1e-3

Wp=0.5*tf([10 1],[10 1e-3])*eye(2); %weight for performance, Integrator moved to 1e-3
Wi=tf([1 0.2],[0.5 1])*eye(2); %weight for uncertainty

%
% OLD method, uses muperf
%
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
N=lft(P,Kinv);
% 
% Nsys=ltisys(N.a,N.b,N.c,N.d,1);  
% delta=ublock(2,1);
% [muup,muupfreq]=muperf(Nsys,delta);   % correct solution = 44.883 at frequency 1.4384

% Using Robust Control Toolbox

udel = [ultidyn('udel',[2 2])]; % uncertainty

Np = lft(udel,N); %Perturbed model

NwLTI = wcgain(Np);  %(ans lower bound 21.0505 at frequency 2.0000e-006)

% Try Freuqnecy domain

omega = logspace(0,0.2,100);      % Very fine gridding between 1 and 1.58
                                  % Note solution is at 1.4384
Npf = frd(Np,omega);
NwFreq = wcgain(Npf);  %(ans lower bound 41.9682 at frequency 1)

% Try directly at frequency 1.4384
Nf = frd(N, 1.4384);
uconst = ucomplexm('utemp',[0 0;0 0]); % uncertainty
NwExactFreq = wcgain(lft(uconst,Nf));  %(ans lower bound 44.8807)
                                  
%%%%%%%%%
% Possible reason - pole zero cancellation between G and K (pole at s = 1/75, 0)
% CRemove pole at 1/75 and change pole in Wp to 10e-4, so no cancellation
% with pole of controller at 10e-3

G1=minreal(ss(G0)); % ss

% Inversed based controller
Kinv1=0.7*tf(1,[1 1e-3])*inv(G0); % Integrator moved to 1e-3

Wp1=0.5*tf([10 1],[10 1e-4])*eye(2); %weight for performance, Integrator moved to 1e-3

% Generalized plant P
systemnames = 'G1 Wp1 Wi';
inputvar    = '[ydel(2); w(2) ; u(2)]';
outputvar   = '[Wi ; Wp1 ; -G1-w]';
input_to_G1  = '[u+ydel]';
input_to_Wp1 = '[G1+w]';
input_to_Wi = '[u]';
sysoutname  = 'P1';
cleanupsysic= 'yes';
sysic;
N1=lft(P1,Kinv1);

NwTry = wcgain(lft(udel,N)); %same answer




