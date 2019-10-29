%Table 2.4. Sample MTALAB program to synthesize H-infinity controller
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% $Id: Tab2_3.m,v 1.2 2004/10/13 09:28:08 kariwala Exp $

clear
s=tf('s');  
G=200/(10*s+1)/(0.05*s+1)^2; % plant is G
M=1.5; wb=10; A= 1e-4; Wp = tf([1/M wb], [1 wb*A]); Wu=1; %weights

% Uncomment to try other weight
%Wp = tf([1/M^0.5 wb], [1 wb*A^0.5]); Wp = Wp*Wp;

%
% Generalized plant P is found with function sysic
%
% systemnames = 'G Wp Wu';
% inputvar = '[ r(1); u(1)]';
% outputvar = '[Wp; Wu; r-G]';
% input_to_G = '[u]';
% input_to_Wp = '[r-G]';
% input_to_Wu = '[u]';
% sysoutname = 'P';
% cleanupsysic = 'yes';
% sysic;
%
% Find H-infinity optimal controller
%
% nmeas=1; nu=1; gmn=0.5; gmx=20; tol=0.001;
%   [khinf,ghinf,gopt] = hinfsyn(P,nmeas,nu,'GMIN',gmn,'GMAX',gmx,...
%                             'TOLGAM',tol,'DISPLAY','on');
[khinf, ghinf, gopt] = mixsyn(G,Wp,Wu,[]);
Marg = allmargin(G*khinf);

Ms = norm(inv(1+khinf*G),inf)
Mt = norm(khinf*G*inv(1+khinf*G),inf)
[Gm,Pm,W180,Wc]=margin(khinf*G)
