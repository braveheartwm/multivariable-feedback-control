% Example 7.9 Produces fig. 7.13
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl7_9.m,v 1.2 2004/04/14 09:37:41 vidaral Exp $
close all; clear all;

% Frequency
omega=logspace(-3,1,201);

% Simple second order plant with RHP-zero at z=2.
G=3*tf([-2 1],[50 15 1]);

% PI controller
Kc1=1.13; ti=12.7;
K1=Kc1*tf([ti 1],[ti 0]);

% Pertubed plant
Gp=4*tf([-3 1],[16 8 1]);

% Frequency response of G and Gp
Gf=freqresp(G,omega);
Gpf=freqresp(Gp,omega);

% Relative differense
li=(Gp-G)/G;

% Uncertainty weigth
wi=tf([10 0.33],[10/5.25 1]);
wif=freqresp(wi,omega);

% Plot
figure(1)
set(cstprefs.tbxprefs,'MagnitudeUnits','log','MagnitudeScale','log') % Setting plot types
bodemag(li,'r',wi,'y--',tf(1,1),'r--',omega);
legend('Relative error','Uncertainty weight')

% Evaluate the complementary sensitivity
L1=G*K1; S1=1/(1+L1); T1=1-S1;
disp('Closed-loop poles of T1 (nominal plant, stable):');disp(pole(T1));
figure(2)
bodemag(1/wi,'--',T1,'r-',tf(1,1),'r--',omega); % Not RS
text(.002,4,'1 / wI'); text(0.2,4,'T1 (NOT RS)');
disp('Check with actual pertubed plant');
L1p=Gp*K1; S1p=1/(1+L1p); T1p=1-S1p;
disp('Closed-loop poles of T1p (perturbation plant, unstable):');disp(pole(T1p));

% Need to reduce the controller gain
Kc2=0.3103; ti=12.7; 
K2=Kc2*tf([ti 1],[ti 0]);
L2=G*K2; S2=1/(1+L2); T2=1-S2;
figure(3)
bodemag(1/wi,'--',T1,'r-',T2,'b-',tf(1,1),'r--',omega)
axis([.001,10,.03,7]);
text(.002,4,'1 / wI'); text(0.2,4,'T1 (NOT RS)');
text (0.2,.07,'T2 (RS)');
