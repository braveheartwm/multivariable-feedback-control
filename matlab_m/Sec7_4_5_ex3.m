% Section 7.4.5 example 3 Produces fig. 7.9 Uncertainty weight
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl7_5.m,v 1.3 2004/04/14 09:37:41 vidaral Exp $
clear all; close all;
w=logspace(-2,2,51);
kmax=3;kmin=2;
k_avg=(kmax+kmin)/2; rk=((kmax-kmin)/2)/k_avg;
theta_max=3;
wi=tf([(1+rk/2)*theta_max rk],[theta_max/2 1]);
wif=frd(wi,w);
w1=w(w<pi/theta_max);
w2=w(w>=pi/theta_max);
li=sqrt(rk*rk+2*(1+rk)*(1-cos(theta_max*w1)));
li=[li ones(size(w2))*(2+rk)];
lif=frd(li,w);

% Figure 7.9
set(cstprefs.tbxprefs,'MagnitudeUnits','abs','MagnitudeScale','log');
% Becuse dB is default
% add this line in matlabrc.m in \toolbox\local to set this at startup. 
% This way you dont have to do this every time bode plots are made
bodemag(lif,'-',wif,'--');
xlabel('Frequency');ylabel('Magnitude');
text(0.5,2,'li');
text(1,1.5,'wi');
axis([0.01 100 0.1 10])

