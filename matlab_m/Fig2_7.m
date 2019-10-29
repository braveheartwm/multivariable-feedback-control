% Plant 2.26,    G(s)=3(-2s+1)/(5s+1)(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Eq2_31.m,v 1.1 2004/01/26 16:37:11 heidisi Exp $

s=tf('s');
g=3*(-2*s+1)/(5*s+1)/(10*s+1);

% frequency response of G
points=301;
w=logspace(-3,1,points);
[mag,pha]=bode(g,w); 
mag=mag(:); pha=pha(:);

% plotting Figure 2.7
subplot(2,1,1)
loglog(w,mag)
axis([1e-2 10 1e-2 10])
ylabel('Magnitude')
xlabel('\omega_{180}')
subplot(2,1,2)
semilogx(w,pha)
axis([1e-2 10 -270 0])
ylabel('Phase')
xlabel('Frequency [rad/s]')
