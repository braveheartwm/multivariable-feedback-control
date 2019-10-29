% Example 7.7 produces Figure 7.7 Uncertainty Regions
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl7_7.m,v 1.2 2004/04/14 09:37:41 vidaral Exp $
close all; clear all;
w=0.5;
% Uncertainty region:
gf=uncreg(w);
plot(real(gf),imag(gf),'--');
grid;axis('square');hold on;

k=2.5;tau=2.5;theta=2.5;
% Norminal Plant option 1:
gf1=k/(j*w*tau+1);
% Uncertainty radius
r1=max(abs(gf1-gf));
phi=0:0.1:2*pi+0.1;
uc1=gf1+r1*exp(-j*phi);
plot(real(uc1),imag(uc1),'-',real(gf1),imag(gf1),'+');

% Norminal Plant option 2:
gf2=gf1*exp(-j*w*theta);
% Uncertainty radius
r2=max(abs(gf2-gf));
uc2=gf2+r2*exp(-j*phi);
plot(real(uc2),imag(uc2),'-',real(gf2),imag(gf2),'o');

% Norminal Plant option 3:
% Optimal plant with minimized radius:
x=fminsearch(@maxrad,[real(gf2);imag(gf2)],[],gf);
gf3=x(1)+i*x(2);
% Uncertainty radius
r3=max(abs(gf3-gf));
uc3=gf3+r3*exp(-j*phi);
plot(real(uc3),imag(uc3),'-',real(gf3),imag(gf3),'*');
