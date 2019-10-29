% Figure 5.6: "Ideal" sensitivity function for a plant with delay
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig5_6.m,v 1.4 2004/02/03 14:12:00 vidaral Exp $
clear all; close all;

% Plot L when the ideal T is a delay
w1 = logspace(-2, 0, 101);
w2 = logspace(0, 1, 1001);
w3 = logspace(1, 2, 10001);
omegax=[w1 w2 w3];
theta=1;

s_g=1-exp(-theta*j*omegax');
z(1:length(omegax'))=1;
mags_g=abs(s_g);

l1=3e-2:1e-1:1/theta;
one(1:length(l1))=1;

loglog(omegax',mags_g,omegax',z,':k',one,l1,':k');
axis([.01 100 .01 10]);
xlabel('Frequency x Delay');
ylabel('Magnitude |S|');
text(1/theta,1.5e-2,'\omega = 1/\theta',     'VerticalAlignment', 'bottom',...
                           'HorizontalAlignment', 'center');
