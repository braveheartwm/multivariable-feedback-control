    
% ExplX3_1.m :  
% Exact frequency response and RGA of MIMO plant with time delay

% This is 3x3 plant where each element is on the form (k, tau, theta)
% where k is the gain, tau the time constant and theta the delay.

% For the numbers used in this example see process T4 of Luyben, 
% which is used as Example 3 by Huang et al. (J.Proc.Control, 15-27, 1994).

% We here want to use the RGA to check which pairing is best.
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: ExplX3_1.m,v 1.3 2004/04/13 12:53:00 vidaral Exp $

clear all; close all;

omega = logspace(-3,1,361);
g11 =tf([0 -1.986], [66.67 1], 'InputDelay', 0.71) ;
g12 =tf([0 5.24], [400 1], 'InputDelay', 60) ;
g13 =tf([0 5.984], [14.29 1], 'InputDelay', 2.24) ;
g21 =tf([0 0.0204], [5 1], 'InputDelay', 4.199) ;
g22 =tf([0 -0.33], [3.904 1], 'InputDelay', 1.883) ;
g23 =tf([0 2.38], [10 1], 'InputDelay', 1.143) ;
g31 =tf([0 0.374], [22.22 1], 'InputDelay', 7.75) ;
g32 =tf([0 -11.3], [35.66 1], 'InputDelay', 14.78) ;
g33 =tf([0 -9.881], [11.35 1], 'InputDelay', 1.59) ;

G=[g11 g12 g13; g21 g22 g23; g31 g32 g33]; 
Gf=freqresp(G,omega);
% Calculate RGA freq by freq, can't invert delay systems.
for i=1:length(omega); 
    Grga(:,:,i)=Gf(:,:,i).*inv(Gf(:,:,i)');
end

figure(1)
for i=1:3
  for j=1:3
    rgaij=Grga(i,j,:) ;
    rgaijabs=abs(rgaij);
    uplot('liv,lm',omega,rgaijabs(:))
    hold on
 end
end
hold off


% Pairing no. 1: Pair on diagonal elements:
i1 = [1 0 0; 0 1 0; 0 0 1];
% Compute (almost) the RGA-number, , see (3.71), but we use singular value;
% peaks of 10 for w > 0.7 rad/s
for i=1:length(omega);
    dGrga(:,:,i)=Grga(:,:,i)-i1;
end
figure (2);
dGrga=max(sum(abs(dGrga))); %inf norm
dGrga=dGrga(:);
uplot('liv,m',omega,dGrga)
% But really frequencies larger than about 0.1 rad/s are of minor importance

% Do same for pairing alternative no. 2:
i2 = [1 0 0; 0 0 1; 0 1 0];
for i=1:length(omega);
    d2Grga(:,:,i)=Grga(:,:,i)-i2;
end
figure (3);
d2Grga=max(sum(abs(d2Grga))); d2Grga=d2Grga(:); %inf norm
uplot('liv,m',omega,d2Grga)

% Conclusion: The last pairing best (consider frequencies up to 0.5 rad/s)
