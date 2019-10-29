%Table 2.1, Peak values and total variation
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Tab2_1.m,v 1.1 2004/01/26 16:37:11 heidisi Exp $

clear all
z=[2.0;1.5;1.0;0.8;0.6;0.4;0.2;0.1;0.01];
tau=1;
disp(sprintf('%12s%14s%10s%12s%12s\n','Zeta','Overshoot','TV ','Mt ','Ms '));
for i=1:9,
  zeta=z(i);
  peaks;
  disp(sprintf('%12.2f%12.2f%12.2f%12.2f%12.2f',zeta,Ov,Tv,Mt,Ms));
end
