%Model (1.76)
%  G(s)=200/(10s+1)(0.05s+1)^2
% Gd(s)=100/(10s+1)
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Eq2_62.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

s=tf('s');
Gd=100/(10*s+1);
G=200/(10*s+1)/(0.05*s+1)^2;