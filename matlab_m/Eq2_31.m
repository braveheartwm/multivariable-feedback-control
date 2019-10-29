% Plant 2.31,    G(s)=3(-2s+1)/(5s+1)(10s+1)
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Eq2_26.m,v 1.1 2004/01/26 16:37:11 heidisi Exp $

s=tf('s');
g=3*(-2*s+1)/(5*s+1)/(10*s+1);

% frequency response of G
points=301;
w=logspace(-3,1,points);
[mag,pha]=bode(g,w); 
mag=mag(:); pha=pha(:);
