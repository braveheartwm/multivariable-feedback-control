%Figure 2.7 Bode plots
%(explained in Example 2.2)

% Plant (2.31), G(s)=3(-2s+1)/(5s+1)(10s+1)
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig2_7.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

clear
s=tf('s');
g=3*(-2*s+1)/(5*s+1)/(10*s+1);

% frequency response of G
points=301;
w=logspace(-3,1,points);
[mag,pha]=bode(g,w); 
mag=mag(:); pha=pha(:);

[gm,pm,wcg,wcp]=margin(g);   % wcg=0.4125

clf;
subplot(211);
loglog(w,mag,[1e-2 wcg],[0.4 0.4],'g:',[wcg wcg],[1e-2 0.4],'g:')
axis([0.01 10 0.01 10]);
text(.015,.7,'0.4');
ylabel('Magnitude')
xlabel('     \omega_{180}')
subplot(212);
semilogx(w,pha,[1e-2 wcg],[-180 -180],'g:',[wcg wcg],[-180 0],'g:')
axis([0.01 10 -300 0]);
text(.015,-160,'-180');
ylabel('Phase')
xlabel('Frequency [rad/s]')




