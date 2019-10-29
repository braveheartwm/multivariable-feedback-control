% Figure 2.2    Frequency response 
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Fig2_2.m,v 1.1 2004/01/26 16:37:11 heidisi Exp $

clear
num=5;den=[10 1];delay=2; 
G=tf(num, den,'OutputDelay',delay); %Input- and OutputDelay equiv. in SISO
w=logspace(-3,1,41)';
[mag,pha]=bode(G,w);
[gm,pm,wcg,wcp]=margin(G); %Wc
clf
subplot(211),loglog(w,mag(:),w,w./w,':',[wcp wcp],[0.01 1],':')
ylabel('Magnitude of g(s)')
xlabel('Frequency')

subplot(212), semilogx(w,pha(:),w,w./w*-180,':')
ylabel('Phase of g(s)')
xlabel('Frequency')
axis([.001,10,-300,0])