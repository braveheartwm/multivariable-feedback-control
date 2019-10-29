% Figure 2.3    Bode plots
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Fig2_3.m,v 1.1 2004/01/26 16:37:11 heidisi Exp $

clear
gain=30;
s=tf('s');
gain=30;
L=gain*(s+1)/(s+0.01)^2/(s+10);

w=logspace(-3,3,100);
[mag,pha]=bode(L,w);

clf;
subplot(211)
loglog(w,mag(:),[1e-3 1e-2 1 10 1e3], [3e4 3e4 3 0.3 3e-5],'g:',...
       [1e-2 1e-2],[1e-5 3e4],'g:',[1 1],[1e-5 3],'g:',[10 10],[1e-5 0.3],'g:') 
ylabel('Magnitude')
xlabel('Frequency')
axis([1e-3,1e3,1e-5,1e5])
text(0.002,5000,'0');text(0.06,100,'-2');text(3,0.1,'-1');text(60,0.001,'-2');

subplot(212) 
semilogx(w,pha(:),[1e-3 1e-2 1e-2 1 1 10 10 1e3],...
                  [0 0 -180 -180 -90 -90 -180 -180],'g:')

ylabel('Phase')
xlabel('Frequency')
axis([1e-3,1e3,-200,0]);
