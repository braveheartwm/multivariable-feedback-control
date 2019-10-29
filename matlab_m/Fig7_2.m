% Figure 7.2 
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig7_2.m,v 1.2 2004/04/14 09:37:41 vidaral Exp $
close all
w=logspace(-2,1,1001);
% inside edge
k=2;tau=3;theta=2;
Gf1=tf(k,[tau 1],'ioDelay',theta);
% outside edge
k=3;tau=2;theta=3;
Gf2=tf(k,[tau 1],'ioDelay',theta);
% Plotting edges
[Re1,Im1]=nyquist(Gf1,w);[Re2,Im2]=nyquist(Gf2,w);
plot(Re1(:),Im1(:),'r--',Re2(:),Im2(:),'r--'); hold on
axis('square');

% Given frequencies.
W=[0.01 0.05 0.2 0.5 1 2 7];
for i=1:length(W),
  w=W(i);
  gf=uncreg(w); % Uncertainty region
  plot(real(gf),imag(gf),'-');drawnow;
end
