% Figure 8.16 Performance and Uncertainty weight
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig8_16.m,v 1.2 2004/04/14 13:16:22 vidaral Exp $

% Performance weight (8.132)
wp = tf([0.5 0.05],[1 1.e-5]); % Approximated integrator
% Uncertainty weight
wi = tf([1 0.2],[0.5 1]);
% Frequency vector
omega = logspace(-3,2,51);
% Plot weights
bodemag(wp,wi);hold on
xlabel('Frequency'); ylabel('Magnitude'); 
text(.03,3,'W_P'), text(10,3,'W_I');
axis([0.001 100 .1 100]);
plot([1e-3 100],[1 1],'k:');
