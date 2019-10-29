%Figure 5.9 Non-causal controllers (code to check example on page 187/188)
% 
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Fig5_9.m,v 1.3 2004/02/03 14:12:00 vidaral Exp $

clear all; close all;

z      = 1;		% Constant z in transfer function.
t_step = 0;		% Unit step at time = t_step.
Gnum   = [-1 z];	% Numerator of transfer function. 
Gden   = [ 1 z];	% Denominator of transfer function.

t = [-5:0.01:5];	% Time vector.
n = length(t);		% Number of time steps.
m = find(t==t_step);	% Find index at step, i.e. step at t(m).

% GENERATE REFERENCE.
r = [zeros(1,m) ones(1,n-m)]; 

% GENERATE u FROM UNSTABLE AND NON-CAUSAL CONTROLLERS.
u_unstable  = [zeros(1,m) 1-2*exp(z*(t(m+1:n)-t(m)))];
u_noncausal = [2*exp(z*(t(1:m)-t(m))) ones(1,n-m)];
u_real=r;

% SIMULATE USING lsim (in Control Systems Toolbox).
[y_unstable,x]  = lsim(Gnum,Gden,u_unstable,t);
[y_noncausal,x] = lsim(Gnum,Gden,u_noncausal,t);
[y_real,x] = lsim(Gnum,Gden,u_real,t);

% PLOT RESULTS.
figure(1)
[Pp] = get(gcf, 'PaperPosition');
set(gcf, 'Units','inches')
set(gcf, 'PaperPosition', [Pp(1) Pp(2) Pp(3) Pp(4)/1.2]);
set(gcf, 'Position', [1 7.0 Pp(3) Pp(4)/1.2]);


subplot(3,2,2)
plot(t,y_unstable)
axis([min(t) max(t) -1 1.1])
title('Output unstable controller:')
xlabel('')

subplot(3,2,1)
plot(t,u_unstable)
axis([min(t) max(t) -10 1])
title('Input unstable controller:')
xlabel('')

subplot(3,2,4)
plot(t,y_noncausal)
axis([min(t) max(t) -1 1.1])
title('Output non-causal controller:')
xlabel('')

subplot(3,2,3)
plot(t,u_noncausal)
axis([min(t) max(t) -0.1 2.1])
title('Input non-causal controller:')
xlabel('')

subplot(3,2,6)
plot(t,y_real)
axis([min(t) max(t) -1 1.1])
title('Output practical controller:')
xlabel('Time [sec]')

subplot(3,2,5)
plot(t,u_real)
axis([min(t) max(t) -0.1 2.1])
title('Input practical controller:')
xlabel('Time [sec]')

[Ap] = get(gca, 'Position');
