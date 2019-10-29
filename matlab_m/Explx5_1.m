% Example5_1.m: Extra example for chapter 5. After remark p. 167
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Explx5_1.m,v 1.3 2004/02/03 14:12:00 vidaral Exp $
clear all; close all;

% 1. System L1 = 2/s(s+1). Illustrates Bode's sensitivity integral
w = logspace(-1.5,1.5,121);
L1=tf(2,[1 1 0]);
S1=1/(L1+1);
%The areas above and below 1 are equal');
% We see that the area at high frequencies also contributes so
% the peak itself need not be very large.

% 2. System L = 2(-0.2s+1) / s(s+1)(0.2s+1)
% |L| is same as above but we have added term (-0.2s+1)/(0.2s+1)
L2=tf([-0.4 2],[0.2 1.2 1 0]);
S2=1/(L2+1);
%[mag,phase]=bode(S2,w);

% Plot both systems
figure(1)
bodemag(S1,'--',S2,'-',tf(1),':',w);
%semilogy(w,squeeze(mag),w,squeeze(mag1),w,z,':k');
xlabel('Frequency [rad/s]');
ylabel('Magnitude');
title('Illustration of Bode sensitivity integral: Equal areas');
axis([0 10 1e-1 1e1])
legend('L_1=2/(s(s+1))','L_2=(-0.4s +2)/(0.2s^3 +1.2s^2 +s)','BR')
text(1.2,1,'L_1'); text(0.9,6,'L_2');
% In BOTH cases the areas of ln|S| below and above 1 must be equal.
% BUT with the RHP-zero the peak must be higher,  
% see Thm.5.2 where only frequencies up to about z=5 count.
