% Example 7.6  
% Multiplicative weight for parametric uncertainty (Figure 7.6)
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl7_6.m,v 1.3 2004/04/14 09:37:41 vidaral Exp $
close all; clear all;
w=logspace(-2,1,151);
T=4; num=[T 0.2]; den=[T/2.5 1];
wI1=tf(num,den);
[mag_wI1,pha_wI1]=bode(wI1,w); loglog(w,mag_wI1(:),'linewidth',2); % Plot
axis([0.01 10, 0.1 3]);
xlabel('Frequency'); ylabel('Magnitude');
hold on
disp('This is the first weight we tried')
drawnow
% We will find that this weight does NOT include all of the eight plants
% considered (the worst case is for case 5 (k_max, \tau_min, \theta_min)

% Nominal model with no delay
k0=2.5; tau=2.5;
g0=tf(k0,[tau 1]); frsp_g0=freqresp(g0,w);
% Generate sampled plants
pars=[2 2.5 3];
for i1=1:3 %uncertainty \theta
   theta=pars(i1);
    for i2=1:3 % uncertainty k
        k=pars(i2); 
        for i3=1:3 %uncertainty \tau
            tau=pars(i3);
            g=tf(k,[tau 1],'OutputDelay',theta); frsp_g=freqresp(g,w);
            lI=abs((frsp_g(:)-frsp_g0(:))./frsp_g0(:)); % Calculating the bound
            loglog(w,lI(:),':');
            disp(sprintf('This is plant with theta=%4.1f,k=%4.1f',theta,k));
            drawnow;
        end
    end
end
% Corrected weight
th=3/3; z1=0.8; z2=0.7;
wI2=tf([th*th 2*z1*th 1],[th*th 2*z2*th 1]);
wI=wI1*wI2; [mag_wI,pha_wI]=bode(wI,w);
disp('Now, the modified weight');
loglog(w,mag_wI(:),'--','linewidth',2)
hold off
