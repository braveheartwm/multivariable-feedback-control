%Example 2.17, Produce Figures 2.30 and 2.31.

%Plant (2.62):  
%  G(s)=200/(10s+1)(0.005s+1)^2 ; Gd(s)=100/(10s+1)
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Expl2_17.m,v 1.2 2004/01/30 13:12:17 espensto Exp $

clear

s=tf('s');
Eq2_62

time=0.01:0.01:3; omega=logspace(-2,2,101);

%weight Wp1
% Require  wb, slope 1 and peak less than 1.5 
M=1.5; wb=10; A= 1e-4; Wp1 = tf([1/M wb], [1 wb*A]); 

%weight Wp2
% Require wb, slope 2 and peak less than 1.5 
M=1.5; wb=10; 
Wp2=(s/sqrt(M)+wb)^2/(s+wb*sqrt(A))^2;
Wu = 1;
systemnames = 'G Wp Wu';
inputvar = '[ r(1); u(1)]';
outputvar = '[Wp; Wu; r-G]';
input_to_G = '[u]';
input_to_Wp = '[r-G]';
input_to_Wu = '[u]';
sysoutname = 'P';
for i=1:2,
  if i==2,
    cleanupsysic = 'yes';
  else
    cleanupsysic = 'no';
  end
  eval(['Wp=Wp',int2str(i),';']);
  sysic;
  nmeas=1; nu=1; gmn=0.667; gmx=20; tol=0.001;
  [K5,ghinf,gopt] = hinfsyn(P,nmeas,nu,'GMIN',gmn,'GMAX',gmx,...
                            'TOLGAM',tol,'DISPLAY','on');
  eval(['K5',int2str(i),'=K5;gopt',int2str(i),'=gopt;']);
end

% RESULTS of H-infinity:
disp(sprintf('%10s%6s%6s%6s%6s%6s%6s','','H_inf','Ms','Mt','GM','PM','Wc'));
for i=1:2,
  eval(['K5=K5',int2str(i),';gopt=gopt',int2str(i),';']);
  L=G*K5;
  S=1/(1+L);
  T=1-S;
  SGd=S*Gd;
  y5=step(SGd,time);
  y5r=step(T,time);
  eval(['y5',int2str(i),'=y5;y5r',int2str(i),'=y5r;S5',int2str(i),'=S;']);
  [Gm,Pm,W180,Wc]=margin(L);
  MS=norm(S,inf,1e-4);
  MT=norm(T,inf,1e-4);
  disp(sprintf('Design%2d: %6.2f%6.2f%6.2f%6.2f%6.1f%6.2f',...
      i,gopt,MS,MT,Gm,Pm,Wc));
end

figure(1); 				%Figure 2.30
w1=1/Wp1; [magw1,phaw1]=bode(w1,omega);
w2=1/Wp2; [magw2,phaw2]=bode(w2,omega);
[mags1,phas1]=bode(S51,omega);
[mags2,phas2]=bode(S52,omega);

%subplot(111)
loglog(omega,mags1(:),'b',omega,mags2(:),'r',omega,magw1(:),...
       '--',omega,magw2(:),'--',omega,omega./omega,':')
axis([.01,100,.0001,10]);
xlabel('Frequency');ylabel('Magnitude')
title('Sensitivity and Inverse of performance weight');
text(0.02,0.0009,'1 / Wp1');text(0.4,.001,'1 / Wp2');
text(0.02, 0.007,'S1');text(.3,.007,'S2');

figure(2); %Figure 2.31
subplot(221)
plot(time, [y5r1,y5r2],time,time./time,':')
axis([0,3,-0.2,1.5]);
xlabel('Time');ylabel('y');
text(0.4,1.35,'y_HINF_2');text(0.4,0.9,'y_HINF_1');
title('Tracking response');
subplot(222)
plot(time, [y51,y52],time,time./time,':')
axis([0,3,-0.2,1.5]);
xlabel('Time');ylabel('y');
title('Disturbance response');
text(0.4,0.5,'y_HINF_2');text(1.4,1.0,'y_HINF_1');
