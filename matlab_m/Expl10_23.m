% Example 10.23 Control structure design and performance verification
% 
% This script produces the results and plots for example 10.19; for the
% detailed-model distillation column, it executes an initial
% controllability analysis, an analysis of decentralized control, and
% finally plots the performance under the effect of disturbances.
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl10_19.m,v 1.10 2004/02/19 20:24:15 vidaral Exp $
clear all; close all;
set( cstprefs.tbxprefs, 'MagnitudeUnits', 'abs', 'MagnitudeScale', 'log' ) ;

% Linear Model
% Distillation-column model with 5 states and disturbances

A = [ -0.005131  0        0       0       0      ;
       0        -0.07366  0       0       0      ;
       0         0       -0.1829  0       0      ;
       0         0        0      -0.4620  0.9895 ;
       0         0        0      -0.9895 -0.4620 ] ;

B = [ -0.629  0.624 ;
       0.055 -0.172 ;
       0.030 -0.108 ;
      -0.186 -0.139 ;
      -1.23  -0.056 ] ;


C= [ -0.7223 -0.5170 0.3386 -0.1633 0.1121 ; 
     -0.8913  0.4728 0.9876  0.8425 0.2186 ] ;

D = [ 0 0 ;
      0 0 ] ;

Bd = [ -0.062 -0.067 ; 
        0.131  0.040 ;
        0.022 -0.106 ; 
       -0.188  0.027 ;
       -0.045  0.014 ] ;
   
% Converting minutes to seconds to get properly scaled results

A = A/60;
B = B/60;
Bd = Bd/60;

G = ss(A,B,C,D);
Gd = ss(A,Bd,C,D) ;

% Change the disturbance magnitudes from 30 to 20% for d1 and
% from 40 to 20% for d2 (this is because the scalings used in this
% example are somewhat different from that used in the paper by
% Hovd and Skogestad (Automatica, 1992) from which the model above
% was obtained). 
% The figures in the book do not use this.

%Dmag = [ 20/30     0 ; 
%             0 20/40 ] ;
%Gd = Gd * Dmag ;

[ Ad,Bd,Cd,Dd] = ssdata(Gd);

% *********************
% Steady-state analysis

G0	= -C*inv(A)*B+D
Gd0	= -Cd*inv(Ad)*Bd+Dd
RGA0 	= G0.*pinv(G0)' 
PRGA0 	= diag(diag(G0))*pinv(G0)     % Performance RGA
CLDG0 	= PRGA0*Gd0                   % Closed-Loop Disturbance Gain
GinvGd0 = inv(G0)*Gd0


% ******************
% Frequency analysis

% w = logspace(-3,1,41);
% Gf = frsp(G,w); Gdf=frsp(Gd,w);

% Drawing figure 10.15

figure(1); clf;
bodemag(Gd(1,1),'-', Gd(1,2),'-', Gd(2,1),'-', Gd(2,2),'-',tf(1),'--', {1e-3/60, 1}); 
title( 'Gd elements' );
axis([0.001/60,10/60,.1,20]);
text(.08/60,0.2,'g_{d11}');
text(.117/60,0.30,'g_{d12}');
text(0.15/60,0.4,'g_{d22}');
text(0.5/60,0.4,'g_{d21}');

% Obtaining transfer functions
Gt = tf(G);
Gdt = tf(Gd);

w = logspace(-5, 2, 301);
Gw = frd(G, w);

invtGw = inv(Gw)';

rga = [ Gw(1,1)*invtGw(1,1), Gw(1,2)*invtGw(1,2); 
        Gw(2,1)*invtGw(2,1), Gw(2,2)*invtGw(2,2)];
    
diagGt = [Gt(1,1),          0 ;
           0,        Gt(2,2)] ;

prga = diagGt*invtGw';
cldg = prga*Gd;


% Drawing figure 10.16
figure(2); clf;
bodemag(cldg(1,1),'-', cldg(1,2),'-', cldg(2,1),'-', cldg(2,2),'-',tf(1),'--', {1e-3/60,1}); 
title( 'CLDG elements' );
text(0.005/60,20,'g~_{d11}'); text(0.005/60,0.4,'g~_{d12}');
text(0.005/60,5,'g~_{d22}'); text(0.02/60,25,'g~_{d21}');
axis([0.001/60,10/60,.1,100]);

% Drawing figure 10.17
figure(3) ; clf;
bodemag(prga( 1, 1),'-', prga( 1, 2 ),'-', prga( 2, 1 ),'-', prga( 2, 2 ),'-',tf(1),'--', { 1e-3/60, 1 });
title( 'PRGA elements' );
text(4/60,1.8,'\gamma_{11}=\gamma_{22}'); text(0.1/60,1.3,'\gamma_{12}');
text(0.6/60,2,'\gamma_{21}');
axis([0.001/60,10/60,.1,100]);

% Defining two single-loop PI controllers
% Time constants are converted from minutes into seconds

k1 = tf(0.261*[3.76*60 1], [3.76*60 0]);
k2 = tf(-0.375*[3.31*60 1], [3.31*60 0]);
K = [k1, 0; 0, k2];

GK =  Gt*K;
S = inv(eye(2)+GK);
SGd = S*Gd;

St = tf(S);
t = 0:0.2:100;     % Time steps in minutes
d = [ones(250, 1)*[1 0];
     ones(251, 1)*[1 1]];
y = lsim(SGd, d, t*60);


% Drawing figure 10.18
figure(4), clf;
plot(t, y); 
xlabel( 'Time [ min ]' );
title( 'Feed composition disturbance' );
text(30,-0.25,'y2'); text(30., 0.17,'y1')
