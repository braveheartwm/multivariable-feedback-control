% Example 4. 5 Controllability of tanks in series.
% Original file from:
% ------------------------------------------------------------------------------
% Simexample.m
% Roy Smith:	very early AM, Thurs, Oct 6th, 1994.
% Implement Sigurd's suggested four tank simulation to show 
% that being state controllable doesn't get you everything.
% Tanks are modeled by a single lag, with unity gain, and are cascaded.
% The input flow is constant, and we vary only the input temperature.
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl4_5.m,v 1.4 2004/04/02 17:33:13 vidaral Exp $
% ------------------------------------------------------------------------------
clear all; close all;

tc = 0.01; % 100 seconds time constant.

%The state-space realtion
A=[-tc,   0,   0,   0;
    tc, -tc,   0,   0;
     0,  tc, -tc,   0;
     0,   0,  tc, -tc];

B=[tc; 0; 0; 0];
C=eye(4,4);
D=0;
fourtanks = ss(A,B,C,D);

%	Check the controllability of the system.
Co = ctrb(A,B);
rCo = rank(Co);
disp(['rank controllability matrix = ' int2str(rCo)])

% Now specify the desired response after h seconds.  
% The initial condition is zero.
% If h is set the large the numerical errors
% coming from the simulation of an unstable system 
% swamp us.

% We allow negative temperatures as we assume that 
% this has been linearized about some point.
t0=0;
tfinal=400;
dt=0.5;
T=t0:dt:tfinal;
Ttotal=t0:dt:500;        %to plot after the final state
Tdiff=tfinal+dt:dt:500;
xh = [1; -1; 1; -1];	% desired temperatures

%calculate the required input from (4.43)
eA = expm(-A*tfinal);
Q = eA*B*B'*eA' - B*B';
X = lyap(A,Q);

% calculate the required input,uspec, by simulating the 
% adjoint, from the final state as the initial condition.

Cadj = B';
Aadj = -A';
Badj = ones(4,1);		% dummy B.  input will be zero anyway.
xadj0 = inv(X)*eA*xh;
Dadj= 0;
Uadj = ss(Aadj,Badj,Cadj,Dadj);

%make it go to zero after the desired input.
x0=[0;0;0;0];
uspec1 = initial(Uadj,xadj0,T);
uspec2 = initial(Uadj,x0,Tdiff);
uspec=vertcat(uspec1,uspec2);
yspec = lsim(fourtanks,uspec,Ttotal);

subplot(211)
plot(Ttotal,uspec,Ttotal,0,':')
xlabel('TIME [seconds]'); ylabel('INPUT SIGNAl')
title(['INPUT TO GIVE DESIRED STATE AT  ' num2str(tfinal) ' SECONDS'])
text(370,90,'T0');

subplot(212)
plot(Ttotal,yspec)
xlabel('TIME [seconds]'); ylabel('STATES')
title(['RESPONSE OF STATES (TANK TEMPERATURES)']);
text(240,16,'T1');text(250,5,'T2');text(290,0,'T3');text(300,-7,'T4');
axis([0 500 -20 20]);

