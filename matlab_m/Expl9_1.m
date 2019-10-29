% Uses the control tool box
% Plant model
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl9_1.m,v 1.2 2004/04/15 08:10:13 vidaral Exp $

clear all; close all;
num = 3*[-2 1]; den = conv([5 1],[10 1]);  % original plant
[a,b,c,d]=tf2ss(num,den);

[n,m]=size(b);   		% finds state-vector and control sizes
p    =size(c,1); 		% finds output vector size

% First design state feedback controller
A = [a zeros(n,m) ; -c zeros(m,m)]; B = [b ; -d]; % augment plant with integrators
Q=[zeros(n,n) zeros(n,m); zeros(m,n) eye(m,m)];   % weight on integrated error 
R=1*eye(m);                   % Input weight - R small gives faster response
Kr=lqr(A,B,Q,R);            % Design optimal state-feedback regulator ;
Krp = Kr(1:m,1:n); Kri = Kr(1:m,n+1:n+m); % Extract integrator and state feedbacks:

% Now design Kalman filter for plant
Bnoise = eye(n);         % Process noise (disturbance) enters directly on states
W      = eye(n);         % Process noise weight
V      = 1*eye(m);       % Measurement noise weight - decreasing makes responsei
			 % faster
			 
[Ketf, Ke] = kalman(ss(a,[b Bnoise],c,[0 0 0]),W,V); % no need to estimate integrator states
%Ke = lqe(a,Bnoise,c,W,V); % obsolete command

% Form 2DOF controller from  [r y]' to u with integrators included
Ac = [ zeros(m,m) zeros(m,n); -b*Kri   a-b*Krp-Ke*c];  % integrators included
Bc = [eye(m) -eye(m); zeros(n,m) +Ke];
Cc = [-Kri -Krp]; Dc = [zeros(m,m) zeros(m,m)];
K = ss(Ac,Bc,Cc,Dc);

Gsim = ss(a,b,[c;c],[d;d]);
Ksim = ss(Ac,Bc,[Cc;Cc],[Dc;Dc]);
sys = lft(Gsim,-Ksim,1,1);

[out,T]=step(sys,0:50);
plot(T,out(:,1),'-',T,out(:,2),'--'); 

% Now form closed loop transfer function from [r] to [y u]'
% (bearing in mind that Dc=0)
% partition controller B-matrix to get Bc1 (r) and Bc2 (y)

Bc1=Bc(:,1:m); Bc2=Bc(:,m+1:m+m);
Cc
Acl = [a b*Cc ; Bc2*c Ac+Bc2*d*Cc ];
Bcl = [zeros(n,m) ; Bc1 ];
Ccl = [c d*Cc ; zeros(m,n) Cc]; % Augment to include control signal u
Dcl = [zeros(m,m); zeros(m,m)];

% Setup simulation

t=linspace(0,50,1000);
hold on;
for j=1:m
	y=step(Acl,Bcl,Ccl,Dcl,j,t);Cc
	plot(t,y(:,j),'-',t,y(:,m+j),'--',[1 50],[1 1],':'); 
	xlabel('Time');
end;
hold off;
axis([0 50 -0.5 2.5])
text(41,1.2,'y')
text(41,0.5,'u')
