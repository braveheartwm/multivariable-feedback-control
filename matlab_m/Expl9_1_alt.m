%
% Script: Expl9_1_alt.m
%
% Author: MC Turner
%         Dept. of Engineering
%         Univ. of Leicester
%
% Date: 13th June 2001
%
% Purpose:
%
% This file solves an LQG problem incorporating integral action.
% The integrator states are not estimated
% Thanks to Carl-Fredrik Lindberg, Uppsala University
%
% Notes to users:
%
%     - includes example 9.1 system
%     - written for Matlab 4.2
%     - weights hard-wired
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp(' General LQG + integral action function ');
disp(' ');
eg=input(' Type 1 for inverse response example ');
disp(' ');

if eg==1

    % Form 'inverse response' plant

    num = 3*[0 -2 1];
    den = conv([5 1],[10 1]);

   [a,b,c,d]=tf2ss(num,den);

else

    disp(' ');
    disp(' Using state-space model (a,b,c,d) in work space... ');
    disp(' ');

end;

% Assumes plant has state-space representation
% xdot=ax + bu
%   y =cx + du

[n,m]=size(b);   % finds state-vector and control sizes
p    =size(c,1); % finds output vector size

if p==m

    % Augment plant with integrator, viz:
    %
    % xdot = ax + bu
    %   y  = cx
    %   xidot = r-cx
    %   u  = xi + utilde
    %
    % This implies
    %
    % [xdot  ] = [ a 0] [x ] + [b ] utilde + [0 ]r
    % [xidot ] = [-c 0] [xi] + [-d ]        +[1 ]
    %
    % Note that the integrator is of order 'm'

    A = [a zeros(n,m) ; -c zeros(m,m)];
    B = [b ; -d];

    % Form weights (weight integrator error only)

    Q=[zeros(n,n) zeros(n,m); zeros(m,n) eye(m,m)]; % Weight on output
    R=eye(m);                   % Weight on input - smaller R gives faster response

    % Design LQR state-feedback;

    Kr=lqr(A,B,Q,R);

    % Extract integrator and plant state feedbacks:

    Krp = Kr(1:m,1:n);
    Kri = Kr(1:m,n+1:n+m);

    % Now design Kalman filter for plant
    % (no need to estimate integrator state)
    %

if eg==1  
    Bnoise = eye(n); % Process noise (disturbance) enters directly on states
else
	Bnoise = input(' Enter disturbance distribution matrix Bd if required ');
       if isempty(Bnoise); Bnoise = eye(n); end;
end;

    W      = eye(n); % Process noise weight
    V      = 10*eye(m);  % Measurement noise weight - decreasing makes response faster
    Ke = lqe(a,Bnoise,c,W,V);

    % Kalman filter for plant is now given by
    %
    % xhat = a xhat +bu +Ke (y - cx)
    %
    % Next form 2DOF controller with integrators included
    %
    % xc = Ac xc + Bc [r y]'
    %  u = Cc xc + Dc [r y]'

    Ac = [ zeros(m,m) zeros(m,n);
           -b*Kri   a-b*Krp-Ke*c];

    Bc = [eye(m) -eye(m); zeros(n,m) +Ke];

    Cc = [-Kri -Krp];

    Dc = [zeros(m,m) zeros(m,m)];

    % Now form closed loop transfer function from [r] to [y u]'
	% (bearing in mind that Dc=0)

    %partition controller B-matrix to get Bc1 (r) and Bc2 (y)

    Bc1=Bc(:,1:m);
    Bc2=Bc(:,m+1:m+m);

    Acl = [a b*Cc ; Bc2*c Ac+Bc2*d*Cc ];

    Bcl = [zeros(n,m) ; Bc1 ];

    Ccl = [c d*Cc ; zeros(m,n) Cc]; % Augment to include control signal u

    Dcl = [zeros(m,m); zeros(m,m)];

    %
    % Setup simulation
    %

    t=linspace(0,50,1000);

    for j=1:m

        y=step(Acl,Bcl,Ccl,Dcl,j,t);
        plot(t,y(:,j),'-'); hold;
        plot(t,y(:,m+j),'--');hold;
        xlabel('Time [sec]');
        grid;
    end;

else

    disp(' Error: for output and input vector must be same length');

end;
