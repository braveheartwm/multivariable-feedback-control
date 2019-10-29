% Uses the Control toolbox
num=3*[-2 1]; den=conv([5 1],[10 1]); % inverse response process
[a,b,c,d] = tf2ss(num,den); % plant is (a,b,c,d)
G = ss(a,b,c,d);
% Model dimensions:
p = size(c,1); % no. of outputs (y)
[n,m] = size(b); % no. of states and inputs (u)
Znm=zeros(n,m); Zmm=zeros(m,m);
Znn=zeros(n,n); Zmn=zeros(m,n);
% 1) Design state feedback regulator
A = [a Znm;-c Zmm]; B = [b;-d]; % augment plant with integrators
Q=[Znn Znm;Zmn eye(m,m)]; % weight on integrated error
R=eye(m); % input weight
Kr=lqr(A,B,Q,R); % optimal state-feedback regulator
Krp=Kr(1:m,1:n);Kri=Kr(1:m,n+1:n+m); % extract integrator and state feedbacks
% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = eye(n); % process noise model (Gd)
W = eye(n); V = 1*eye(m); % process and measurement noise weight
Estss = ss(a,[b Bnoise],c,[0 0 0]);
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter gain
% 3) Form overall controller
Ac=[Zmm Zmn;-b*Kri a-b*Krp-Ke*c]; % integrators included
Bcr = [eye(m); Znm]; Bcy = [-eye(m); Ke];
Cc = [-Kri -Krp]; Dcr = Zmm; Dcy = Zmm;
Klqg2 = ss(Ac,[Bcr Bcy],Cc,[Dcr Dcy]); % Final 2-DOF controller from [r y]' to u
Klqg = ss(Ac,-Bcy,Cc,-Dcy); % Feedback part of controller from -y to u
% Simulation
sys1 = feedback(G*Klqg,1); step(sys1,50); % 1-DOF simulation
sys = feedback(G*Klqg2,1,2,1,+1); % 2-DOF simulation
sys2 = sys*[1; 0]; hold; step(sys2,50);
allmargin(G*Klqg)% Uses the Control toolbox
num=3*[-2 1]; den=conv([5 1],[10 1]); % inverse response process
[a,b,c,d] = tf2ss(num,den); % plant is (a,b,c,d)
G = ss(a,b,c,d);
% Model dimensions:
p = size(c,1); % no. of outputs (y)
[n,m] = size(b); % no. of states and inputs (u)
Znm=zeros(n,m); Zmm=zeros(m,m);
Znn=zeros(n,n); Zmn=zeros(m,n);
% 1) Design state feedback regulator
A = [a Znm;-c Zmm]; B = [b;-d]; % augment plant with integrators
Q=[Znn Znm;Zmn eye(m,m)]; % weight on integrated error
R=eye(m); % input weight
Kr=lqr(A,B,Q,R); % optimal state-feedback regulator
Krp=Kr(1:m,1:n);Kri=Kr(1:m,n+1:n+m); % extract integrator and state feedbacks
% 2) Design Kalman filter % don’t estimate integrator states
Bnoise = eye(n); % process noise model (Gd)
W = eye(n); V = 1*eye(m); % process and measurement noise weight
Estss = ss(a,[b Bnoise],c,[0 0 0]);
[Kess, Ke] = kalman(Estss,W,V); % Kalman filter gain
% 3) Form overall controller
Ac=[Zmm Zmn;-b*Kri a-b*Krp-Ke*c]; % integrators included
Bcr = [eye(m); Znm]; Bcy = [-eye(m); Ke];
Cc = [-Kri -Krp]; Dcr = Zmm; Dcy = Zmm;
Klqg2 = ss(Ac,[Bcr Bcy],Cc,[Dcr Dcy]); % Final 2-DOF controller from [r y]' to u
Klqg = ss(Ac,-Bcy,Cc,-Dcy); % Feedback part of controller from -y to u
% Simulation
hold
sys1 = feedback(G*Klqg,1); step(sys1,50); % 1-DOF simulation
sys = feedback(G*Klqg2,1,2,1,+1); % 2-DOF simulation
sys2 = sys*[1; 0]; hold; step(sys2,50);
allmargin(G*Klqg)
Ms=norm(sys1,inf,1e-4)
Mt=norm(1-sys1,inf,1e-4)