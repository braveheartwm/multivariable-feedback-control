%Example 6.3
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl6_3.m,v 1.4 2004/04/14 08:13:34 vidaral Exp $

clear all; close all;
format short e

% Define zero system.
RHPZ = 2;
n1   = tf( [1 -RHPZ], [0.1 1]);
n2   = tf( [1  RHPZ], [0.1 1]);
Nz=append(n1,n2);

% Define pole system, this system is unstable.
RHPP = 3; 
d1   = tf( [1], [1 -RHPP]);
d2   = tf( [1], [1  RHPP]);
Dp = append(d1,d2);

% Setup for H-infinity
Wu = eye(2);   % Equal weight on both manipulated inputs.
epsv = 1e-4; wb = 0.5; M = 2;
wp = tf( [1/M wb], [1 epsv]);
Wp = append(wp,wp);
Wp = minreal( Wp );
nmeas = 2; nu = 2; gam = 0.5; gmx = 20; tol = 0.001;


%**************  ALL ANGLES ********************
angles = [0 30 60  90];
Yz=[];
Phi=[];
Cc=[];
Snorm=[];
Tnorm=[];
Gamma=[];
Wp=ss(Wp);
Wu=ss(Wu);
%    systemnames = 'G Wp Wu';
%    inputvar    = '[r(2); u(2)]';
%    outputvar   = '[Wp; Wu; r-G]';
%    input_to_G  = '[u]';
%    input_to_Wp = '[r-G]';
%    input_to_Wu = '[u]';
%    sysoutname  = 'Psys';


for i=1:4
% angle of rotation matrix
angle = angles(i);

% Rotation matrix
za = angle*pi/180;
Tz = [cos(za) -sin(za)
      sin(za)  cos(za) ];

% overall plant
G = Dp*Tz*Nz;
[A,B,C,D]=ssdata(G);  

% zero direction
[U,Sv,V]=svd(C*inv(RHPZ*eye(size(A))-A)*B+D);
yz = U(:,2); 
Yz=[Yz yz];

% pole direction
epsp1 = 0.001;  %should try different values and make sure yp is the same
[U,Sv,V]=svd(C*inv((RHPP+epsp1)*eye(size(A))-A)*B+D);
yp =  U(:,1);

% Compute inner product and angle
inp = yz'*yp;
phi = acos(abs(inp))*180/pi;
Phi = [Phi phi];

% Compute c according to (6.8)
c=sqrt((1-inp^2)+(RHPZ+RHPP)^2/(RHPZ-RHPP)^2*inp^2);
Cc=[Cc c];

% % H-infinity controller design
% if i==4,
%    cleanupsysic= 'yes';
% else
%    cleanupsysic= 'no';
% end
%    sysic;
% [K, clp, gamfin] = hinfsyn(Psys, nmeas, nu);
[K,clp,gamfin] = mixsyn(G,Wp,Wu,[]);
Gamma=[Gamma gamfin];

% Perform simulation
L=G*K;
S=inv(eye(2)+L);
T=eye(2)-S;

timestep=0:0.01:10;
U1=ones(1,size(timestep,2));
U2=-ones(1,size(timestep,2));
U=vertcat(U1,U2);
[y,time]=lsim(T,U,timestep);

ms=norm(S,inf);
mt=norm(T,inf);
Snorm=[Snorm ms];
Tnorm=[Tnorm mt];

%Plot response
subplot(2,2,i)
plot(time,y(:,1),time,y(:,2),time,-1,':k',time,1,':k')
axis( [0 5 -2 2] )
if i>2
 xlabel('Time');
end

end

%Table on p.260
disp(sprintf('\nAlpha%6d%6d%6d%6d',angles));
disp(sprintf('  Yz %6d%6.2f%6.2f%6.0f\n     %6d%6.2f%6.2f%6d',Yz'));
disp(sprintf(' Phi %6d%6.1f%6.1f%6.0f',Phi));
disp(sprintf('   c %6.1f%6.2f%6.2f%6.1f',Cc));
disp(sprintf('  Ms %6.2f%6.2f%6.2f%6.2f',Snorm));
disp(sprintf('  Mt %6.2f%6.2f%6.2f%6.2f',Tnorm));
disp(sprintf('Gamma%6.2f%6.2f%6.2f%6.2f',Gamma));

%Figure 6.1
subplot(2,2,1), title('phi = 0 degrees'), 
subplot(2,2,2), title('phi = 71 degrees'), 
subplot(2,2,3), title('phi = 83 degrees'), 
subplot(2,2,4), title('phi = 90 degrees'), 