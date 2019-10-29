%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 3.7.1 Spinning satellite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Sec3_7_1.m,v 1.3 2004/02/09 12:01:24 vidaral Exp $

%State-space model (3.78)
a = 10;
A = [0  a; -a 0]; B = [1 0 ; 0 1]; C = [1 a; -a 1];  D = [0 0; 0 0];
G=ss(A,B,C,D);

%Figure 3.6 (a)
clf;
w = logspace(-2,2,120);
sigma(G,w);     % plotting singular values
text(3e0,4e-1,'\sigma_{min}(S)');
text(3e0,1e1,'\sigma_{max}(S)');
I2 = eye(2);
K=I2;

S=inv(eye(2)+G*K); % Sensitivity
T=eye(2)-S;	% Complementary Sensitivity

H_inf_S=norm(S,'inf');
H_inf_T=norm(T,'inf');
disp(sprintf('H_inf norm S: %0.5g',H_inf_S));
disp(sprintf('H_inf norm T: %0.5g',H_inf_T));

