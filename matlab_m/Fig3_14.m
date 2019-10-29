%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section 3.7.2 DISTILLATION PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%   The example process is from:
%   Skogestad, Morari and Doyle,``Robust Control of Ill-Conditioned
%   Plants: High-Purity Distillation'', IEEE Atomat. Control 33,
%   1092-1105 (1988).
%
% Copyright 2005 Sigurd Skogestad & Ian Postlethwaite
% Multivariable Feedback Control 2nd Edition
% $Id: Sec3_7_2.m,v 1.3 2004/02/09 12:01:24 vidaral Exp $

%Plant (3.82): distilation column, G(s)=1/(75s+1)[87 -86.4;108.2 -109.6]
G0 = [87.8 -86.4; 108.2 -109.6];
[U,sval,V]=svd(G0);
dyn = tf(1,[75 1]); 
G=dyn*G0;

Lambda=rga(G0);
disp(sprintf('\nRGA(0)=[%6.1f%6.1f]\n       [%6.1f%6.1f]',Lambda'));

% Inverse-based controller (3.84)
dynk = 0.7*tf([75 1],[1 1.e-6]); 
Kinv = (dynk*inv(G0));
I2=eye(2); 

% Nominal
S=minreal(inv(I2+G*Kinv));
T=minreal(I2-S);

% With 20% INPUT uncertainty
Unc = [1.2 0; 0 0.8];
Lu=G*Unc*Kinv;
Su=minreal(inv(I2+Lu));
Tu=minreal(I2-Su);

[Sp,wSp]=norm(S,inf); [Sup,wSup]=norm(Su,inf); % Peak of S; 1 nominally and 14.2 with uncertainty
[Tp,wTp]=norm(T,inf); [Tup,wTup]=norm(Tu,inf); % Peak of T; 1 nominally and 14.2 with uncertainty
disp(sprintf('						 %6.5s %6.5s', '(S/T)', 'w'));
disp(sprintf('						 %12.121212121212121212121212s','----------------'));
disp(sprintf('Peak of singular value S and w nominally	: %6.5g %6.5g',Sp,wSp));
disp(sprintf('Peak of singular value T and w nominally	: %6.5g %6.5g',Tp,wTp));
disp(sprintf('Peak of singular value S and w robust \t	: %6.5g %6.5g',Sup,wSup));
disp(sprintf('Peak of singular value T and w robust \t	: %6.5g %6.5g',Tup,wTup));

% TIME  simulation
Kr=tf(1,[5 1]); % 5 min filter on reference change
Tr=T*Kr;

u1=[1*ones(1001,1) 0*ones(1001,1)];
t=[0:0.1:100];

y=lsim(Tr,u1,t);
u=lsim(Kinv*S*Kr,u1,t);

Tru=Tu*Kr;
yu=lsim(Tru,u1,t);
uu=lsim(Kinv*Su*Kr,u1,t);

% Plotting
subplot(211); p=plot(t,y(:,1),'-',t,y(:,2),'-',t,yu(:,1),'--',t,yu(:,2),'--');  %Figure 3.12
text(37,2.3,'Nominal plant: Solid Line ');text(37,2,'Perturbed plant: Dashed Line');
axis([0 60 -0.2 2.7]);text(30,1.2,'y1'),text(30,0.2,'y2');title('OUTPUTS');
subplot(212);plot(t,u,'-',t,uu,'--'),title('INPUTS');figure(1);xlabel('TIME (min)');
