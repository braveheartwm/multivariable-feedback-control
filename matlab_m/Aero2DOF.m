
%  The first part is the aero-engine (jusr for testing)
clear
load aero1
G = G5; 
				
% build weight Wpbar=I/s (where I is a 3*3 identity matrix)
sys1=zpk([],0,1);
Wpbar=blkdiag(sys1,sys1,sys1);  
% WEIGHT CROSSOVER EQUAL TO 1 rad/second.
% build weight W2=I
W2=eye(3);
% build W2*G*Wpbar and find the align gain Wa (at 7 rad/sec)
% W2GWpbar=mmult(W2,G,Wpbar);
W2GWpbar=W2*G*Wpbar;
[a,b,c,d]=ssdata(W2GWpbar);
ff=c*inv(sqrt(-1)*7*eye(size(a))-a)*b+d;	
Wa=align(ff);
% build Wg, W1 and the final shaped plant W2*G*W1
Wg=blkdiag(1,2.5,.3);
W1=Wpbar*Wa*Wg;

Gs=W2*G*W1;
[As,Bs,Cs,Ds] = ssdata(Gs);

% build the reference model Tref
sys1=tf(1,[.018 1]);sys2=tf(1,[.008 1]);sys3=tf(1,[.2 1]);
Tref=blkdiag(sys1,sys2,sys3);
[Ar,Br,Cr,Dr] = ssdata(Tref);

K=hinf2dof(Gs,Tref);
% but get extra states, since states from Gs come twice
%
% Final gamma value is: 2.9091
%


