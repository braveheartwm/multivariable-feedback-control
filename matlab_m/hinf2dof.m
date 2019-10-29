function K=hinf2dof(Gs, Tref)
%HINF2DOF  synthesizes the H_inf 2-DOF controller in (9.80) 
% Uses MATLAB RCAST toolbox
% Usage: K=hinf2dof(Gs, Tref);
% INPUTS:  Shaped plant Gs and reference model Tref
% OUTPUT: Two degrees-of-freedom controller K
% 
% Coprime factorization of Gs 
[As,Bs,Cs,Ds] = ssdata(Gs);
[Ar,Br,Cr,Dr] = ssdata(Tref);
[nr,nr] = size(Ar); [lr,mr] = size(Dr);
[ns,ns] = size(As); [ls,ms] = size(Ds);
Rs = eye(ls)+Ds*Ds.'; Ss = eye(ms)+Ds'*Ds;
A1 = (As - Bs*inv(Ss)*Ds'*Cs); 
R1 = Cs'*inv(Rs)*Cs; Q1 = Bs*inv(Ss)*Bs'; 
% [Z1, Z2, fail, reig_min] = ric_schr([A1' -R1; -Q1 -A1]); Zs = Z2/Z1;
% Alt. Robust Control toolbox: 
[Z1,Z2,eig,zerr,wellposed,Zs] = aresolv(A1',Q1,R1);
%
% Choose rho=1 (Designer's choice) and
% build the generalized plant P in (9.87)
% 
rho=1;
A = blkdiag(As,Ar);
B1 = [zeros(ns,mr) ((Bs*Ds')+(Zs*Cs'))*inv(sqrtm(Rs));
      Br zeros(nr,ls)];
B2 = [Bs;zeros(nr,ms)];
C1 = [zeros(ms,ns+nr);Cs zeros(ls,nr);rho*Cs -rho*rho*Cr];
C2 = [zeros(mr,ns+nr);Cs zeros(ls,nr)];
D11 = [zeros(ms,mr+ls);zeros(ls,mr) sqrtm(Rs);-rho*rho*Dr rho*sqrtm(Rs)];
D12 = [eye(ms);Ds;rho*Ds];
D21 = [rho*eye(mr) zeros(mr,ls);zeros(ls,mr) sqrtm(Rs)];
D22 = [zeros(mr,ms);Ds];
B = [B1 B2]; C = [C1;C2]; D = [D11 D12;D21 D22];
P = ss(A,B,C,D);
% Alternative: Use sysic to generate P from Figure 9.21 
% but may get extra states, since states from Gs may enter twice. 
%
% Gamma iterations to obtain H-infinity controller
%
[l1,m2] = size(D12); [l2,m1] = size(D21); 
nmeas = l2; ncon = m2; gmin = 1; gmax = 5; gtol = 0.01;
[K, Gnclp, gam] = hinfsyn(P, nmeas, ncon, gmin, gmax, gtol);
% Alt. Robust toolbox, use command: hinfopt
%
% [K, Gnclp, gam] = hinfopt(P, nmeas, ncon, gtol, gmax, gmin)
