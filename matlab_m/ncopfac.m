function [N,M]=ncopfac(a,b,c,d)
%NCOPFAC    [M,N]=ncopfac(A,B,C,D)
% Find the normalized coprime factors of system [A,B,C,D] using (4.26).
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: ncopfac.m,v 1.3 2004/04/02 17:36:26 vidaral Exp $
           
S=eye(size(d'*d))+d'*d; 
R=eye(size(d*d'))+d*d';
A1 = a-b*inv(S)*d'*c;
R1 = c'*inv(R)*c;
[R1s,R1err]=sqrtm(R1);
Q1 = b*inv(S)*b';
[Z,L,G]=care(A1',R1s,Q1);  %Solve Ricatti eq

H = -(b*d' + Z*c')*inv(R);
A = a + H*c;
Bn = b + H*d;
Bm = H;
C = inv(sqrtm(R))*c;
Dn = inv(sqrtm(R))*d;
Dm = inv(sqrtm(R));
N = ss(A,Bn,C,Dn);
M = ss(A,Bm,C,Dm);
