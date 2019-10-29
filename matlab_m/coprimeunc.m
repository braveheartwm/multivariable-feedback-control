%function [Ac,Bc,Cc,Dc,gammin] = coprimeunc(a,b,c,d,gamrel)
%
% Finds the controller which optimally robustifies a given
% shaped plant in terms of tolerating maximum coprime uncertainty.
% Used in the McFarlane-Glover H-infinity loopshaping procedure.
% --- Uses the robust control toolbox ---
%
%    a,b,c,d:     State-space description of (shaped) plant  
%    gamrel:      Final gamma used is gamrel*gammin [default: gamrel=1.1]
%    Ac,Bc,Cc,Dc: State-space description of "robustifying" controller
%                 assuming positive feedback
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite

% $Id: coprimeunc.m,v 1.2 2004/04/15 08:10:13 vidaral Exp $

function [Ac,Bc,Cc,Dc,gammin] = coprimeunc(a,b,c,d,gamrel)
if nargin <4, 
   disp('usage:  [Ac,Bc,Cc,Dc] = coprimeunc(a,b,c,d,gamrel)'); return; end
if nargin <5, gamrel=1.1; end

% Find Normalized Coprime factors of the shaped plant
S=eye(size(d'*d))+d'*d;
R=eye(size(d*d'))+d*d';
Rinv=inv(R);Sinv=inv(S);

A1 = (a-b*Sinv*d'*c); R1 =S; B1=b; Q1 = c'*Rinv*c;
[X,XAMP,G,REP]=care(A1,B1,Q1,R1);
if REP == -1
    disp('The Hamiltonian matrix has jw-axis eigenvalues')
elseif REP == -2
    disp('There is no finite stabilizing solution X')
else
    sprintf('X: Frobenius norm of relative residual= %0.5g',REP)
end    

A2 = A1'; Q2 = b*Sinv*b'; B2=c'; R2 = R;
[Z,ZAMP,G,REP]=care(A2,B2,Q2,R2);
if REP == -1
    disp('The Hamiltonian matrix has jw-axis eigenvalues')
elseif REP == -2
    disp('There is no finite stabilizing solution X')
else
    sprintf('Z: Frobenius norm of relative residual= %0.5g',REP);
end   

% display optimal gamma
XZ = X*Z; gammin=sqrt(1+max(eig(XZ)))

% Use higher gamma 
gam=gamrel*gammin; gam2 = gam*gam; gamconst = (1-gam2)*eye(size(XZ)); 
Lc = gamconst + XZ; Li = inv(Lc'); Fc = -Sinv*(d'*c+b'*X);
Ac = a + b*Fc + gam2*Li*Z*c'*(c+d*Fc);
Bc = gam2*Li*Z*c';
Cc = b'*X; 
Dc = -d';
%---------------------------------------------------------------------
