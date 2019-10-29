% Example 4.13, p.138 - Compute pole and zero directions
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Expl4_13.m,v 1.5 2004/04/02 17:34:41 vidaral Exp $
%
%
% G(s) = 1/(s+2) [ s-1    4    ]      (eq. 4.68)
%                [ 4.5  2(s-1) ]

clear all; close all;
G=tf({[1 -1], [4]; [4.5], [2 -2]}, {[1 2], [1 2]; [1 2], [1 2]});
Gr=minreal(ss(G)); % Minimum realization

%---------------------------------------------------------------------
% First compute pole and zero directions the "crude" way using SVD
% This method is NOT generallay recommended, e.g., with several poles
% and zeros at the same location.
%---------------------------------------------------------------------
p=pole(Gr)        % one stable   (LHP) pole at p = -2
z=zero(Gr)        % one unstable (RHP) zero at z = 4
pp = p*(1+1.e-4); % move the pole a little away before evaluating G(p)
Gp=evalfr(Gr,pp);
[u,s,v]=svd(Gp);  % Take SVD of G(p)
yp = u(:,1)       % Output pole direction is [-0.5547  0.8321]'
up = v(:,1)       % Input  pole direction is [-0.6     0.8   ]'

Gg=evalfr(Gr,-i*4.76434); r=rank(Gg) % r is Normal rank of G
Gz = evalfr(G,-i*z); 
Gz = evalfr(G,z); 
[u,s,v]=svd(Gz);
yz = u(:,r)  % Output zero direction is [-0.8321 0.5547]'
uz = v(:,r)  % Input  zero direction is [-0.8    0.6   ]'

%---------------------------------------------------------------------
% Then compute the poles zeros from the state-space description.
% This is the PREFERRED way.
%---------------------------------------------------------------------

% Pole directions first; see (4.37) and (4.42)
[A,B,C,D]=ssdata(Gr);
eig(A)                       % the poles ...
[T,P] = eig(A);
YP = C*T                     % associated pole output vectors
YPn = YP / norm(YP)          %   normalize to get pole output directions
eig(A')                      % the poles again ...
[Q,P] = eig(A'); UP = B'*Q   % associated pole input vectors
UPn= UP / norm(UP)           %   normalize to get pole input directions

% Now the zero directions, see (4.70) to (4.71)
% NOTE: Functions ozde and izde are available on the book's web-page
[Zy,YZ,Xy] = ozde(Gr)        % ozde: function for zero output directions
[Zu,UZ,Xu] = izde(Gr)        % izde: function for zero input directions


