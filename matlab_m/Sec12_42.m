% Example of an LMI fesibility problem (Lyapunov equation)
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Lyap_LMI.m,v 1 2004/26/04 kariwala Exp $

A = randn(4,4);

setlmis([])
P = lmivar(1,[size(A,1) 1]);   % Specfiy the matrix variable

Lyap = newlmi                  
lmiterm([Lyap 1 1 P],1,A,'s')  % Specifying the LMI AP + P'A < 0
lmiterm([Lyap 1 2 0],0)
lmiterm([Lyap 2 2 P],-1,1)     % P > 0

LMIsys = getlmis;               % Obtaining the system of LMIs

[tmin,xfeas] = feasp(LMIsys);  % Solving the feasibility problem
% Feasible (A is stable) iff tmin < 0 

Pfeas = dec2mat(LMIsys,xfeas,P); % Obtaining the feasible P

evalsys = evallmi(LMIsys,xfeas); % Validation of results
[Lhs,Rhs]=showlmi(evalsys,1); % Feasible if Lhs < 0 


