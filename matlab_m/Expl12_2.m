% Example of an Linear objective minimization problem (H-infinity norm)
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
% $Id: Hinf_LMI.m,v 1 2004/26/04 kariwala Exp $
clear all;
clc;

G = rss(5,2,2);
[A,B,C,D] = ssdata(G);

setlmis([])
P = lmivar(1,[size(A,1) 1]);   % Specfiy the matrix variables
gamma = lmivar(1,[1 1]);

HinfLMI = newlmi                  % New LMI 

% Only terms above the diagonal need to be sepcified
lmiterm([HinfLMI 1 1 P],1,A,'s')  % PA + A'P 
lmiterm([HinfLMI 1 2 P],1,B)      % PB
lmiterm([HinfLMI 1 3 0],C')       % C'
lmiterm([HinfLMI 2 2 gamma],-1,1) %-gamma.I
lmiterm([HinfLMI 2 3 0],D')       %D'
lmiterm([HinfLMI 3 3 gamma],-1,1) %-gamma.I

Ppos = newlmi                  % New LMI 
lmiterm([Ppos 1 1 P],-1,1)     % P > 0

LMIsys = getlmis;               % Obtaining the system of LMIs

c = mat2dec(LMIsys,zeros(size(A,1),size(A,1)),1);

options = [1e-5,0,0,0,0]; % Relative accuracy of solution

[copt,xopt] = mincx(LMIsys,c, options);  % Solving the minimization problem
 
Popt = dec2mat(LMIsys,xopt,P);  %Obtaining the optimal P


