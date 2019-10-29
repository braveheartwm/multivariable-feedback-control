clear all;
clc
A = [1 -2 4;4 -5 6;7 -8 -9]
setlmis([])
P = lmivar(1,[size(A,1) 1]);
Lyap = newlmi;
lmiterm([Lyap 1 1 P],1,A,'s')
lmiterm([Lyap 1 2 0],1)
lmiterm([Lyap 1 1 P],-1,1)
LMIsys = getlmis;
[tmin,xfeas] = feasp(LMIsys)

% Feasible (A is stable) if tmin < 0

eig(A)

