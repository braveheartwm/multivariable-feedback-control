
% function padetf = padetf(theta,n)
%
% Pade approximation of order n of delay theta (actually this is 
% n first-order Pade approximations in series). 
%
% This gives a transfer function approximation.
% To compute the exact frequency response of a scalar delay use delay.m
%
% S. Skogestad ; Jan. 1996 / May 2007

function padetf = padetf(theta,n)

if (nargin == 0) | (nargin > 3),
   disp('usage: padetf = padetf(theta,n)')
   return
end

padetf=1;
delay = tf([-theta/(2*n) 1],[theta/(2*n) 1]);
for i=1:n padetf=padetf*delay; 
end


