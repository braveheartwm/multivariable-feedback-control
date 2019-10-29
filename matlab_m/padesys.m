
% function padesys = padesys(theta,n)
%
% Pade approximation of order n of delay theta (actually this is 
% n first-order Pade approximations in series). 
%
% This gives a transfer function approximation.
% To compute the exact frequency response of a scalar delay use delay.m
%
% S. Skogestad ; Jan. 1996 / corrected May 2007 and renamed from pade to padesys to avoid confilit lti/pade.m

function padesys = padesys(theta,n)

if (nargin == 0) | (nargin > 3),
   disp('usage: padesys = padesys(theta,n)')
   return
end

padesys=1;
delay = nd2sys([-theta/(2*n) 1],[theta/(2*n) 1]);
for i=1:n padesys=mmult(padesys,delay); 
end


