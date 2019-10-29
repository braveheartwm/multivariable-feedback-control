
% function delay_g = delay(theta,omega,nu,epp)
%
% This computes the EXACT frequency response of a scalar delay 
% To find a Pade transfer function approximation use pade.m
%
% delay_g = delay(theta,omega,nu,epp)
% delay_g : frequency respones of  a delay
%           delay_g = exp(-s*theta)
% theta : time delay
% omega : frequency vector
% nu    : dimension of varying matrix     (default = 1)
% epp   : small negative A-matrix element (default = 0)

% Author: Petter Lundstrom
%

function delay_g = delay(theta,omega,nu,epp)

if (nargin == 0) | (nargin > 4),
   disp('usage: delay_g = delay(theta,omega,nu,epp)')
   return
end
if nargin < 4,
   epp = 0;
end
if nargin < 3,
   nu = 1;
end

if epp > 0,
   epp = -epp;
end

i = pck(epp,1,1,0);
tmp = frsp(i,omega,theta);
delay_g = tmp;
for j=2:nu,
   delay_g = daug(delay_g,tmp);
end

