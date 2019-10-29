
% function fod = fod(k,tau,theta,n)
%
% Generates state-space realization of first-order with delay process
%  k     - gain
%  tau   - time constant
%  theta - delay
%  n     - order used for Pade approximation of delay
%
% This is the chemical engineers favorite process...
% S. Skogestad ; Mar. 1997 / corrected May 2007 (divide by 2n and not n)

function  fod = fod(k,tau,theta,n)

if (nargin == 0) | (nargin > 5),
   disp('usage:  fod = fod(k,tau,theta,n)')
   return
end

pade=1;
delay = nd2sys([-theta/(2*n) 1],[theta/(2*n) 1]);
for i=1:n pade=mmult(pade,delay); 
end

fod = mmult( nd2sys(k,[tau 1]), pade); 


