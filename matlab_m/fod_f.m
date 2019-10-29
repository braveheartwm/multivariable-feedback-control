
% function fod_f = fod_f(k,tau,theta,omega)
%
% Generates EXACT frequency response first-order with delay process
%  k     - gain
%  tau   - time constant
%  theta - delay
%  omega - frequency vector
%
% This is the chemical engineers favorite process...
% S. Skogestad ; Mar. 1997

function  fod_f = fod_f(k,tau,theta,omega)

if (nargin == 0) | (nargin > 5),
   disp('usage:  fod_f = fod_f(k,tau,theta,omega)')
   return
end

epp=0;
i = pck(epp,1,1,0);
delay_f = frsp(i,omega,theta);

time_f = frsp( nd2sys(k,[tau 1]), omega);
fod_f = mmult( time_f, delay_f); 




