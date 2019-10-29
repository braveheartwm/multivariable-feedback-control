function RGA = vrga(Gf);
%VRGA    RGA=vrga(G) returns a frequency-varying matrix (in mu-toolbox)
% 	 of the RGA of an input frequency-varying matrix G.
%
%   $Id: vrga.m,v 1.1 2004/01/16 14:05:46 standber Exp $

if (nargin == 0) | (nargin > 1),
   disp('usage: vrga(mat)')
   return
end
RGA = feval('.*',Gf,pinv(transpose(Gf)));
%Alternative:
%RGA=veval('rga',Gf);


