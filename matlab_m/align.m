function A=align(V)
% A=ALIGN(V) returns a constant matrix A which is the real alignment
%       of the INVERSE of the complex input matrix V, i.e. 
%
%            A=arg min_{A,theta}||V*A-diag(exp(j*theta)||
%
% References:
% 1. Maciejowski, J.M., Multivariable Feedback Design, 1989, pp. 145--148
% 2. Edmunds, J., and Kouvaritakis, B., Extensions of the frame alignment
%   technique and their used in the characteristic locus design method.
%   Internation J. Control, 1979, 29, (5), pp.787--796
%
% By Yi Cao, 1 May 1996, University of Exeter
%
% Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
%
% $Id: align.m,v 1.1 2004/01/22 18:48:49 zenith Exp $

if ( nargin == 0 ) | ( nargin > 1 ),
   disp( 'usage: mat_inv_real = align(mat)' )
   return
end

D = pinv( real( V' * V )  ) ;
A = D * real( V' * diag( exp( j * angle( diag( V * D *  V.' ) ) / 2 ) ) ) ;