function [ hankel_sys, error_bound ] = modred_hankel( sys, elim )
% modred_hankel: Calculates the optimal Hankel-norm approximation
%
% This function, given a system in state-space form, and the number
% of the first of the Hankel singular values to be discarded with,
% calculates the optimal Hankel-norm approximation of the system, 
% as shown in Theorem 11.2 in the book. The data has to be entered
% as follows:
% 
%       sys: system to be reduced, as an ss variable;
%       elim: index of first Hankel singular value to be discarded.
%
% The resulting system will be returned in state-space form:
%
%       hankel_sys: Hankel-optimal reduced system
%       error_bound: the error bound according to formulae 11.29 and 11.30
%
% See pp. 530-533 in the book for the theory.
%
% ***Remark: for use only with square and stable transfer functions.
% 
% Copyright 1996-2004 Sigurd Skogestad & Ian Postlethwaite
% $Id: modred_hankel.m,v 1.6 2004/02/23 13:10:57 zenith Exp $


% Checking validity of input

if class( sys ) ~= 'ss'
    error( 'The function did not receive a system variable as first input' ) ;
end

if ( elim < 1 ) | ( round( elim ) ~= elim )
    error( 'Only positive integer values are acceptable for the first state to be discarded')
end

if size( sys.a, 1 ) < elim
    error( 'The original system has less than %i states', elim ) ;
end

% Balancing the state-space realization

[ sys, SIG ] = balreal(sys) ;

% If elim is provided as a vector, take the smallest value

elim = min( elim ) ;


% Check if, and if so how many times, the "elim-th" value is repeated

l = 1 ;
sigmaelim = SIG( elim ) ;   % First Hankel singular value to be discarded

for i = ( elim + 1 ):max( size( SIG ) )
    if sigmaelim == SIG( i )
        l = l + 1 ;
    else
        break
    end
end


% Index sets for Hankel singular values different from and equal to
% sigmaelim; most of the times indexeq has only one element

indexdiff = [ 1:(elim-1), ( elim + l ):max( size( SIG ) ) ] ;
indexeq   = [ elim:( elim + l -1 ) ] ;


% Rearrange Hankel singular values, build SIGMA1 matrix

SIG1 = diag( SIG( indexdiff ) ) ;


% Partition system matrices

A11 = sys.a( indexdiff, indexdiff ) ;
B1  = sys.b( indexdiff, : ) ;
B2  = sys.b( indexeq, : ) ;
C1  = sys.c( :, indexdiff ) ;
C2  = sys.c( :, indexeq ) ;

% Precalculate other useful matrices

U        = - C2 * B2 / ( B2 * B2' ) ;
GAMMA    = SIG1^2 - sigmaelim^2 * eye( size( SIG1 ) ) ;
INVGAMMA = inv( GAMMA ) ;

% Calculate "chapeau" (hat) state-space model

A_chapeau = INVGAMMA * ( sigmaelim^2 * A11' + SIG1 * A11 * SIG1 - sigmaelim * C1' * U * B1' ) ;
B_chapeau = INVGAMMA * ( SIG1 * B1 + sigmaelim * C1' * U ) ;
C_chapeau = C1 * SIG1 + sigmaelim * U * B1';
D_chapeau = sys.d - sigmaelim * U ;


% Extract indeces of stable and unstable eigenvalues

[ V, Eigen ] = eig( A_chapeau ) ;

stab   = vec2ind( real( diag( Eigen ) ) <  0 ) ;
unstab = vec2ind( real( diag( Eigen ) ) >= 0 ) ;

% Switch to block-diagonal form, removing imaginary parts

[ V, Eigen ] = cdf2rdf( V, Eigen ) ;

B = inv( V ) * B_chapeau ;
C = C_chapeau * V ;

% Create stable subsystem

Ag = Eigen( stab, stab ) ;
Bg = B( stab, : ) ;
Cg = C( :, stab ) ;
Dg = D_chapeau ;

% Return reduced system and error bound

hankel_sys = ss( Ag, Bg, Cg, Dg ) ;

error_bound = sigmaelim + sum( SIG( ( elim + l ):end ) ) ;