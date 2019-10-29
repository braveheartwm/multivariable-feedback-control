% function [Z, Y, X] = ozde(G,epp)
%
%   Inputs:  G   - system matrix in mu-tools format.
%            EPP - tolerance, see below, default value EPS.
%   Outputs: Z   - zeros.
%            Y   - output zero directions, stored as column vectors.
%            X   - state directions, stored as column vectors.
%            Each column X(:,i) and U(:,i) corresponds to the zeros Z(i).
%
%   This is a modifaction of szeros.m written by Kjetil Havre
%
%   This Output-Zero-Direction-"through genaralized Eigenvalue" decomposition
%   OZDE function is a modification of the szeros function contained in mu 
%   - toolbox. The modification consists of returning the output zero 
%   directions and the state zero directions in addition to the zeros.
%   The output zero directions are defined as:
%              y'*G(z) = 0,
%   where s = z is a zero of G(s) and ' denotes conjugate transposed.
%   This is done by solving the generalized eigenvalue problem of the 
%   transposed system:
%
%           | x' | y' | * | A-Iz |  B  |  =  |  0 | 0  |
%                         |------------|
%                         |  C   |  D  |
%
%   are solved by genaralized eigenvalues:
%
%           | A"-Iz |  C"  | * | xi |  =  | 0 |
%           |--------------|   |----|     |---|
%           |  B"   |  D"  |   | yi |     | 0 |
%
%   x = conj(xi); y = conj(yi);
%
%   The prime ' denotes the conjugate transposed and " denotes the ordinary 
%   transposed.
%
%   OZDE finds the transmission zeros z of a SYSTEM matrix. Occasionally, 
%   large zeros are included which may actually be at infinity. Solving for
%   the transmission zeros of a system involves two generalized eigenvalue 
%   problems. EPP (optional) defines if the difference between two generalized
%   eigenvalues is small. OZDE also finds the output y and the state x 
%   directions of the zeros.
%
%   The output zero directions are stored as column vectors in Y, and the y's 
%   are normalized so that Y'*Y = I. The state zero directions are stored as
%   columns in X. The degree of freedom to normalize the generalized eigen-
%   vectors is used to normalize the y part of the vectors. So the length of 
%   x is not equal to one. Each column in Y and X corresponds to the element
%   in Z with same place.
%
%   For systems with more outputs than inputs the output zero direction
%   is not a complete basis for the left nullspace of G(z). 
%   Zeros with multiplicity greater than one (rare cases which may 
%   occure in non-minimal realizations), may (not sure) cause wrong 
%   directions.
%
%   Comments, corrections and malfunctions, can be e-mailed to: 
%          havre@kjemi.unit.no or skoge@kjemi.unit.no
%
%   See also: EIG, SZEROS, IZDE and SPOLES. 
%   Algorithm based on Laub & Moore 1978 paper, Automatica
%

%   Note that when the number of outputs is larger than the number of inputs,
%   the output zero direction is not complete. A bit clearer: If z is a zero 
%   of a non-square plant with number of outputs greater than number of inputs,
%   the zero direction is not a line but a surface. As an example, consider 
%   G(s) with dimensions 3x2 (3 rows and 2 columns). Let s=z be a zero of 
%   G(s) such that Yz'*G(z) = [0 0]; Since s=z is a zero then the rank of
%   G(z) has to be less than the normal rank of G(s), which at maximum can 
%   be 2. This implies that rank of G(z) must be less than 2.
%   y is element in the three dimensional field of real numbers. Since the
%   rank of G(z) is maximum one the zero direction is a actually a subspace
%   in this three dimensional field of real numbers given by two basis 
%   vectors. Since this function only gives one zero direction for a given 
%   zero this direction does not describe the zero space on the output 
%   completely. Two basis vectors are requiered. 
%
%   Modification for square systems was made by:      Kjetil Havre 24/4-1995.
%   Modification for non square systems was made by:  Kjetil Havre 2/5 -1995.
%   Including the state zero dircetions was made by:  Kjetil Havre 15/5-1995.
%   Modified so that first element of U(:,i) is real: Kjetil Havre  3/2-1996. 
%  Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
%  $Id: ozde.m,v 1.2 2004/01/19 14:52:11 aske Exp $

function [Z, Y, X] = ozde(sys,epp)
  if nargin < 1
    disp('usage:  [Z,Y,X] = ozde(G) ')
    return
  end
  if nargin == 1
    epp = eps;
  end
  [ny,nu,nx]=size(sys);
  if class(sys) == 'tf'|'ss'|'zpk'|'frd'
    [a,b,c,d] =ssdata(sys);
    if nx == 0
      disp('SYSTEM has no states')
    end
    sysu = [a b; c d];
    sysu = sysu';
    
% find generalized eigenvalues of a square system matrix

    if ny == nu
      x = zeros(nx+nu,nx+nu);
      x(1:nx,1:nx) = eye(nx);
      [vech, ev] = eig(sysu, x);

      % Add something here to check for not a number or infinite.
      % remove corresponding vectors.
      % Also remove top part corresponding to the states
      % and normalize bottom part. KH 21/4 - 1995.

      z = diag(ev);                       % Extract the eigenvalues.
      kc=0;                               % Counter for eigenvalues.

      for k=1:max(size(z)),
         logic = ~isnan(z(k)) & finite(z(k));
         if logic
            kc= kc+1;
            Z(kc,1) = z(k);
            vech2(:,kc) = vech(:,k);
         end
      end
% Split in x and y.
      vx = vech2(1:nx, : );
      vy = vech2(nx+1:nx+ny,:);
% Normalize columns.
      [nvr, nvc] = size( vy ); 
      for i=1:nvc,
          nrmy   = norm(vy(:,i));
	  if nrmy > 1000*epp
	     vx(:,i) = vx(:,i)/nrmy;
	     vy(:,i) = vy(:,i)/nrmy;
	  else
	     vy(:,i) = zeros(ny,1);
	  end
          Inz = find( abs(vy(:,i) ) > 1000*epp );
	  if isempty(Inz) == 0
             X(:,i) = conj(vx(:,i) * exp( -angle(vy(Inz(1),i))*sqrt(-1) ) );
             Y(:,i) = conj(vy(:,i) * exp( -angle(vy(Inz(1),i))*sqrt(-1) ) );
	  else
	     X(:,i) = conj(vx(:,i));
	     Y(:,i) = conj(vy(:,i));
	  end
      end
      
    else                                  % Non-square systems 
      nrm = norm(sysu,1);
      if nu > ny
        x1 = [ sysu  nrm*(rand(nx+nu,nu-ny)-.5)];
        x2 = [ sysu  nrm*(rand(nx+nu,nu-ny)-.5)];
      else
        x1 = [ sysu; nrm*(rand(ny-nu,nx+ny)-.5)];
        x2 = [ sysu; nrm*(rand(ny-nu,nx+ny)-.5)];
      end
      [x]= zeros(size(x1));
      x(1:nx,1:nx) = eye(nx);
      
      [v1h z1h] = eig(x1,x);            % Compute the genaralized eigenvalues 
      [v2h z2h] = eig(x2,x);            % for the two augumented systems.

      z1h2 = diag( z1h );
      z2h2 = diag( z2h );
      z2 = z2h2(~isnan(z2h2) & finite(z2h2));
      kc=0;                             % Counter for eigenvalues.

      for k=1:max(size(z1h2)),
         logic = ~isnan(z1h2(k)) & finite(z1h2(k));
         if logic
            kc= kc+1;
            z1(kc,1) = z1h2(k);
            vech2(:,kc) = v1h(:,k);
         end
      end
      nz = length(z1);
      vech3 = [];
      Z = [];
      for i=1:nz,
        if any(abs(z1(i)-z2) < nrm*sqrt(epp))
          Z = [Z; z1(i)];
          vech3 = [vech3 vech2(:,i)];
	end
      end      
      if isempty(vech3)
         Z = []; Y = []; X = [];
	 return;
      end
% Split columns in x and y.
      vx = vech3(1:nx,:);
      vy = vech3(nx+1:nx+ny,:);
% Normalize columns.
      [nvr, nvc] = size( vy ); 
      for i=1:nvc,
          nrmy   = norm(vy(:,i));
	  if nrmy > 1000*epp
	     vx(:,i) = vx(:,i)/nrmy;
	     vy(:,i) = vy(:,i)/nrmy;
	  else
	     vy(:,i) = zeros(ny,1);
	  end
          Inz = find( abs(vy(:,i) ) > 1000*epp );
	  if isempty(Inz) == 0
             X(:,i) = conj(vx(:,i) * exp( -angle(vy(Inz(1),i))*sqrt(-1) ) );
             Y(:,i) = conj(vy(:,i) * exp( -angle(vy(Inz(1),i))*sqrt(-1) ) );
	  else
	     X(:,i) = conj(vx(:,i));
	     Y(:,i) = conj(vy(:,i));
	  end
      end
      
    end

  else
     error('matrix is not a SYSTEM matrix')
     return
  end 
%
% Copyright MUSYN INC 1991,  All Rights Reserved
% Copyright MUSYN INC 1995,  All Rights Reserved
%
