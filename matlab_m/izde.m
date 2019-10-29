%   [Z, U, X] = izde(G,EPP);
%     
%   Inputs:  G   - system matrix in mu-tools format.
%            EPP - tolerance, see below, default value EPS.
%   Outputs: Z   - zeros.
%            U   - input zero directions, stored as column vectors.
%            X   - state directions, stored as column vectors.
%            Each column X(:,i) and U(:,i) corresponds to the zeros Z(i).
%
%   This is a modifaction of szeros.m written by Kjetil Havre
%
%   This Input-Zero-Direction-"through genaralized Eigenvalue" decomposition
%   IZDE function is a modification of the szeros function contained in mu 
%   - toolbox. The modification consists of returning the input zero 
%   directions and the state zero directions in addition to the zeros. The
%   input zero directions u are defined as:
%              G(z)*u = 0,
%   where s = z is a zero of G(s). This is done by solving the generalized 
%   eigenvalue problem:
%
%           | A-Iz |  B  | * | x |  =  | 0 |
%           |------------|   |---|     |---|
%           |  C   |  D  |   | u |     | 0 |
%
%   IZDE finds the transmission zeros z of a SYSTEM matrix. Occasionally, 
%   large zeros are included which may actually be at infinity. Solving
%   for the transmission zeros of a system involves two generalized eigen-
%   value problems. EPP (optional) defines if the difference between two
%   generalize eigenvalues is small. IZDE also finds the input u and the 
%   state x directions of the zeros.
%
%   The input zero directions are stored as column vectors in U, and each  
%   of the columns are normalized. The state zero directions are stored as 
%   columns in X. The degree of freedom to normalize the generalized eigen-
%   vector is used to normalize the u part. So the length of x is not equal 
%   to one. Each column in U and X corresponds to the element in Z with 
%   same place.
%
%   For systems with more inputs than outputs the input zero direction
%   is not a complete basis for the nullspace of G(z). 
%   Zeros with multiplicity greater than one (rare cases which may 
%   occure in non-minimal realizations), may (not sure) cause wrong 
%   directions.
%
%   Comments, corrections and malfunctions, can be e-mailed to:
%          havre@kjemi.unit.no or skoge@kjemi.unit.no
%
%   See also: EIG, SZEROS, OZDE and SPOLES. 
%   Algorithm based on Laub & Moore 1978 paper, Automatica
%
   
%   Note that when the number of inputs is larger than the number of outputs,
%   the input zero direction is not complete. A bit clearer: If z is a zero of 
%   a non-square plant with number of inputs greater than number of outputs,
%   the zero direction is not a line but a surface. As an example, consider 
%   G(s) with dimensions 2x3 (2 rows and 3 columns). Let s=z be a zero of 
%   G(s) such that G(z)*u = [0;0]; Since s=z is a zero then the rank of
%   G(z) has to be less than the normal rank of G(s), which at maximum can 
%   be 2. This implies that rank of G(z) must be less than 2.
%   u is element in the three dimensional field of real numbers. Since the
%   rank of G(z) is maximum one the zero direction is a actually a subspace
%   in this three dimensional field of real numbers given by two basis 
%   vectors. Since this function only gives one zero direction for a given 
%   zero this direction does not describe the input zero space completely.
%   Two basis vectors are requiered. 
%
%   Modification for square systems was made by:      Kjetil Havre 14/5-1995.
%   Modification for non square systems was made by:  Kjetil Havre 14/5-1995.
%   Inclusion of state directions was made by:        Kjetil Havre 14/5-1995.
%   Modified so that first element of U(:,i) is real: Kjetil Havre  3/2-1996. 
%   Copyright 1996-2003 Sigurd Skogestad & Ian Postlethwaite
%   $Id: izde.m,v 1.2 2004/01/19 14:52:11 aske Exp $

function [Z, U, X] = izde(sys,epp)
  if nargin < 1
    disp('usage:  [Z,Y,X] = ozde(G) ')
    return
  end
  if nargin == 1
    epp = eps;
  end
  [ny,nu,nx]=size(sys);
  if class(sys) == 'tf'|'ss'|'zpk'|'frd'
    [a,b,c,d] = ssdata(sys);
    if nx == 0
      disp('SYSTEM has no states')
    end
    sysu = [a b; c d];
    
% find generalized eigenvalues of a square system matrix

    if ny == nu
      x = zeros(nx+nu,nx+nu);
      x(1:nx,1:nx) = eye(nx);
      [vech, ev] = eig(sysu, x);

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
      
% Split x and u.
      vx = vech2(1:nx, : );
      vu = vech2(nx+1:nx+nu,:);
% Normalize columns.
      [nvr, nvc] = size( vu ); 
      for i=1:nvc,
          nrmu = norm(vu(:,i));
	  if nrmu > 1000*epp
	     vx(:,i) = vx(:,i)/nrmu;
	     vu(:,i) = vu(:,i)/nrmu;
	  else
	     vu(:,i) = zeros(nu,1);
	  end
	  Inz = find( abs(vu(:,i)) > 1000*epp );
	  if isempty(Inz) == 0
             U(:,i) = vu(:,i) * exp( -angle(vu(Inz(1),i))*sqrt(-1) );
	     X(:,i) = vx(:,i) * exp( -angle(vu(Inz(1),i))*sqrt(-1) );
	  else
	     X(:,i) = vx(:,i);
	     Y(:,i) = vu(:,i);
	  end
      end

    else                                  % Non-square systems 
      nrm = norm(sysu,1);
      if nu < ny
        x1 = [ sysu  nrm*(rand(nx+ny,ny-nu)-.5)];
        x2 = [ sysu  nrm*(rand(nx+ny,ny-nu)-.5)];
      else
        x1 = [ sysu; nrm*(rand(nu-ny,nx+nu)-.5)];
        x2 = [ sysu; nrm*(rand(nu-ny,nx+nu)-.5)];
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
      
% Split in ux and xz
      if isempty( vech3 ) 
         Z = []; U = []; X = [];
	 return;
      end
      vx = vech3(1:nx, : );
      vu = vech3(nx+1:nx+nu,:);
% Normalize columns.
      [nvr, nvc] = size( vu ); 
      for i=1:nvc,
          nrmu = norm(vu(:,i)); 
	  if nrmu > 1000*epp
	     vx(:,i) = vx(:,i)/nrmu;
	     vu(:,i) = vu(:,i)/nrmu;
	  else
	     vu(:,i) = zeros(nu,1);
	  end
	  Inz = find( abs(vu(:,i)) > 1000*epp );
	  if isempty(Inz) == 0
             U(:,i) = vu(:,i) * exp( -angle(vu(Inz(1),i))*sqrt(-1) );
	     X(:,i) = vx(:,i) * exp( -angle(vu(Inz(1),i))*sqrt(-1) );
	  else
	     X(:,i) = vx(:,i);
	     Y(:,i) = vu(:,i);
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
