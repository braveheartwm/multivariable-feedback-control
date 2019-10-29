%   [Po, Yp, Xpo, Spo] = opde(G,epp,RF)
%
%   The ``Output-Pole-Direction-Eigenvalue'' opde-function computes the 
%        output pole directions through the use of  egenvalue decomposition.
%   NOTE:  Avoids the change in ordering of eigenvalues that may occur if 
%        the input and output directions are computed separately.
%
%   Outputs: Po   - Vector containing the poles.
%            Yp   - Matrix containing the output pole vectors/directions (normalized = length 1), 
%                  column i corresponds to pole Po(i).
%            Xpo   - Output state directions  = right eigenvectors of A (non-normalized).
%                    NOTE: These are scaled such that pole vectors = pole directions
%            Spo   - Scalars used to scale the state directions (eigenvectors).
%                 
%   Inputs:  G   - System matrix on mutools format: G=pck(A,B,C,D)
%            epp - (optional) Tolerance, default value is 1e-12.
%                  It is used to check the norm of the output directions.
%                  If the norm of output direction is less than EPP,
%                  zero is stored in the corresponding output direction.
%            RF  - (optional) If FLG_RF==1 then the first element of 
%                  the output direction is made real.
%
%   Written by: Kjetil Havre, 12/9-1996, e-mail: kjetil@ife.no
%   Requires subroutines: zerm.m and sorte.m
%
%   See also:  IPDE, PDFST, IZDE, OZDE, ZDSVD and PDSVD.
%
%   Reference: Havre K. and S. Skogestad, 1996, ``Effect of RHP Zeros and
%              Poles on Performance in Multivariable Systems''.
%
   function [P,Y,X,S] = opde(G,epp,rf);

   P = []; Y = []; X = []; 
   [mt,ny,nu, nx] = minfo(G);
   S = ones(nx,1);
   
   if strcmp(mt, 'syst')==0
      disp( 'System matrix G is required, usage: [] = opde(G,epp)' );
      return
   end

   if nargin < 3
     rf = 0;
   end
   if nargin < 2
      epp = 1e-12;
   end
   [A, B, C, D] = unpck(G);
   
   if norm(D) > 1/epp
      disp('Warning: Directions may be inaccurate due to large effect from D.');
   end
   
   [V, D] = eig(A,'nobalance');   P = diag(D);
   [P,Is] = sorte(P); V = V(:,Is);
   Hm = A*V-V*diag(P);
   if norm(Hm) > 100*epp
      disp('Warning: inaccurate eigenvalue computation.')
      Hm = Hm
   end
   

   Npd = max(size(P));
   for i=1:Npd	 
      X(:,i) = V(:,i);
      if abs(P(i)-real(P(i))) < epp 
         P(i) = real(P(i));
	 if isreal(X(:,i)) == 0 % Rotate eigenvector.
	    inz = find(abs(X(:,i))>epp); 
	    if isempty(inz) == 0     % Rotate eigenvector.
	       vh = angle(X(inz(1),i));
               X(:,i)=X(:,i)*exp(-sqrt(-1)*vh);
	    end
            X(:,i) = zerm(X(:,i),epp/10);
         end
      end
      
      yph = C*X(:,i); nyph = norm(yph);
      
      if nyph > epp
         Y(:,i) = yph/nyph;
	 X(:,i) = X(:,i)/nyph;	 
	 S(i,1) = 1/nyph;
	 if (isreal(Y(:,i))==0) & (rf==1)   % Rotate output direction. 
            inz = find(abs(Y(:,i))>epp); 
	    if isempty(inz) == 0
	       vh = angle(Y(inz(1),i));
               Y(:,i) = Y(:,i)*exp(-sqrt(-1)*vh); 
               X(:,i)=X(:,i)*exp(-sqrt(-1)*vh);
	       S(i,1) = S(i,1)*exp(-sqrt(-1)*vh);
            end
         end
      else
         Y(:,i) = zeros(ny,1);    % Norm is small unobservable.
	 S(i,1) = 0;
      end
   end
   P = zerm(P,epp/10); X = zerm(X,epp/10); Y = zerm(Y,epp/10);
end

   

