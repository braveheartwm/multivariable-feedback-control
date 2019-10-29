%   [Pi, Up, Xpi, Spi] = ipde(G,EPP,RF)
%
%   The ``Input-Pole-Direction-Eigenvalue'' ipde-function computes the 
%        output pole directions through the use of  egenvalue decomposition.
%   NOTE:  Avoids the change in ordering of eigenvalues that may occur if 
%        the input and output directions are computed separately.
%
%   Outputs: Pi   - Vector containing the poles.
%            Up   - Matrix containing the iinput pole vectors/directions (normalized = length 1), 
%                  column i corresponds to pole Pi(i).
%            Xpi   - Input state directions  = left eigenvectors of A (non-normalized).
%                    NOTE: These are scaled such that pole vector = pole directions
%            Spi   - Scalars used to scale the state directions (eigenvectors).
%
%   Inputs:  G   - System matrix on mutools format: G=pck(A,B,C,D)
%            EPP - (optional) Tolerance, default value is 1e-12.
%                  It is used to check the norm of the input directions.
%                  If the norm of input direction is less than EPP,
%                  zero is stored in the corresponding input direction.
%            RF  - (optional) If FLG_RF==1 then the first element of 
%                  the input direction is made real.
%
%   Written by: Kjetil Havre, 12/9-1996, e-mail: kjetil@ife.no
%   Requires subroutines: zerm.m and sorte.m
%
%   See also:  OPDE, PDFST, IZDE, OZDE, ZDSVD and PDSVD.
%
%   Reference: Havre K. and S. Skogestad, 1996, ``Effect of RHP Zeros and
%              Poles on Performance in Multivariable Systems''.
%
   function [P,U,X,S] = ipde(G,epp,rf);

   P = []; U = []; X = []; 
   [mt,ny,nu, nx] = minfo(G);
   S = ones(nx,1);
   
   if strcmp(mt, 'syst')==0
      disp( 'System matrix G is required, usage: [] = ipde(G,epp)' );
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
   
   [V, D] = eig(A.','nobalance');   
   P = diag(D); X = V;
   [P,Is] = sorte(P); X = X(:,Is);
   Hm = A.'*X-X*diag(P);
   if norm(Hm) > 100*epp
      disp('Warning: inaccurate eigenvalue computation.')
      Hm = Hm
      Test = 5
   end
   X = conj(X);
   Hm = X'*A-diag(P)*X';
   if norm(Hm) > 100*epp
      disp('Warning: inaccurate left-eigenvalue computation.')
      Hm = Hm
   end
   
   Npd = max(size(P));
   for i=1:Npd	 
      if abs(P(i)-real(P(i))) < epp 
         P(i) = real(P(i));
	 if isreal(X(:,i)) == 0 % Rotate eigenvector.
	    inz = find(abs(X(:,i))>epp); 
	    if isempty(inz) == 0     
	       vh = angle(X(inz(1),i));
               X(:,i)=X(:,i)*exp(-sqrt(-1)*vh);
	    end
            X(:,i) = zerm(X(:,i),epp/10);
         end
      end
      
      uph = B'*X(:,i); nuph = norm(uph);
      
      if nuph > epp
         U(:,i) = uph/nuph;
	 X(:,i) = X(:,i)/nuph;	 
	 S(i,1) = 1/nuph;
	 if (isreal(U(:,i))==0) & (rf==1)   % Rotate output direction. 
            inz = find(abs(U(:,i))>epp); 
	    if isempty(inz) == 0
	       vh = angle(U(inz(1),i));
               U(:,i) = U(:,i)*exp(-sqrt(-1)*vh); 
               X(:,i)=X(:,i)*exp(-sqrt(-1)*vh);
               S(i,1) = S(i,1)*exp(-sqrt(-1)*vh);
            end
         end
      else
         U(:,i) = zeros(nu,1);    % Norm is small uncontrollable.
	 S(i,1) = 0;
      end
   end
   P = zerm(P,epp/10); X = zerm(X,epp/10); U = zerm(U,epp/10);
end

   













