%   [I, O, XI, XO, err, SI, SO] = pdfst(G, P, SST, NDL)
%
%   The ``Pole-Direction-From-State'' pdfst-function computes the 
%   input and output pole directions through the use of Singular 
%   Value Decomposition of (pI-A) where p is a pole.
%   If smallest singular value of (pI-A) is smaller than SST then 
%   the output pole direction yp can be computed as yp = C*xr,
%   where xr is the eigenvektor of A corresponding to p, A*xr = p*xr. 
%   We note that xr is similar to the input singular vektor 
%   corresponding to the zero singular value.
%   The input pole direction up is computed from the transposed system 
%   i.e.  up = B.'*xl, where xl is the solution to the eigenvalue 
%   problem A^T*xl = p*xl, which corresponds to the output singular 
%   vector with zero singular value of (pI-A).
%   To calculate the zero direction in the state space, the function 
%   uses singular decomposition of (pI-A) and not the the eigen value 
%   formulation. The pole directions are stored as column vectors in 
%   I and O when the smallest singular value is less than SST.
%   If the smallest singular value is larger than SST, the zero vectors
%   are stored in the columns of I and O.
%   A warning message is printed if also the second singular value is 
%   greater than SST. 
%
%   Note! PDFST is numerically more robust than PDSVD.
%   
%   
%   Inputs:  G   - system matrix.
%            P   - vector of poles.
%            SST - tolerance (optional) for singular values, 
%                  default value is 1e-12.
%            NDL - if the norm of the D matrix is large, the directions
%                  computed in PDFST may be wrong. A warning message is
%                  displayed if || D ||_2 >= NDL. Default value is 
%                  1000/SST.
%
%   Outputs: I   - Matrix containing the input directions, 
%                  column i corresponds to pole P(i).
%            O   - Matrix containing the output directions, 
%                  column i corresponds to pole P(i).
%            XI  - Input state directions.
%            XO  - Output state directions.
%            err - Error indicator:
%                     err =  0, all is OK.
%                     err =  1, No state direction with gain smaller than 
%                               SST.
%                     err = -1, One basis vector is not enough to describe 
%                               the subspace of states with gains smaller 
%                               than SST.
%                     err = -2, One basis vector is not enough
%                               to describe the subspace with 
%                               gains larger than EPP.
%            SI  - Scaler used on inputs.
%            SO  - Scaler used on outputs.
%
%   Written by: Kjetil Havre, 30/11-1996, e-mail: kjetil@ife.no
%
%   See also:  IZDE, OZDE, ZDSVD and PDSVD.
%
%   Reference: Havre K. and S. Skogestad, 1995, ``Effect of RHP Zeros and
%              Poles on Performance in Multivariable Systems''.
%
   function [I, O, XI, XO, err, SI, SO] = pdfst(G, P, sst, ndl)

   I = []; O = []; XI = []; XO = []; err = 0;
   [mt,ny,nu, nx] = minfo(G);
   
   if strcmp(mt, 'syst')==0
      disp( 'System matrix G is required, usage: [I, O, XI, XO, err] = pdfst(G, P, SST, NDL)' );
      return
   end
   if nargin < 2
       P = spoles(G);
%      disp( 'Usage: [I, O, XI, XO, err] = pdsvd(G, P, SST, NDL)' );
%      return
   end   
   if nargin < 3
      sst = 1e-12;
   end
   if nargin < 4
      ndl = 1000/sst;
   end
   [A, B, C, D] = unpck(G);
   if norm(D) > ndl
      disp('Warning: Directions may be inaccurate due to large effect from D.');
      err = -2;
   end
   
   Pc = spoles(G);

   Npd = max(size(P));
   
   for i=1:Npd
      Ipc = find( abs(Pc-P(i)) < 10*sst );
      if isempty(Ipc)
         disp(['Warning: Pole: ', num2str(P(i)),' not in system'])
      else
         [U, S, V] = svd( A-eye(nx)*P(i) );
         if S(nx,nx) < 100*sst
             if nx > 1
                 if S(nx-1,nx-1) < sst
                     disp(['Warning: Second smallest singular value is smaller than SST: ', num2str(sst),','])
                     disp(['         S(',int2str(nx-1),',',int2str(nx-1),'): ', num2str(S(nx-1,nx-1)),'.'])
                     err = -1;
                 end
             end
%	 
%            [Vr, Dr] = eig(A);   Lr = diag(Dr)
%            [Vl, Dl] = eig(A.'); Ll = diag(Dl)
%            Ir = find( abs(Lr-P(i)) < sst)
%            Il = find( abs(Ll-P(i)) < sst)
%            XO(:,i) = Vr(:,Ir);
%            XI(:,i) = conj(Vl(:,Il));
%
%  Using eigenvalues does not work since EIG returns wrong eigenvectors 
%  when multiple eigenvalues occures. I have experienced this.
%  Kjetil Havre 20/3 - 1996.
%
            if isreal(P(i))
               XO(:,i) = real(V(:,nx));
               XI(:,i) = real(U(:,nx));
            else
               XO(:,i) = V(:,nx);
               XI(:,i) = U(:,nx);	 
           end
           
           Amp  = A-eye(nx)*P(i);
%           Hpv = A*XO(:,i) - P(i)*XO(:,i)
           xon  = norm(Amp*XO(:,i));
           xin  = norm(XI(:,i)'*Amp);
           
           if xon > 100*sst
               disp('Her')
               disp(['Warning: Norm || (A-p*I)*xo ||_2 is larger than 100*sst: ',...
                       num2str(xon(1,1)), '.']);
           end
           
           if xin > 100*sst
               disp(['Warning: Norm || xi^H*(A-p*I) ||_2 is larger than 10*sst: ',...
                       num2str(xin(1,1)), '.']);
           end

%
%   Kjetil Havre 16/12 - 1995.
%
	 
           yph = C*XO(:,i); nyph = norm(yph);
           if nyph > sst
               O(:,i) = yph/norm(yph);
               XO(:,i) = XO(:,i)/norm(yph);	 
               if isreal(O(:,i)) == 0
                   inz = find(abs(O(:,i))>sst); 
                   vh = angle(O(inz(1),i));
                   O(:,i) = O(:,i)*exp(-sqrt(-1)*vh); 
                   XO(:,i)=XO(:,i)*exp(-sqrt(-1)*vh);
               end
           else
               O(:,i) = yph;  
           end
           
           uph=B'*XI(:,i); nuph = norm(uph);
%
%        Comment: Almost always is the B matrix real, however when
%                 facorizing first on RHP pole then the B matrix 
%                 becomes complex if the pole is complex.
%        in order that up^H*G^{-1}(p) = 0 we need to define
%        u_p = B^H x_p
%        This make on a difference when B is complex.
%        Kjetil Havre 15/2 - 1996
%
	       if nuph > sst
               I(:,i) = uph/norm(uph); 
               XI(:,i)=XI(:,i)/norm(uph);
               if isreal(I(:,i)) == 0
                   inz=find(abs(I(:,i))>sst); 
                   vh = angle(I(inz(1),i));
                   I(:,i)=I(:,i)*exp(-sqrt(-1)*vh); 
                   XI(:,i)=XI(:,i)*exp(-sqrt(-1)*vh);
               end
           else
               I(:,i) = uph;
           end
           
           wp = P(i)/sqrt(-1)+max(10*eps,0.01*sst);
           Gp = frsp(G, wp);
           ypn = vnorm( mmult(vpinv(Gp), O(:,i)) );
           upn = vnorm( mmult(I(:,i)',vpinv(Gp)) );
           
           if ypn(1,1) > 100*sst
               disp(['Warning: Norm || G^{-1}(p)*yp ||_2 is larger than 100*sst: ',...
                       num2str(ypn(1,1)), '.']);
           end
           
           if upn(1,1) > 100*sst
               isrB = isreal(B);
               if isrB
                   %	          disp('B is real')
               else 
                   disp('B is not real')
                   B = B
               end
               
               disp(['Warning: Norm || up^H*G^{-1}(p) ||_2 is larger than 100*sst: ',...
                       num2str(upn(1,1)), '.']);
           end
           
       else
           disp(['Error:   Matrix (I*p-A) is not singular for i=',...
                   int2str(i),', P(i)=', num2str(P(i)), '.'] );
           disp(['         Smallest singular value is: ', num2str(S(nx,nx)),'. Storing zero vectors.'] );
           O  = [O zeros(ny,1)]; 
           I  = [I zeros(nu,i)];
           XO(:,i) = zeros(nx,1);  
           XI(:,i) = zeros(nx,1);   
           err = 1;
       end
   end
end

   
