function [skm, bnds,rowd,sens,rowp,rowg] = skrpmu(matin,blk,tol,opt,maxit)
%
%  [skm, bnds,rowd,sens,rowp,rowg] = skrpmu(matin,blk,tol,opt)
%
%  Function skrpmu computes skewed mu by scaling the part of the plant
%  corresponding to the last perturbation block so that mu(Km*M) = 1,
%  where Km = diag(I(1), I(2), ... , I(nblk-1), km*I(np) )
%  I(i) denotes the identity of delta block i. 
%  It is assumed that diag(I(1), I(2), ..., I(nblk-1)) is square.
%  The last last delt block can be nonsquare.
%  The skrpmu function iterates on km until th LOWER bound from the 
%  mu(Km*M) calculation is equal to one within the specified tolerance tol.
%  Nessecary requirement is that mu(M11) < 1 for all frequencies with 
%  perturbation equal to the n-1 first blocks.
%
%  Inputs:
%     matin       varying matrix for the frequenct respons.
%     blk         description of perturbation block, same as for mu.
%                 The size of matin must be consistent with the specified
%                 block perturbation blk.
%     tol         tolerance for exiting iteration, default 1e-6.
%     opt         options to be passe on to mu-Tools function mu.
%     maxit       maximum number of iteration, default value 20.
%
%  Outputs:
%     skm = 1/km, skewd mu with last perturbation held constant.
%     bnds,       bnds from the mu calculation, both lower and upper 
%                 should be equal to one.
%     rowd,       rowd from last mu calculation.
%     sens,       sens from last mu calculation.
%     rowp,       rowp from last mu calculation.
%     rowg,       rowg from last mu calculation.
%
%  Note! skewed mu calculation by iterating on mu calculations in 
%  mu tools is very slow. However, it is easy to implement.
%  More efficient algorithms for calculation of skewed mu
%  are likely to be provided by mu tools in next version.
%  
%  Written by Kjetil Havre, e-mail: kjetil@ife.no after proposal  
%  by Sigurd Skogestad, e-mail: skoge@kjemi.unit.no.
%

    if nargin < 5
       maxit = 20;
    end    
    if nargin < 4
       opt = 'lu';
    end
    if nargin < 3
       tol = 1e-6;
    end
    if nargin < 2
       disp( '   Error in skrpmu: To few input arguments.')
       disp( '      [skm, bnds,rowd,sens,rowp,rowg] = skewedRPmu(matin,blk)' )
       return
    end
    

    [nblk,dum] = size(blk);
    blkp = ptrs(blk);
    if dum~=2
        error('Invalid Block structure');
        return
    end

    blkRP = blk(nblk, :);
    nI2 = min(blkRP);
    
    I2 = eye(nI2); 
%     R2=abv(I2,zeros(blkRP(1)-nI2,nI2));
%     R2=sbs(R2, zeros(blkRP(1), blkRP(2)-nI2));
    R2=[I2; zeros(blkRP(1)-nI2,nI2)];
    R2=[R2 zeros(blkRP(1), blkRP(2)-nI2)];
    
    nr1 = sum(blk(1:nblk-1,1));
    R1 = eye(nr1);
    Z1 = zeros(nr1, blkRP(2));
    Z2 = zeros(blkRP(1), nr1);
    
    iv = getiv(matin);
    np = max(size(iv));
    ma  = xtract(matin, iv(1));
    [b, rd, s, rp, rg] = mu(ma,blk, opt); 
    fp = b(1) - 1;xp = 1; 
    x  = 1/b(1);
    
    skm = zeros(np+1, 2);
    [nr, nc] = size( b); bnds=zeros(np+1,nc);
    [nr, nc] = size(rd); rowd=zeros(np+1,nc);
    [nr, nc] = size( s); sens=zeros(np+1,nc);
    [nr, nc] = size(rp); rowp=zeros(np+1,nc);
    if isempty(rg)
       rowg = [];
    else
       [nr, nc] = size(rg); rowg=zeros(np+1,nc);
    end
    
    for i=1:np
       disp(['Frequency point ', int2str(i)] );
       if i > 1
          ma = xtract(matin, iv(i));
          Km = [R1 Z1; Z2 xp*R2];
          [b, rd, s, rp, rg] = mu(mmult(Km,ma),blk,opt); fp = b(1)-1;
       end
       
       Km = [R1 Z1; Z2 x*R2];
       [b, rd, s, rp, rg] = mu(mmult(Km,ma),blk,opt); 
       f = b(1) - 1;
       mr = 1; antit = 1;
       while (abs(f) > tol) & (mr==1)
          xn = x-f*(x-xp)/(f-fp); xn = abs(xn);
	  Km = [R1 Z1; Z2 xn*R2];
	  m = mmult(Km,ma);
          [b, rd, s, rp, rg] = mu(m,blk,opt);
	  fp = f; f = b(1)-1;
	  xp = x; x = xn;
	  if antit >= maxit
	     disp('Warning: Maximum number of iterations, exiting without convergence.')
             mr = 0;
	  end
	  antit = antit + 1;
       end
       disp(['Frequency point ', int2str(i),'. Solution for  km = ', num2str(x), ', residual = ', num2str(f)])
       skm(i,:)   = [1/x iv(i)];
       bnds(i,:) =  b(1,:);
       rowd(i,:) = rd(1,:);
       sens(i,:) =  s(1,:);
       rowp(i,:) = rp(1,:);
       if isempty(rowg) == 0
           rowg(i,:) = rg(1,:);
       end
    end
    skm(np+1,:) =  [np Inf];  
    [nr, nc] = size( b); bnds(np+1,:) = zeros(1, nc); bnds(np+1, nc-1:nc) = [np Inf];
    [nr, nc] = size(rd); rowd(np+1,:) = zeros(1, nc); rowd(np+1, nc-1:nc) = [np Inf];
    [nr, nc] = size( s); sens(np+1,:) = zeros(1, nc); sens(np+1, nc-1:nc) = [np Inf];
    [nr, nc] = size(rp); rowp(np+1,:) = zeros(1, nc); rowp(np+1, nc-1:nc) = [np Inf];
    if isempty(rowg) == 0
       [nr, nc] = size(rg); rowg(np+1,:) = zeros(1, nc); rowg(np+1, nc-1:nc) = [np Inf];
    end
        
