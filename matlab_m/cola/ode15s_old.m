function [tout,yout,stats] = ode15s(ydot,tspan,y0,options)
%ODE15S	Solve stiff differential equations, variable order method.
%	[T,Y] = ODE15S('ydot',TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates
%	the system of first order differential equations y' = ydot(t,y) from
%	time T0 to TFINAL with initial conditions Y0.  Function ydot(t,y)
%	must return a column vector.  Each row in solution matrix Y
%	corresponds to a time returned in column vector T.  To obtain
%	solutions at the specific times T0, T1, ..., TFINAL (all increasing
%	or all decreasing), use TSPAN = [T0 T1 ... TFINAL].
%	
%	[T,Y] = ODE15S('ydot',TSPAN,Y0,OPTIONS) solves as above with default
%	integration parameters replaced by values in OPTIONS, an argument
%	created with the ODESET function.  See ODESET for details.  Commonly
%	used options are scalar relative error tolerance 'rtol' (1e-3 by
%	default) and vector of absolute error tolerances 'atol' (all
%	components 1e-6 by default).
%	
%	It is possible to specify tspan, y0 and options in ydot.  If TSPAN
%	or YO is empty, or if ODE15S is invoked as ODE15S('ydot'), ODE15S
%	calls [tspan,y0,options] = ydot([],[]) to obtain any values not
%	supplied at the command line.  TYPE CHM6EX to see how this is coded.
%	
%	The Jacobian matrix J(t,y) is critical to the reliability and
%	efficiency of the integration.  If J(t,y) is constant and/or sparse,
%	use the 'constantJ' and/or 'sparseJ' options (see B5EX, BRUSSEX).
%	If ydot is coded so that ydot([t1 t2 ...],[y1 y2 ...]) returns
%	[ydot(t1,y1) ydot(t2,y2) ...], setting 'vectorized' true may speed
%	up the computation of J (see VDPEX).  If an M-file function that
%	evaluates analytically J(t,y) is available, use the 'analyticJ'
%	option (see VDPJAC, BRUSSJAC).
%	
%	As an example, the command
%	
%	    ode15s('vdpex',[0 3000],[2 0]);
%	
%	solves the system y' = vdpex(t,y) with the default relative error
%	tolerance 1e-3 and the default absolute tolerance of 1e-6 for each
%	component.  When called with no output arguments, as in this
%	example, ODE15S calls the default output function ODEPLOT to plot
%	the solution as it is computed.
%	
%	ODE15S also solves problems M(t)*y' = ydot(t,y) with a mass matrix
%	M(t) that is nonsingular and (usually) sparse.  Use the 'mass'
%	option to supply the name of a function of t that returns M(t) (see
%	FEM1EX).  If M(t) is a constant matrix M, use the 'constantM'
%	option, or supply M as the value of 'mass' (see FEM2EX).
%	
%	See also
%	    other ODE solvers:   ODE23S, ODE45, ODE23, ODE113
%	    options handling:    ODESET, ODEGET
%	    output functions:    ODEPLOT, ODEPHAS2, ODEPHAS3
%	    ydot examples:       VDPEX, BRUSSEX, B5EX, CHM6EX
%	    ydot Jacobians:      VDPJAC, BRUSSJAC
%	    Jacobian functions:  NUMJAC, COLGROUP
%	    mass matrices:	 FEM1EX, FEM1MASS, FEM2EX, FEM2MASS

%	ODE15S is a quasi-constant step size implementation in terms of
%	backward differences of the Klopfenstein-Shampine family of
%	Numerical Differentiation Formulas of orders 1-5.  The natural
%	"free" interpolants are used.  Local extrapolation is not done.  By
%	default, Jacobians are generated numerically.  Details are to be
%	found in The MATLAB ODE Suite, L. F. Shampine and M. W. Reichelt,
%	Rept. 94-6, Math. Dept., SMU, Dallas, TX, 1994.

%	Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%	Copyright (c) 1984-95 by The MathWorks, Inc.

false = 0;
true = ~false;

nsteps = 0; 				% stats
nfailed = 0; 				% stats
nfevals = 0; 				% stats
npds = 0; 				% stats
ndecomps = 0; 				% stats
nsolves = 0; 				% stats

if nargin == 1
  tspan = []; y0 = []; options = [];
elseif nargin == 2
  y0 = []; options = [];
elseif nargin == 3
  options = [];
end

% Get default tspan and y0 from ydot if none are specified.
if (length(tspan) == 0) | (length(y0) == 0)
  [ydot_tspan,ydot_y0,ydot_options] = feval(ydot,[],[]);
  if length(tspan) == 0
    tspan = ydot_tspan;
  end
  if length(y0) == 0
    y0 = ydot_y0;
  end
  if length(options) == 0
    options = ydot_options;
  else
    options = odeset(ydot_options,options);
  end
end

% Test that tspan is internally consistent.
tspan = tspan(:);
ntspan = length(tspan);
if ntspan == 1
  t = 0;
  next = 1;
else
  t = tspan(1);
  next = 2;
end
tfinal = tspan(ntspan);
if t == tfinal
  error('The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t);
if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan-1)) <= 0)
  error('The entries in tspan must strictly increase or decrease.');
end

y = y0(:);
neq = length(y);
zerosneq = zeros(neq,1);

% Get options, and set defaults.
rtol = odeget(options,'rtol',1e-3);
atol = odeget(options,'atol',1e-6);
atol = atol(:);
if (rtol <= 0) | any(atol <= 0)  
  error('Tolerance rtol and all entries of atol must be positive.');
end
if rtol < 100 * eps 
  rtol = 100 * eps;
  fprintf('Warning: rtol has been increased to %e\n',rtol);
end
if length(atol) == 1
  atol = atol + zerosneq;
elseif length(atol) ~= neq
  error(['Vector atol must be same length as y0 vector (' num2str(neq) ').']);
end
threshold = atol ./ rtol;

userhmax = abs(odeget(options,'hmax',0.1*abs(tfinal-t)));
hmax = min(userhmax, abs(tfinal-t));

hmin = 4 * eps * abs(t);
userhmin = abs(odeget(options,'hmin',hmin));
hmin = min(max(hmin, userhmin), hmax);

stopfun = odeget(options,'stopfun');
if nargout ~= 0
  outfun = odeget(options,'outfun');
else
  outfun = odeget(options,'outfun','odeplot');
end

refine = odeget(options,'refine',1);

printstats = odeget(options,'printstats',false);

vectorized = odeget(options,'vectorized',false);
constantJ = odeget(options,'constantJ',false);
Js = odeget(options,'sparseJ');
analyticJ = odeget(options,'analyticJ');
if length(analyticJ) == 0
  notanalyticJ = true;
else
  notanalyticJ = false;
end

mass = odeget(options,'mass');
if length(mass) ~= 0
  if isstr(mass)
    Mt = feval(mass,t);
    constantM = odeget(options,'constantM',false);
  else
    Mt = mass;
    constantM = odeget(options,'constantM',true);
    if ~constantM
      error('Mass matrix has constant value, and yet constantM is set false.');
    end
  end
  [L,U] = lu(Mt);
else
  Mt = sparse((1:neq)',(1:neq)',1,neq,neq);
  L = Mt;
  U = Mt;
  constantM = true;
end
if constantM
  Mtnew = Mt;
end

maxk = odeget(options,'maxorder',5);
bdf = odeget(options,'bdf',false);

% Initialize the output function.
if length(outfun) == 0
  haveoutfun = false;
else
  haveoutfun = true;
  feval(outfun,[t tfinal],y,'init');
end

% Allocate memory if we're generating output.
if nargout ~= 0
  if (ntspan == 1) & (refine == 0) 	% only 1 output at tfinal
    nout = 0;
  else
    if ntspan > 2 			% output only at tspan points
      tout = zeros(ntspan,1);
      yout = zeros(ntspan,neq);
    else 				% alloc in chunks
      chunk = max(ceil(128 / neq),refine);
      tout = zeros(chunk,1);
      yout = zeros(chunk,neq);
    end
    nout = 1;
    tout(nout) = t;
    yout(nout,:) = y.';
  end
end

% Initialize method parameters.
G = [1; 3/2; 11/6; 25/12; 137/60];
if bdf
  alpha = [0; 0; 0; 0; 0];
else
  alpha = [-37/200; -1/9; -0.0823; -0.0415; 0];
end
invGa = 1 ./ (G .* (1 - alpha));
erconst = alpha .* G + (1 ./ (2:6)');
difU = [ -1, -2, -3, -4,  -5; 		% difU is its own inverse!
          0,  1,  3,  6,  10;
          0,  0, -1, -4, -10;
          0,  0,  0,  1,   5;
          0,  0,  0,  0,  -1 ];
maxK = 1:maxk;
[kJ,kI] = meshgrid(maxK,maxK);
difU = difU(maxK,maxK);
maxit = 4;

% Set the output flag.
if (ntspan > 2) | (refine == 0)
  outflag = 1; 				% output only at tspan points
elseif refine == 1
  outflag = 2; 				% computed points, no refinement
else
  outflag = 3; 				% computed points, with refinement
  ones1r = ones(1,refine-1);
  S = 1 / refine;
  S = cumsum(S(ones1r)) - 1;
  p = cumprod((S(ones(maxk,1),:)+kI(:,1:refine-1)-1) ./ kI(:,1:refine-1));
end

% Compute the partial derivatives dfdt and dfdy.  The error of
% BDF1 is 0.5*h^2*y''(t).  Use this to determine the optimal h.
f0 = feval(ydot,t,y);
[m,n] = size(f0);
if n > 1
  error(['Function ' ydot '(t,y) must return a column vector.'])
elseif m ~= neq
  error(['Vector ' ydot '(t,y) must be same length as initial conditions.']);
end
F0 = U \ (L \ f0);
wt = abs(y) + threshold;
absh = min(hmax, abs(tspan(next) - t));
rh = norm(F0 ./ wt,inf) / (0.8 * sqrt(rtol));
if absh * rh > 1
  absh = 1 / rh;
end
absh = max(absh, hmin);
h = tdir * absh;

tdel = (t + tdir*min(sqrt(eps)*max(abs(t),abs(t+h)),absh)) - t;
f1 = feval(ydot,t+tdel,y);
nfevals = nfevals + 2; 			% stats
dfdt = (f1 - f0) ./ tdel;

if notanalyticJ
  [dfdy,fac,g,nF] = numjac(ydot,t,y,f0,atol,[],vectorized,Js,[]);
  nfevals = nfevals + nF; 		% stats
else
  dfdy = feval(analyticJ,t,y);
end
npds = npds + 1; 			% stats
JMcurrent = true;

absh = min(hmax, abs(tspan(next) - t));
rh = 1.25 * sqrt(0.5*norm((U \ (L \ (dfdt + dfdy*F0))) ./ wt,inf) / rtol);
if absh * rh > 1
  absh = 1 / rh;
end
absh = max(absh, hmin);
h = tdir * absh;

% Initialize.
hnew = absh;
k = 1; 					% start at order 1 with BDF1
K = 1; 					% K = 1:k
newhk = false;

dif = zeros(neq,maxk+2);
dif(:,1) = h * F0;

hinvGak = h * invGa(k);
nconhk = 0; 				% steps taken with current h and k
[L,U] = lu(Mt - hinvGak * dfdy);
ndecomps = ndecomps + 1; 		% stats
havrate = false;

% THE MAIN LOOP

last = false; 				% is it the last step?
while ~last
  
  hmin = max(userhmin, 4 * eps * abs(t));
  if 1.1*hnew >= abs(tfinal - t)
    hnew = abs(tfinal - t);
    newhk = true;
    last = true;
  end

  if newhk
    newhk = false;
    
    difRU = cumprod((kI - 1 - kJ*(hnew/absh)) ./ kI) * difU;
    dif(:,K) = dif(:,K) * difRU(K,K);
    absh = hnew;
    h = tdir * absh;

    hinvGak = h * invGa(k);
    nconhk = 0;
    [L,U] = lu(Mt - hinvGak * dfdy);
    ndecomps = ndecomps + 1; 		% stats
    havrate = false;
  end
  
  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true; 			% no failed attempts
  while true
    
    gotynew = false; 			% is ynew evaluated yet?
    while ~gotynew

      % Compute the constant terms in the equation for ynew.
      psi = (dif(:,K) * G(K)) * invGa(k);

      % Predict a solution at t+h.
      tnew = t + h;
      if k == 1
	pred = y + dif(:,1);
      else
	pred = y + sum(dif(:,K).').';
      end
      ynew = pred;
      
      % The difference, difkp1, between pred and the final accepted 
      % ynew is equal to the backward difference of ynew of order
      % k+1. Initialize to zero for the iteration to compute ynew.
      difkp1 = zerosneq; 
      invwt = 1 ./ (max(abs(y), abs(ynew)) + threshold);
      minnrm = 100*eps*norm(ynew .* invwt,inf);

      % Mtnew is required in the RHS function evaluation.
      if ~constantM
	Mtnew = feval(mass,tnew);
      end

      % Iterate with simplified Newton method
      tooslow = false;
      for iter = 1:maxit
        del = U \ (L \ (hinvGak*feval(ydot,tnew,ynew) -  Mtnew*(psi+difkp1)));
        newnrm = norm(del .* invwt,inf);
        difkp1 = difkp1 + del;
        ynew = pred + difkp1;
 
        if newnrm <= minnrm
          gotynew = true;
          break;
	elseif iter == 1
          if havrate
            errit = newnrm * rate / (1 - rate) ;
	    if errit <= 0.05*rtol 	% More stringent when using old rate
              gotynew = true;
              break;
            end
	  else
	    rate = 0;
	  end
	elseif newnrm > 0.9*oldnrm
	  tooslow = true;
	  break;
	else
	  rate = max(0.9*rate, newnrm / oldnrm);
	  havrate = true;                 
	  errit = newnrm * rate / (1 - rate);
	  if errit <= 0.5*rtol             
	    gotynew = true;
	    break;
	  elseif iter == maxit            
	    tooslow = true;
	    break;
	  elseif 0.5*rtol < errit*rate^(maxit-iter)
	    tooslow = true;
	    break;
	  end
	end
	
	oldnrm = newnrm;
      end 				% end of Newton loop
      nfevals = nfevals + iter; 	% stats
      nsolves = nsolves + iter; 	% stats
      
      if tooslow
	nfailed = nfailed + 1; 		% stats
	if absh <= hmin
	  fprintf('Convergence failure at %e with step size of %e\n',t,h);
	  if nargout ~= 0
	    tout = tout(1:nout);
	    yout = yout(1:nout,:);
	    stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
	  end
	  return;
	end
	% Speed up the iteration by forming J(t,y) or reducing h.
	if ~JMcurrent
	  if ~constantJ
	    if notanalyticJ
	      f0 = feval(ydot,t,y);
	      [dfdy,fac,g,nF] = numjac(ydot,t,y,f0,atol,fac,vectorized,Js,g);
	      nfevals = nfevals + nF + 1; % stats
	    else
	      dfdy = feval(analyticJ,t,y);
	    end
	    npds = npds + 1; 		% stats
	  end
	  JMcurrent = true;
	else
	  hnew = max(0.3 * absh, hmin);
	  last = false;

	  difRU = cumprod((kI - 1 - kJ*(hnew/absh)) ./ kI) * difU;
	  dif(:,K) = dif(:,K) * difRU(K,K);
	  absh = hnew;
	  h = tdir * absh;
	  
	  hinvGak = h * invGa(k);
	  nconhk = 0;
	end
	[L,U] = lu(Mt - hinvGak * dfdy);
	ndecomps = ndecomps + 1; 	% stats
	havrate = false;
      end   
    end     % end of while loop for getting ynew
    
    % difkp1 is now the backward difference of ynew of order k+1.
    err = norm(difkp1 .* invwt,inf) * erconst(k);
    
    if err > rtol 			% Failed step
      nfailed = nfailed + 1; 		% stats
      if absh <= hmin
	fprintf('Step failure at %e with a minimum step size of %e\n', t, h);
	if nargout ~= 0
	  tout = tout(1:nout);
	  yout = yout(1:nout,:);
	  stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
	end
	return;
      end
      
      if nofailed
	nofailed = false;
	hnew = absh * max(0.1, 0.833*(rtol/err)^(1/(k+1))); % 1/1.2
	if k > 1
	  errkm1 = norm((dif(:,k) + difkp1) .* invwt,inf) * erconst(k-1);
	  hkm1 = absh * max(0.1, 0.769*(rtol/errkm1)^(1/k)); % 1/1.3
	  if hkm1 > hnew
	    k = k - 1;
	    K = 1:k;
	    hnew = hkm1;
	  end
	end
      else
	hnew = 0.5 * absh;
      end
      hnew = max(hnew, hmin);
      last = false;

      difRU = cumprod((kI - 1 - kJ*(hnew/absh)) ./ kI) * difU;
      dif(:,K) = dif(:,K) * difRU(K,K);
      absh = hnew;
      h = tdir * absh;

      hinvGak = h * invGa(k);
      nconhk = 0;
      [L,U] = lu(Mt - hinvGak * dfdy);
      ndecomps = ndecomps + 1; 		% stats
      havrate = false;
      
      if haveoutfun
	if feval(outfun,tnew,ynew,'failed') == 1
	  if nargout ~= 0
	    tout = tout(1:nout);
	    yout = yout(1:nout,:);
	    stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
	  end
	  return;
	end
      end
      
    else 				% Successful step
      break;
      
    end
  end

  dif(:,k+2) = difkp1 - dif(:,k+1);
  dif(:,k+1) = difkp1;
  for j = k:-1:1
    dif(:,j) = dif(:,j) + dif(:,j+1);
  end
  
  if nargout ~= 0
    if outflag == 2 			% computed points, no refinement
      nout = nout + 1;
      if nout > length(tout)
	tout = [tout; zeros(chunk,1)];
	yout = [yout; zeros(chunk,neq)];
      end
      tout(nout) = tnew;
      yout(nout,:) = ynew.';
    elseif outflag == 3 		% computed points, with refinement
      i = nout+1:nout+refine-1;
      nout = nout + refine;
      if nout > length(tout)
	tout = [tout; zeros(chunk,1)]; 	% requires chunk >= refine
	yout = [yout; zeros(chunk,neq)];
      end
      tout(i) = tnew + h*S';
      yout(i,:) = (ynew(:,ones1r) + dif(:,K) * p(K,:)).';
      tout(nout) = tnew;
      yout(nout,:) = ynew.';
    elseif outflag == 1 		% output only at tspan points
      oldnout = nout;
      while true
	if next > ntspan 
	  break;
	elseif tdir * (tnew - tspan(next)) < 0
	  break;
	end
	nout = nout + 1; 		% tout and yout are already allocated
	tout(nout) = tspan(next);
	s = (tspan(next) - tnew) / h;
	yout(nout,:) = (ynew + dif(:,K) * cumprod((s+K-1) ./ K).').';
	next = next + 1;
      end
    end
    
    if haveoutfun
      if outflag == 3
	for i = nout-refine+1:nout-1
	  feval(outfun,tout(i),yout(i,:).','interp');
	end
      elseif outflag == 1
	for i = oldnout+1:nout
	  feval(outfun,tout(i),yout(i,:).','interp');
	end
      end	
      if feval(outfun,tnew,ynew) == 1
	tout = tout(1:nout);
	yout = yout(1:nout,:);
	stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
	return;
      end
    end
    
  elseif haveoutfun
    if outflag == 3 			% computed points, with refinement
      for i = 1:refine-1
	yinterp = ynew + dif(:,K) * cumprod((S(i)+K-1) ./ K).';
	feval(outfun,tnew + h*S(i),yinterp,'interp');
      end
    elseif outflag == 1 		% output only at tspan points
      while true
	if next > ntspan 
	  break;
	elseif tdir * (tnew - tspan(next)) < 0
	  break;
	end
	s = (tspan(next) - tnew) / h;
	yinterp = ynew + dif(:,K) * cumprod((s+K-1) ./ K).';
	feval(outfun,tspan(next),yinterp,'interp');
	next = next + 1;
      end
    end
    if feval(outfun,tnew,ynew) == 1
      return;
    end
  end
  
  nconhk = min(nconhk+1,maxk+2);
  if nconhk >= k + 2
    hopt = min(hmax,absh/max(0.1,1.2*(err/rtol)^(1/(k+1))));
    kopt = k;
    if k > 1
      errkm1 = norm(dif(:,k) .* invwt,inf) * erconst(k-1);
      hkm1 = min(hmax,absh/max(0.1,1.3*(errkm1/rtol)^(1/k)));
      if hkm1 > hopt 
	hopt = hkm1;
	kopt = k - 1;
      end
    end
    if k < maxk
      errkp1 = norm(dif(:,k+2) .* invwt,inf) * erconst(k+1);
      hkp1 = min(hmax,absh/max(0.1,1.4*(errkp1/rtol)^(1/(k+2))));
      if hkp1 > hopt 
	hopt = hkp1;
	kopt = k + 1;
      end
    end
    if hopt > absh
      hnew = hopt;
      if k ~= kopt
	k = kopt;
	K = 1:k;
      end
      newhk = true;
    end
  end
  
  % Advance the integration one step.
  t = tnew;
  y = ynew;
  nsteps = nsteps + 1; 			% stats
  JMcurrent = constantJ & constantM;
  if ~constantM
    Mt = Mtnew;
  end
  
end

if tdir * (tfinal - tnew) > 0
   fprintf('Singularity likely near %g.\n',tnew)
end

if printstats 				% print cost statistics
  fprintf('%g successful steps\n', nsteps);
  fprintf('%g failed attempts\n', nfailed);
  fprintf('%g calls to ydot\n', nfevals);
  fprintf('%g partial derivatives\n', npds);
  fprintf('%g LU decompositions\n', ndecomps);
  fprintf('%g solutions of linear systems\n', nsolves);
end

if nargout ~= 0
  tout = tout(1:nout);
  yout = yout(1:nout,:);
  stats = [nsteps; nfailed; nfevals; npds; ndecomps; nsolves];
end
