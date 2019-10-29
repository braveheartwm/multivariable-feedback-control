function options = odeset(arg1,arg2,arg3,arg4,arg5,arg6,arg7, ...
    arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18)
%ODESET	Build or change an options argument for an ODE suite integrator.
%	OPTIONS = ODESET('OptionName',OptionValue) creates an integrator
%	options argument OPTIONS in which the named option has the specified
%	value.  OPTIONS = ODESET('Name1',Value1,'Name2',Value2,...) sets
%	multiple option values with a single statement.  For Boolean
%	options, a value need not be supplied, i.e. merely specifying the
%	name of the option will cause it to be set to 1.  Only as many
%	leading characters of an option name as necessary to uniquely
%	identify it need be typed, and case is ignored.
%	
%	OPTIONS = ODESET(OLDOPTS,'Name',Value,...) first sets OPTIONS equal
%	to the integrator options argument OLDOPTS, and then changes the
%	values of any named options in the returned OPTIONS argument.
%	
%	OPTIONS = ODESET(OLDOPTS,OLDOPTS2) combines the options arguments
%	OLDOPTS and OLDOPTS2, with the non-empty option values in OLDOPTS2
%	taking precedence.
%	
%	With no input arguments, ODESET displays all option names and their
%	possible values.
%	
%	The available options depend on the integrator.  All the codes
%	provide for the following.  (Note there are option synonyms.)
%	
%	'rtol' or 'tol' or 'reltol', a positive scalar, is a relative error
%	tolerance that applies to all components of the solution vector.
%	
%	'atol' or 'abstol', a positive vector, are absolute error tolerances
%	that apply to the corresponding components of the solution vector.
%	If a scalar value is specified, it applies to all components.
%	
%	'hmin' or 'minstep', a positive scalar, is the magnitude of the
%	smallest step size that the integrator is permitted to use.
%	
%	'hmax' or 'maxstep', a positive scalar, is the magnitude of the
%	largest step size that the integrator is permitted to use.
%	
%	'outfun', a string, is the name of an installable output function.
%	
%	'refine', an integer, helps determine the quantity of output
%	produced.  If the tspan vector input of the integrator contains more
%	than two elements or if 'refine' has value 0, solutions are produced
%	only at timepoints specified in tspan.  Otherwise, if 'refine' is N
%	(with N > 0), solutions are produced at the initial time t0 and at
%	the end of each time step.  In addition, interpolated solutions are
%	produced at N-1 equally spaced points in the span of each time step.
%	
%	'printstats' or 'stats', a Boolean, specifies whether statistics
%	about the cost of the integration should be displayed.
%	
%	'stopfun' or 'eventfun' or 'gstop', a string, is the name of an
%	event location function.  (NOT YET IMPLEMENTED).
%	
%	The codes for stiff problems benefit from optional information about
%	the Jacobian dF/dy.  The options are not exclusive:
%	
%	'constantJ', a Boolean, states whether the Jacobian dF/dy is
%	constant (see B5EX).
%	
%	'analyticJ', a string, is the name of a function for evaluating
%	analytically the Jacobian dF/dy (see VDPJAC, BRUSSJAC).
%	
%	'sparseJ', a sparse matrix of 0's and 1's, indicates the sparsity
%	pattern of the Jacobian (see BRUSSEX).
%	
%	'vectorized,' a Boolean, should be set true when integrating F(t,y)
%	if F(t,[y1 y2 ...]) returns [F(t,y1) F(t,y2) ...].  This may speed
%	up numerical Jacobian computation.
%	
%	Two options are specific to ODE15S:
%	
%	'maxorder', an integer, is the maximum order formula allowed.
%	
%	'bdf' or 'gear', a Boolean, states whether the Backward Difference
%	Formulas (Gear's methods) are to be used instead of the default
%	Numerical Differentiation Formulas.
%	
%	With the following options, ODE15S and ODE23S have been extended to
%	solve problems of the form M*y' = F(t,y) or M(t)*y' = F(t,y).
%	
%	'mass', a matrix or a string, is either a constant matrix M, or is
%	the name of a function of t that returns the value of mass matrix M.
%	
%	'constantM', a Boolean, specifies whether or not M is a constant
%	matrix.  Only ODE15S can solve problems with a time-varying M. 
%	
%	See also ODEGET, ODE45, ODE23, ODE113, ODE15S, ODE23S.

%	Mark W. Reichelt and Lawrence F. Shampine, 5/6/94
%	Copyright (c) 1984-95 by The MathWorks, Inc.

% Print out possible values of options.
if nargin == 0
  fprintf('\tanalyticJ: [ string ]\n');
  fprintf('\tatol: [ positive scalar or vector ]\n');
  fprintf('\tbdf: [ Boolean (1 | 0) ]\n');
  fprintf('\tconstantJ: [ Boolean (1 | 0) ]\n');
  fprintf('\tconstantM: [ Boolean (1 | 0) ]\n');
  fprintf('\thmax: [ positive scalar ]\n');
  fprintf('\thmin: [ positive scalar ]\n');
  fprintf('\tmass: [ string or matrix ]\n');
  fprintf('\tmaxorder: [ positive integer ]\n');
  fprintf('\toutfun: [ string ]\n');
  fprintf('\trefine: [ non-negative integer ]\n');
  fprintf('\trtol: [ positive scalar ]\n');
  fprintf('\tsparseJ: [ sparse matrix ]\n');
  fprintf('\tstats: [ Boolean (1 | 0) ]\n');
  fprintf('\tstopfun: [ string ]\n')
  fprintf('\tvectorized: [ Boolean (1 | 0) ]\n');
  fprintf('\n');
  return;
end

names = [
    'analyticj '
    'atol      '
    'bdf       '
    'constantj '
    'constantm '
    'hmax      '
    'hmin      '
    'mass      '
    'maxorder  '
    'outfun    '
    'refine    '
    'rtol      '
    'sparsej   '
    'stats     '   
    'stopfun   '
    'vectorized'
    ];
[m,n] = size(names);
booleans = [3 4 5 14 16]; 		% indices of Booleans

% row index into names array, external string
externals = [
    2,  'abstol     '
    1,  'analyticj  '
    2,  'atol       '
    3,  'bdf        '
    8,  'capacitance'
    4,  'constantj  '
    5,  'constantm  '
    15, 'eventfun   '
    3,  'gear       '
    15, 'gstop      '
    6,  'hmax       '
    7,  'hmin       '
    8,  'mass       '
    9,  'maxorder   '
    6,  'maxstep    '
    7,  'minstep    '
    10, 'outfun     '
    14, 'printstats '
    11, 'refine     '
    12, 'reltol     '
    12, 'rtol       '
    13, 'sparsej    '
    14, 'stats      '
    15, 'stopfun    '
    12, 'tol        '
    16, 'vectorized '
    ];
[me,ne] = size(externals);

options = [];
i = 1;
while i <= nargin
  eval(['arg = arg' num2str(i) ';']);

  if isstr(arg)
    break;
  end
  
  len = length(arg);
  if len == 0
    % Do nothing.
    
  elseif len <= 3
    % Support some old SIMULINK style syntax, OPTIONS = [RTOL,HMIN,HMAX].
    if len >= 1
      if arg(1) ~= 0
	options = odeset(options,'rtol',arg(1));
      end
    end
    if len >= 2
      if arg(2) ~= 0
	options = odeset(options,'hmin',arg(2));
      end
    end
    if len >= 3
      if arg(3) ~= 0
	options = odeset(options,'hmax',arg(3));
      end
    end

  elseif length(options) == 0
    options = arg;
    
  else
    for j = 1:m
      jnans = [0 find(isnan(arg))];
      val = arg((jnans(j)+1):(jnans(j+1)-1));
      
      if length(val) ~= 0
	inans = [0 find(isnan(options))];
	before = options(1:inans(j));
	after = options(inans(j+1):length(options));
	options = [before val after];
      end
    end
  end

  i = i + 1;
end

if length(options) == 0
  options = NaN + zeros(1,m);
end

% A finite state machine to parse name-value pairs or Booleans.
expectval = 0; 			% start out expecting a name
while i <= nargin
  eval(['arg = arg' num2str(i) ';']);
    
  if ~expectval
    if ~isstr(arg)
      error(['expected option name, got ' arg(1) '...']);
    else
      arg = lower(arg);
    end
      
    j = strmatch(arg,externals(:,2:ne));
    if length(j) == 0
      error(['unrecognized option name, ' arg]);
    elseif length(j) > 1
      msg = sprintf('ambiguous option name %s\n',arg);
      s = externals(j(1),2:ne);
      s(s == ' ') = [];
      msg = [msg '[' s];
	for k = j(2:length(j))'
	  s = externals(k,2:ne);
	  s(s == ' ') = [];
	  msg = [msg ',' s];
	end
	msg = [msg ']'];
      error(msg);
    end
    inans = [0 find(isnan(options))];
    before = options(1:inans(externals(j,1)));
    after = options(inans(externals(j,1)+1):length(options));
      
    expectval = 1; 			% we usually expect a value next
    if any(externals(j,1) == booleans) 	% but if it's a Boolean flag
      if i < nargin
	eval(['arg = arg' num2str(i+1) ';']);
	if isstr(arg)
	  expectval = 0;
	end
      else
	expectval = 0;
      end
    end
    if ~expectval
      options = [before 1 after];
    end
    i = i + 1;
    
  else
    if length(arg) ~= 0
      if externals(j,1) == 8 		% mass
	% Tag mass, and pack it if it's a matrix.
	if isstr(arg)
	  arg = [0, real(arg)];
	else
	  [mm,nm] = size(arg);
	  if issparse(arg)
	    [is,js,s] = find(arg);
	    arg = [1; mm; nm; length(is); is; js; s].';
	  else
	    arg = [2; mm; nm; arg(:)].';
	  end
	end
	
      elseif externals(j,1) == 13 	% sparsej
	[ms,ns] = size(arg);
	s = find(arg);
	arg = [ms; ns; s].';
	
      elseif isstr(arg)
	arg = real(arg);

      end
    end
    
    options = [before arg after];
    expectval = 0;
    i = i + 1;
      
  end
end

if expectval
  error(sprintf('expected value for option ''%s''',arg));
end
