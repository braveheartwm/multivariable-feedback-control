function o = odeget(options,name,default)
%ODEGET	Extract options from an options argument created with ODESET.
%	O = ODEGET(OPTIONS,NAME) extracts the value of the option specified
%	by string NAME from integrator options argument OPTIONS, returning
%	an empty matrix if the option value is not specified in OPTIONS.
%	Only as many leading characters of an option name as necessary to
%	uniquely identify it need be typed, and case is ignored.
%	
%	O = ODEGET(OPTIONS) displays all option names and their current
%	values for argument OPTIONS.  [] is a valid OPTIONS argument.
%	
%	O = ODEGET(OPTIONS,NAME,DEFAULT) returns O = DEFAULT if the named
%	option is not specified in OPTIONS.
%	
%	See also ODESET, ODE45, ODE23, ODE113, ODE15S, ODE23S.

%	Mark W. Reichelt and Lawrence F. Shampine, 3/1/94
%	Copyright (c) 1984-95 by The MathWorks, Inc.

if nargin == 0
  error('Not enough input arguments');
end

% Support some old SIMULINK style syntax, OPTIONS = [RTOL,HMIN,HMAX].
len = length(options);
if (0 < len) & (len <= 3)
  arg = options;
  options = odeset([],'rtol',[]);
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
end

% Print out current values of options.
if nargin == 1
  fprintf('\tanalyticJ = %s\n',setstr(odeget(options,'analyticJ')));
  
  v = odeget(options,'atol');
  if length(v) > 1
    fprintf('\tatol = [%g',v(1));
    fprintf(' %g',v(2:length(v)));
    fprintf(']\n');
  else
    fprintf('\tatol = %s\n',sprintf('%.1g',v));
  end
  
  fprintf('\tbdf = %s\n',sprintf('%d',odeget(options,'bdf')));
  fprintf('\tconstantJ = %s\n',sprintf('%d',odeget(options,'constantJ')));
  fprintf('\tconstantM = %s\n',sprintf('%d',odeget(options,'constantM')));
  fprintf('\thmax = %s\n',sprintf('%g',odeget(options,'hmax')));
  fprintf('\thmin = %s\n',sprintf('%g',odeget(options,'hmin')));
  
  v = odeget(options,'mass');
  if length(v) == 0
    fprintf('\tmass =\n');
  elseif isstr(v)
    fprintf('\tmass = %s\n',v);
  else
    [r,c] = size(v);
    if issparse(v)
      fprintf('\tmass = [ (sparse %d by %d) ]\n',r,c);
    else
      fprintf('\tmass = [ (%d by %d) ]\n',r,c);
    end
  end
    
  fprintf('\tmaxorder = %s\n',sprintf('%d',odeget(options,'maxorder')));
  fprintf('\toutfun = %s\n',setstr(odeget(options,'outfun')));
  fprintf('\trefine = %s\n',sprintf('%d',odeget(options,'refine')));
  fprintf('\trtol = %s\n',sprintf('%.1g',odeget(options,'rtol')));

  v = odeget(options,'sparseJ');
  if length(v) == 0
    fprintf('\tsparseJ =\n');
  else
    [r,c] = size(v);
    fprintf('\tsparseJ = [ (sparse %d by %d) ]\n',r,c);
  end
  
  fprintf('\tstats = %s\n',sprintf('%d',odeget(options,'stats')));
  fprintf('\tstopfun = %s\n',setstr(odeget(options,'stopfun')));
  fprintf('\tvectorized = %s\n',sprintf('%d',odeget(options,'vectorized')));
  fprintf('\n');
  return;
end

if length(options) == 0
  if nargin == 3
    o = default;
  else
    o = [];
  end
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
strings = [1 10 15];

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

j = strmatch(lower(name),externals(:,2:ne));
if length(j) == 0
  error(['odeget: unrecognized option name, ' name]);
elseif length(j) > 1
  error(['odeget: ambiguous option name, ' name]);
end

inans = [0 find(isnan(options))];
i = externals(j,1);
o = options(inans(i)+1:inans(i+1)-1);

if (length(o) ~= 0) & any(i == strings)
  o = setstr(o);
elseif i == 2 				% atol
  o = o(:);
elseif i == 8 				% mass
  len = length(o);
  if len > 0
    if o(1) == 0 			% string
      o = setstr(o(2:len));
    elseif o(1) == 1 			% sparse matrix
      if len > 4
	I = o(4+(1:o(4)));
	J = o(4+o(4)+(1:o(4)));
	s = o(4+2*o(4)+(1:o(4)));
      else
	I = [];
	J = [];
	s = [];
      end
      o = sparse(I,J,s,o(2),o(3));
    elseif o(1) == 2 			% matrix
      if len > 3
	o = reshape(o(4:len),o(2),o(3));
      else
	o = [];
      end
    end
  end
elseif i == 13 				% sparseJ
  len = length(o);
  if len > 0
    onew = sparse([],[],[],o(1),o(2),len-2);
    if len > 2
      s = o(3:len);
      onew(s) = ones(size(s));
    end
    o = onew;
  end
end

if (nargin == 3) & (length(o) == 0)
  o = default;
end
