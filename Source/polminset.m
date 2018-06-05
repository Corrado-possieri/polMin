function options = polminset(varargin)
%POLMINSET Create/alter ODE OPTIONS structure.
%   OPTIONS = POLMINSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. 
%
%POLMINSET PROPERTIES
%   
%maxIter - maximum number of iterations  [ positive integer {1e3} ]
%   
%maxRep - Maximum number of zero incremental operatitons  [ positive integer {1e1}]
%
%AbsTol - Absolute error tolerance  [ positive scalar {1e-3} ]
%   
%RelTol -  Relative error tolerance  [ positive scalar {1e-3} ]
%   
%verbose - verbosity of the method  [ {0} | 1 ]
%   
%p - probability of transverse directions  [ scalar in [0,1] {0.5} ]
%
%dist - probability mass function for the coordinate directions [ vector {'uniform'} ]

if (nargin == 0) && (nargout == 0)
  fprintf('         maxIter: [ positive integer {1e3} ]\n');
  fprintf('          maxRep: [ positive integer {1e1} ]\n');
  fprintf('          AbsTol: [ positive scalar {1e-3} ]\n');
  fprintf('          RelTol: [ positive scalar {1e-3} ]\n'); 
  fprintf('         verbose: [ 1 | {0} ]\n');
  fprintf('               p: [ scalar in [0,1] {0.5} ]\n');
  fprintf('            dist: [ vector  {uniform} ]\n');
  fprintf('\n');
  return;
end

Names = [
    'maxIter         '
    'maxRep          '
    'AbsTol          '
    'RelTol          '
    'verbose         '
    'p               '            
    'dist            '   
    ];
m = size(Names,1);
names = lower(Names);

options.maxIter = 1e3;
options.maxRep = 1e1;
options.AbsTol = 1e-3;
options.RelTol = 1e-3;
options.verbose = 0;
options.p = 0.5;
options.dist = [];

i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg) || (isstring(arg) && isscalar(arg)) % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(message('MATLAB:polminset:NoPropNameOrStruct', i));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end
% Convert string arguments and options.
for ii = 1:nargin
    if isstring(varargin{ii}) && isscalar(varargin{ii})
        varargin{ii} = char(varargin{ii});
    end
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error(message('MATLAB:polminset:ArgNameValueMismatch'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error(message('MATLAB:polminset:NoPropName', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(message('MATLAB:polminset:InvalidPropName', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
            matches = deblank(Names(j(1),:));
        for k = j(2:length(j))'
                matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
        end
            error(message('MATLAB:polminset:AmbiguousPropName',arg,matches));
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(message('MATLAB:polminset:NoValueForProp', arg));
end


end