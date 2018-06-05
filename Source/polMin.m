function [minVal, xOpt] = polMin(f,g,x,vars,options)
%POLMIN --- solve a polynomial optimization problem 
%
% POLMIN attempts to solve the polynomial optimization problem 
%   min F(X) subject to G(X) <= 0
% by using the transverse block descend optimization method 
%
% [MINVAL, XOPT] = POLMIN(F) attempts to solve the unconstrained
% polynomial optimization problem 
%   min F(X)
% The function returns the estimated minimum value MINVAL of F and the 
% estimated minimum XOPT. 
%
%
% [MINVAL, XOPT] = POLMIN(F,G) attempts to solve the constrained
% polynomial optimization problem 
%   min F(X) subject to G(X) <= 0
% The function returns the estimated minimum value MINVAL of F
% and the estimated minimum XOPT.
%
% [MINVAL, XOPT] = POLMIN(F,G,X0) attempts to solve the constrained
% polynomial optimization problem 
%   min F(X) subject to G(X) <= 0
% The function returns the estimated minimum value MINVAL of F and the 
% estimated minimum XOPT.
%
% [MINVAL, XOPT] = POLMIN(F,G,X0,V) attempts to solve the constrained
% polynomial optimization problem 
%   min F(X) subject to G(X) <= 0
% starting at X0. The vector V specifies the order of the symbolic 
% variables. The function returns the estimated minimum value MINVAL 
% of F and the estimated minimum XOPT.
%
% [MINVAL, XOPT] = POLMIN(F,G,X0,V,OPTIONS) attempts to solve the 
% constrained polynomial optimization problem 
%   min F(X) subject to G(X) <= 0
% starting at X0. The vector V specifies the order of the symbolic 
% variables. The function returns the estimated minimum value MINVAL 
% of F and the estimated minimum XOPT.
% OPTIONS specifies the parameters of the method as
% options.maxIter  (max number of iterations,               defult = 1e3)
% options.maxRep   (max repetions without increment,        default = 1e1)
% options.AbsTol   (absolute tolerance,                     default = 1e-3)
% options.RelTol   (relative tolerance,                     default = 1e-3)
% options.verbose  (verbosity of the method,                defaul = 0)
% options.p        (probability of transverse directions,   default = 0,5)
% options.dist     (distribution of directions,             default = 1/n*ones(1,n))
% The function returns the estimated minimum value MINVAL of F and the 
% estimated minimum XOPT.


switch nargin
    case 1
        p = 0.5;
        g = -1;
        vars = symvar(f);
        x = zeros(length(vars),1);
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 2
        p = 0.5;
        vars = symvar(f);
        x = zeros(length(vars),1);
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 3
        p = 0.5;
        vars = symvar(f);
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 4
        p = 0.5;
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 5
        p = options.p;
        n = length(x);
end

if isempty(options.dist)
    options.dist = 1/n*ones(1,n);
end

if length(x) ~= length(options.dist)
    error 'x and options.dist must have the same dimensions';
end

Pnorm=[0 options.dist]/sum(options.dist);
Pcum=cumsum(Pnorm);

count = 0;
ncount = 0;

valFun = zeros(options.maxIter,1);

subsfun = @(x,s) genScalarProb(f, g, vars, x, s);

for ii = 1:options.maxIter
    r = rand(1);
    if r < p
        R = rand(1);
        [~,ind] = histc(R,Pcum); 
        eee = eye(n);
        s = eee(:,ind);
    else    
        s = randn(n,1);
        s = s./norm(s);
    end
    [ScalPol, ScalConstr] = subsfun(x, s);
    [minu, valFun(ii),mval] = minScalPol(ScalPol, ScalConstr, options);
    x = x + s*minu;
    
    if ii > 1 && valFun(ii) >= valFun(ii-1) - options.AbsTol
        count = count + 1;
    else
        count = 0;
    end
    
    if ii > 1 &&  abs((valFun(ii)-valFun(ii-1))/min(valFun(ii),valFun(ii-1))) < options.RelTol
        ncount = ncount + 1;
    else
        ncount = 0;
    end
    
    if options.verbose == 1
        fprintf('\n|  %3i    |',ii);
        fprintf('|    %+.4e    |', valFun(ii));
        fprintf('|    %+.4e    |', eval(max(subs(ScalConstr,symvar(ScalConstr),minu))));
        fprintf('|    %+.4e    |', eval(max(subs(g,vars,x'))));
        if mval == Inf
            fprintf(' *');
        end
    end
    
    if count == options.maxRep || ncount == options.maxRep
        fprintf('\n\nminimum found\n')
        break;
    end
end

if ii >= options.maxIter
    fprintf('\n\ncomputation interrupted\n')
end

minVal = valFun(ii);
xOpt = x + s*minu;

end