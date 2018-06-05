function [const, xfeas] = findFeas(pg,px,pvars,options)
%FINDFEAS --- find a feasible point  
%
% FINDFEAS attempts to determine a point X such that
%   G(X) <= 0
% by using the transverse block descend optimization method 
%
% [CONST, XFEAS] = FINDFEAS(G) attempts to compute a feasible
% point XFEAS by solving the entry-wise inequality G(X) - CONST <= 0
%
% [CONST, XFEAS] = FINDFEAS(G,X0) aattempts to compute a feasible
% point XFEAS by solving the entry-wise inequality G(X) - CONST <= 0
% by starting the iterations at X0.
%
% [CONST, XFEAS] = FINDFEAS(G,X0,V) aattempts to compute a feasible
% point XFEAS by solving the entry-wise inequality G(X) - CONST <= 0
% by starting the iterations at X0. The vector V specifies the order of 
% the symbolic variables. 
%
% [CONST, XFEAS] = FINDFEAS(G,X0,V,OPTIONS) aattempts to compute a feasible
% point XFEAS by solving the entry-wise inequality G(X) - CONST <= 0
% by starting the iterations at X0. The vector V specifies the order of 
% the symbolic variables. OPTIONS specifies the parameters of the method as
% options.maxIter  (max number of iterations,               defult = 1e3)
% options.maxRep   (max repetions without increment,        default = 1e1)
% options.AbsTol   (absolute tolerance,                     default = 1e-3)
% options.RelTol   (relative tolerance,                     default = 1e-3)
% options.verbose  (verbosity of the method,                defaul = 0)
% options.p        (probability of transverse directions,   default = 0,5)
% options.dist     (distribution of directions,             default = 1/n*ones(1,n))

syms('e','real')

fslac = e;

switch nargin
    case 1
        p = 0.5;
        pvars = symvar(pg);
        px = zeros(length(pvars),1);
        valg0 = max(subs(pg, pvars, px'));
        x = [px; valg0];
        vars = [pvars, e];
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 2
        p = 0.5;
        pvars = symvar(pg);
        valg0 = max(subs(pg, pvars, px'));
        x = [px; valg0];
        vars = [pvars, e];
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 3
        p = 0.5;
        valg0 = max(subs(pg, pvars, px'));
        x = [px; valg0];
        vars = [pvars, e];
        n = length(x);
        options.maxIter = 1e3;
        options.maxRep = 1e1;
        options.AbsTol = 1e-3;
        options.RelTol = 1e-3;
        options.verbose = 0;
        options.dist = 1/n*ones(1,n);
    case 4
        p = options.p;
        valg0 = max(subs(pg, pvars, px'));
        x = [px; valg0];
        vars = [pvars, e];
        n = length(x);
end

if isempty(options.dist)
    options.dist = 1/n*ones(1,n);
end

Pnorm=[0 options.dist]/sum(options.dist);
Pcum=cumsum(Pnorm);

count = 0;
ncount = 0;

m = length(pg);
g = [-2*abs(valg0), pg];
for ii = 1:m+1
    g(ii) = g(ii) - e;
end

valFun = zeros(options.maxIter,1);

subsfun = @(x,s) genScalarProb(fslac, g, vars, x, s);

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
    
    if valFun(ii) <= options.AbsTol
        fprintf('\n\nfeasible point found\n')
        break;
    end
end

if ii >= options.maxIter
    fprintf('\n\ncomputation interrupted\n')
end

const = valFun(ii);
xfeas = eval(x(1:end-1));

end