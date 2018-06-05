function [minu, minval,mval] = minScalPol(pol,constr,options)
%MINSCALPOL --- minimize a scalar polynomial problem
% 
% [MINU, MINVAL] = MINSCALPOL(POL,COSTR,options) computes the minimum value 
% MINVAL and the corresponding point MINU of the univariate polynomial 
% optimization problem
%  min POL(X) subject to CONSTR(X) <= 0

coeffPol = sym2poly(pol);
deg = length(coeffPol)-1;
coeffDPol = (deg:-1:1).*(coeffPol(1:end-1));

rDPol = roots(coeffDPol);
rDPol = rDPol(imag(rDPol)==0);

crictPoints = [0; rDPol];

m = length(constr);
if m ~= 0
    for ii = 1:m
        conII = sym2poly(constr(ii));
        rGII = roots(conII);
        rGII = rGII(imag(rGII)==0);
        crictPoints = [crictPoints; rGII];
    end
end

ell = length(crictPoints);

critVals = zeros(ell,1);
for ii = 1:ell
    critVals(ii) = polyval(coeffPol,crictPoints(ii));
end

[sortedVals, inds] = sort(critVals);
sortedPoints = crictPoints(inds);

for ii = 1:ell
    mval = max(eval(subs(constr,symvar(constr),sortedPoints(ii))));
    if mval <= options.AbsTol
        minval = sortedVals(ii);
        minu = sortedPoints(ii);
        break;
    end
end

if mval > options.AbsTol
    mval = Inf;
    minu = 0;
    minval = eval(subs(pol,symvar(pol),0));
end

% fMinF = coeffPol;
% fMinF(end) = fMinF(end)-sortedVals(1)+2*(1+abs(sortedVals(1)));
% rootFMinF = roots(fMinF);
% rootFMinF = rootFMinF(imag(rootFMinF) == 0);
% if ~isempty(rootFMinF)
%     nRoots = length(rootFMinF);
%     for jj = 1:nRoots
%         flag = 0;
%         for gg = 1:m
%             conII = sym2poly(constr(gg));
%             valg = polyval(conII,rootFMinF(jj));
%             if valg > 0
%                 flag = 1;
%                 break;
%             end
%         end
%         if flag == 0
%             error 'the problem is unbounded'
%         end
%     end
% end
% 
% flag = 0;

% if ii == ell && flag == 1
%     valuesg = zeros(1,m);
%     for gg = 1:m
%         conII = sym2poly(constr(gg));
%         valuesg(gg) = polyval(conII,0);
%     end
%     disp(valuesg)
%     error 'the problem is unfeasible'
% end

end

