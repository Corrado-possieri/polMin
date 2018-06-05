function [ScalPol,ScalConstr] = genScalarProb(pol,constr,vars,x,dir)
%GENSCALARPROB --- restriction of a scalar polynomial problem to a line
% 
% [SCALPOL,SCALCONSTR] = GENSCALARPROB(POL,CONSTR,VARS,X,DIR) computes
% the restricion of the multivariate optimization problem
%  min POL(X) subject to CONSTR(X) <= 0
% to the line obtained by substituting to the variables in VARS the line
%  X + u*DIR, 
% where u is an auxiliary scalar symbolic variable. 
% The restriction of the problem to the line is given by
%  min SCALPOL(u) subject to SCALCONSTR(u) <= 0

syms('u','real');

newVal = x + dir*u;

ScalPol = subs(pol, vars, newVal');

ScalConstr = subs(constr, vars, newVal');


end