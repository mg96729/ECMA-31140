function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 21);

T = model2e.static_resid_tt(T, y, x, params);

T(15) = getPowerDeriv(y(1),params(3),1);
T(16) = getPowerDeriv(y(2),params(6),1);
T(17) = getPowerDeriv(T(4)*T(6),params(1),1);
T(18) = (T(7)*T(4)*params(6)*getPowerDeriv(y(2),params(6)-1,1)-T(5)*T(4)*T(16)*T(17))/(T(7)*T(7));
T(19) = getPowerDeriv(y(3),1-params(3),1);
T(20) = (-(getPowerDeriv(1-y(3),1-params(6),1)));
T(21) = (T(7)*T(3)*T(20)-T(5)*T(17)*T(6)*T(20))/(T(7)*T(7));

end
