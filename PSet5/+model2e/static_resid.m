function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model2e.static_resid_tt(T, y, x, params);
end
residual = zeros(4, 1);
lhs = y(2)+y(1);
rhs = y(4)*T(1)*T(2)+y(1)*(1-params(4));
residual(1) = lhs - rhs;
lhs = T(8);
rhs = T(8)*params(2)*T(11);
residual(2) = lhs - rhs;
lhs = 0;
rhs = (params(6)-1)*T(6)/T(13)+T(1)*(1-params(3))*y(4)*T(8)*T(14);
residual(3) = lhs - rhs;
lhs = log(y(4));
rhs = log(y(4))*params(5)+x(1);
residual(4) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
