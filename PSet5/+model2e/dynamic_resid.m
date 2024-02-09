function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = model2e.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(4, 1);
lhs = y(4)+y(3);
rhs = T(2)*T(3)+y(1)*(1-params(4));
residual(1) = lhs - rhs;
lhs = T(9);
rhs = T(16)*T(20);
residual(2) = lhs - rhs;
lhs = 0;
rhs = (params(6)-1)*T(7)/T(22)+(1-params(3))*y(6)*T(9)*T(23)*T(24);
residual(3) = lhs - rhs;
lhs = log(y(6));
rhs = params(5)*log(y(2))+x(it_, 1);
residual(4) = lhs - rhs;

end
