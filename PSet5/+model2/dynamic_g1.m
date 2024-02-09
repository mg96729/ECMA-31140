function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = model2.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(4, 10);
g1(1,1)=(-(1-params(3)+T(3)*y(6)*getPowerDeriv(y(1),params(2),1)));
g1(1,3)=1;
g1(1,4)=1;
g1(1,5)=(-(T(2)*getPowerDeriv(y(5),1-params(2),1)));
g1(1,6)=(-(T(1)*T(3)));
g1(2,3)=(-(T(4)*T(7)*params(2)*y(9)*getPowerDeriv(y(3),params(2)-1,1)));
g1(2,4)=(-params(5))/(y(4)*y(4));
g1(2,7)=(-(T(8)*params(1)*(-params(5))/(y(7)*y(7))));
g1(2,8)=(-(T(4)*T(6)*getPowerDeriv(y(8),1-params(2),1)));
g1(2,9)=(-(T(4)*T(7)*params(2)*T(5)));
g1(3,3)=(-(T(10)*(1-params(2))*y(6)*params(5)/y(4)*getPowerDeriv(y(3),params(2),1)));
g1(3,4)=(-(T(10)*T(9)*(1-params(2))*y(6)*(-params(5))/(y(4)*y(4))));
g1(3,5)=(-((-(1-params(5)))/((y(5)-1)*(y(5)-1))+(1-params(2))*y(6)*params(5)/y(4)*T(9)*getPowerDeriv(y(5),(-params(2)),1)));
g1(3,6)=(-(T(10)*T(9)*(1-params(2))*params(5)/y(4)));
g1(4,2)=(-(params(4)*1/y(2)));
g1(4,6)=1/y(6);
g1(4,10)=(-1);

end
