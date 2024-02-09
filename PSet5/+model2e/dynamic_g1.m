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
    T = model2e.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(4, 10);
g1(1,1)=(-(1-params(4)+T(3)*y(6)*getPowerDeriv(y(1),params(3),1)));
g1(1,3)=1;
g1(1,4)=1;
g1(1,5)=(-(T(2)*getPowerDeriv(y(5),1-params(3),1)));
g1(1,6)=(-(T(1)*T(3)));
g1(2,3)=(-(T(16)*T(19)*params(3)*y(9)*getPowerDeriv(y(3),params(3)-1,1)));
g1(2,4)=T(27);
g1(2,7)=(-(T(20)*params(2)*(T(15)*T(11)*params(6)*getPowerDeriv(y(7),params(6)-1,1)-T(12)*T(11)*getPowerDeriv(y(7),params(6),1)*T(28))/(T(15)*T(15))));
g1(2,5)=(T(8)*T(4)*T(29)-T(6)*T(26)*T(7)*T(29))/(T(8)*T(8));
g1(2,8)=(-(T(20)*params(2)*(T(15)*T(10)*T(30)-T(12)*T(28)*T(13)*T(30))/(T(15)*T(15))+T(16)*T(18)*getPowerDeriv(y(8),1-params(3),1)));
g1(2,9)=(-(T(16)*T(19)*params(3)*T(17)));
g1(3,3)=(-(T(24)*(1-params(3))*y(6)*T(9)*getPowerDeriv(y(3),params(3),1)));
g1(3,4)=(-((T(22)*(params(6)-1)*T(25)-(params(6)-1)*T(7)*T(21)*T(5)*T(25)*T(26))/(T(22)*T(22))+T(24)*T(23)*(1-params(3))*y(6)*T(27)));
g1(3,5)=(-((-((params(6)-1)*T(7)*(T(21)*T(26)*T(7)*T(29)+T(8)*(-(getPowerDeriv(1-y(5),params(6),1))))))/(T(22)*T(22))+T(24)*T(23)*(1-params(3))*y(6)*(T(8)*T(4)*T(29)-T(6)*T(26)*T(7)*T(29))/(T(8)*T(8))+(1-params(3))*y(6)*T(9)*T(23)*getPowerDeriv(y(5),(-params(3)),1)));
g1(3,6)=(-(T(24)*T(23)*(1-params(3))*T(9)));
g1(4,2)=(-(params(5)*1/y(2)));
g1(4,6)=1/y(6);
g1(4,10)=(-1);

end
