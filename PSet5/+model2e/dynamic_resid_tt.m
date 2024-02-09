function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 24);

T(1) = y(1)^params(3);
T(2) = y(6)*T(1);
T(3) = y(5)^(1-params(3));
T(4) = params(6)*y(4)^(params(6)-1);
T(5) = (1-y(5))^(1-params(6));
T(6) = T(4)*T(5);
T(7) = y(4)^params(6);
T(8) = (T(5)*T(7))^params(1);
T(9) = T(6)/T(8);
T(10) = params(6)*y(7)^(params(6)-1);
T(11) = (1-y(8))^(1-params(6));
T(12) = T(10)*T(11);
T(13) = y(7)^params(6);
T(14) = T(11)*T(13);
T(15) = T(14)^params(1);
T(16) = params(2)*T(12)/T(15);
T(17) = y(3)^(params(3)-1);
T(18) = params(3)*y(9)*T(17);
T(19) = y(8)^(1-params(3));
T(20) = 1+T(18)*T(19)-params(4);
T(21) = (1-y(5))^params(6);
T(22) = T(8)*T(21);
T(23) = y(3)^params(3);
T(24) = y(5)^(-params(3));

end
