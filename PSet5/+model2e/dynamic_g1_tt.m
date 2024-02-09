function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 30);

T = model2e.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(25) = getPowerDeriv(y(4),params(6),1);
T(26) = getPowerDeriv(T(5)*T(7),params(1),1);
T(27) = (T(8)*T(5)*params(6)*getPowerDeriv(y(4),params(6)-1,1)-T(6)*T(5)*T(25)*T(26))/(T(8)*T(8));
T(28) = getPowerDeriv(T(14),params(1),1);
T(29) = (-(getPowerDeriv(1-y(5),1-params(6),1)));
T(30) = (-(getPowerDeriv(1-y(8),1-params(6),1)));

end
