var k, c;

parameters alpha, beta, delta, sigma, k_s, c_s, k_0;

alpha = 0.36;
beta  = 0.95;
delta = 0.025;
sigma = 2;

k_s = ((alpha*beta)/(1+beta*(delta-1)))^(1/(1-alpha));
c_s = k_s^alpha-delta*k_s;

k_0 = 0.1*k_s;

model;
k + c = (k(-1)^alpha) + (1-delta)*k(-1);
c^(-sigma) = beta*((c(+1))^(-sigma))*(alpha*(k)^(alpha-1)+(1-delta));
end;

initval;
k = k_0;
end;

endval;
k = k_s;
c = c_s;
end;

resid;

perfect_foresight_setup(periods = 200);
perfect_foresight_solver;

rplot k;

