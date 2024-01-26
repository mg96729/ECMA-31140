var k, c, l;

parameters alpha, beta, delta, gamma, k_s, c_s, l_s, k_0;

alpha = 0.36;
beta  = 0.95;
delta = 0.025;

l_s = 1/3;
k_s = (l_s)*((1/alpha)*(1/beta-1+delta))^(1/(alpha-1));
c_s = (k_s^alpha)*(l_s^(1-alpha))-delta*k_s;
gamma = (1/c_s)*(1-alpha)*k_s^alpha*l_s^-alpha;

k_0 = 0.1*k_s;

model;
k + c = (l^(1-alpha))*(k(-1)^alpha) + (1-delta)*k(-1);
1/c = beta*(1/(c(+1)))*(alpha*l(+1)^(1-alpha)*k^(alpha-1)+1-delta);
gamma=(1/c)*(1-alpha)*(k(-1))^alpha*l^-alpha;
end;

initval;
k = k_0;
end;

endval;
k = k_s;
c = c_s;
l = l_s;
end;

%resid;

perfect_foresight_setup(periods = 200);
perfect_foresight_solver;

%rplot k;


