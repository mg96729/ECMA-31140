var k, c, z;
varexo e;

parameters beta, sigma, alpha, delta, rho,z_s k_s, c_s, z_0;

beta = 0.95;
sigma = 2;
alpha = 0.36;
delta = 0.025;
z_s = 1;
rho = 0.95;

k_s = ((z_s*alpha)^-1)*(1/beta-(1-delta))^(1/(alpha-1));
c_s = z_s*k_s^alpha-delta*k_s;

z_0 = 0.1*z_s;

model;
c = z*k(-1)^alpha + (1-delta)*k(-1)-k;
1 = beta*(c/c(+1))^sigma*(alpha*z(+1)*k^(alpha-1)+1-delta);
log(z) = rho*log(z(-1)) + e;
end;

initval;
z=1;
k=k_s;
c=c_s;
end;

steady(solve_algo=2);

check;

shocks;
var e = 0.1^2;
end;

stoch_simul(irf=100, order=1, pruning);

%rplot k;

