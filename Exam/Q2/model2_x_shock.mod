var k, c, x, z;
varexo ex,ex;

parameters beta, alpha, delta, rho, gamma, sigma, z_s, k_s, c_s, l_s;

alpha = 0.36;
delta = 0.025;
rho = 0.95;
beta = 0.98;
sigma = 2;

z_s=1;
x_s=1;
k_s = ((1/beta-1+delta)/alpha)^(1/(alpha-1));
c_s = k_s^alpha-delta*k_s;


model;
c+x*k = z*(k(-1)^alpha) + x*(1-delta)*k(-1);
c^(-sigma)= beta*(c(+1)^(-sigma))*(alpha*z(+1)*k^(alpha-1)+x(+1)*(1-delta));
log(x) = rho*log(x(-1)) + ex;
log(z) = rho*log(z(-1)) + ez;
end;

initval;
z=1;
x=1
k=k_s;
c=c_s;
l=l_s;
end;

steady(solve_algo=1);

check;

shocks;
var ex = 0.01^2;
end;

stoch_simul(irf=100, order=1, pruning);

%rplot k;

