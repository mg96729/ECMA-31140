var k, c, l, z;
varexo e;

parameters beta, alpha, delta, rho, gamma, z_s, k_s, c_s, l_s;

alpha = 0.36;
delta = 0.025;
rho = 0.95;
beta = 0.98;

%we use the gamma found in the previous part:
gamma = 0.38514680483592400690846286701209;

z_s=1;
l_s = 1/3;
k_s = ((1/beta-1+delta)/(alpha*l_s^alpha))^(1/(alpha-1));
c_s = 0.8564;


model;
c+k = z*(k(-1)^alpha)*(l^(1-alpha)) + (1-delta)*k(-1);
gamma/c = beta*(gamma/c(+1))*(alpha*z(+1)*k^(alpha-1)*(l(+1)^(1-alpha))+1-delta);
0=(1-gamma)/(l-1) + (gamma/c)*z*(1-alpha)*(k^alpha)*(l^(-alpha));
log(z) = rho*log(z(-1)) + e;
end;

initval;
z=1;
k=k_s;
c=c_s;
l=l_s;
end;

steady(solve_algo=1);

check;

shocks;
var e = 0.1^2;
end;

stoch_simul(irf=100, order=1, pruning);

%rplot k;

