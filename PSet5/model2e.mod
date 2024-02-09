var k, c, l, z;
varexo e;

parameters sigma, beta, alpha, delta, rho, gamma, z_s, k_s, c_s, l_s;

alpha = 0.36;
delta = 0.025;
rho = 0.95;
beta = 0.98;

load params_model;
set_param_value('sigma',sigma);

z_s=1;
l_s = 1/3;
k_s = ((1/beta-1+delta)/(alpha*l_s^alpha))^(1/(alpha-1));
c_s = k_s^alpha*l_s^(1-alpha)-delta*k_s;

gamma = (c_s*l_s^alpha)/(k_s^alpha - alpha*k_s^alpha + c_s*l_s^alpha - k_s^alpha*l_s + alpha*k_s^alpha*l_s);


model;
c+k = z*(k(-1)^alpha)*(l^(1-alpha)) + (1-delta)*k(-1);
(c^(gamma-1)*gamma*(1-l)^(1-gamma))/(c^gamma*(1-l)^(1-gamma))^sigma = beta*((c(+1)^(gamma-1)*gamma*(1-l(+1))^(1-gamma))/(c(+1)^gamma*(1-l(+1))^(1-gamma))^sigma)*(alpha*z(+1)*k^(alpha-1)*(l(+1)^(1-alpha))+1-delta);
0=(c^gamma*(gamma-1))/((c^gamma*(1-l)^(1-gamma))^sigma*(1-l)^gamma)+(c^(gamma-1)*gamma*(1-l)^(1-gamma))/(c^gamma*(1-l)^(1-gamma))^sigma*z*(1-alpha)*(k^alpha)*(l^(-alpha));
log(z) = rho*log(z(-1)) + e;
end;

initval;
z=1;
k=k_s;
c=c_s;
l=l_s;
end;

shocks;
var e = 0.1^2;
end;

stoch_simul(irf=100, order=1, pruning);

%rplot k;

