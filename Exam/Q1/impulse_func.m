function [C,K] = impulse_funcs(sig)

beta = 0.98;
sigma = sig;
alpha = 0.36;
delta = 0.025;
z = 1;
rho = 0.95;

syms c k;

u = (c^(1-sigma))/(1-sigma);
f = z*(k^alpha)+(1-delta)*k;

%the derivatives we need
uc = diff(u,c);
ucc = diff(u,c,2);
fz = k^alpha;
fkz = alpha*(k^(alpha-1));
fkk = diff(f,k,2);

%Steady State Values:
eqn = diff(f,k) == 1/beta;
k_s = vpasolve(eqn,k);
c_s = k_s^alpha - delta*k_s;

%Now substitute the ss values into the derivatives:
u_c = subs(uc,c,c_s);
u_cc = subs(ucc,c,c_s);
f_z = subs(fz,k,k_s);
f_kz = subs(fkz,k,k_s);
f_kk = subs(fkk,k,k_s);

%Now define matrix A:
A = [1/beta -1;-u_c*f_kk/u_cc 1+beta*u_c*f_kk/u_cc];

%Now find Eigenvalues and vectors of A:
[V,D] = eig(A);

%Sort by eigenvalue ascending:
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
V = V(:,ind);

V_l = inv(V)

B = [f_z; -beta*u_c*(f_kk*f_z+f_kz*rho)/u_cc];

C = V_l*B;

a_1 = -V_l(2,1)/V_l(2,2);
a_2 = -C(2,1)/((Ds(2,2)-rho)*V_l(2,2));

b_1 = 1/beta+V_l(2,1)/V_l(2,2);
b_2 = f_z+C(2,1)/((Ds(2,2)-rho)*V_l(2,2));

syms k_t z_t
C = a_1*k_t+a_2*z_t;
K = b_1*k_t+b_2*z_t;