import scipy.optimize as opt
from numpy import exp

s = 2
a = 0.36
b = 0.98
d = 0.025

def f(x):
    return x**a+(1-d)*x

def du(x):
    return x**(-s)

def df(x):
    return a*(x**(a-1))+(1-d)*x

def e(v) :
    list_of_eq = []
    for i in range(0, 101):
        list_of_eq.append(du(f(v[i]-v[i+1]))/(b*du(f(v[i+1])-v[i+2])) - df(v[i+1]))
    list_of_eq.append(v[101])
    return list_of_eq

solution = opt.fsolve(e,[0]*101) # fsolve(equations, X_0)
print(solution)
