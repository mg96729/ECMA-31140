#Code modified from https://macroeconomics.github.io/dynamic-programming-in-python.html

from __future__ import division
import pylab
import matplotlib.pyplot as plt
from numpy import interp
import numpy as np

from scipy.optimize import fminbound

class LinInterp:
    "Provides linear interpolation in one dimension."

    def __init__(self, X, Y):
        """Parameters: X and Y are sequences or arrays
        containing the (x,y) interpolation points.
        """
        self.X, self.Y = X, Y

    def __call__(self, z):
        """Parameters: z is a number, sequence or array.
        This method makes an instance f of LinInterp callable,
        so f(z) returns the interpolation value(s) at z.
        """
        if isinstance(z, int) or isinstance(z, float):
            return interp ([z], self.X, self.Y)[0]
        else:
            return interp(z, self.X, self.Y)


def U(c,sigma=1):
    '''This function returns the value of utility when the CRRA
    coefficient is sigma. I.e.
    u(c,sigma)=(c**(1-sigma)-1)/(1-sigma) if sigma!=1
    and
    u(c,sigma)=ln(c) if sigma==1
    Usage: u(c,sigma)
    '''
    if sigma!=1:
        u=(c**(1-sigma)-1)/(1-sigma)
    else:
        u=np.log(c)
    return u

def F(K, r=0.05):
    '''
    Cobb-Douglas production function
    F(K,L)=K^alpha L^(1-alpha)
    '''
    return (1+r)*K

def Va(k,beta=.95):
    B=1/(1-beta)
    return B*np.log(k)

beta=.95
sigma=1
delta=1


# Grid of values for state variable over which function will be approximated
gridmin, gridmax, gridsize = 0.1, 15, 1000
grid = np.linspace(gridmin, gridmax**1e-1, gridsize)**10

# Maximize function V on interval [a,b]
def maximum(V, a, b):
    return float(V(fminbound(lambda x: -V(x), a, b)))
# Return Maximizer of function V on interval [a,b]
def maximizer(V, a, b):
    return float(fminbound(lambda x: -V(x), a, b))

# The following two functions are used to find the optimal policy and value functions using value function iteration
# Bellman Operator
def bellman(w):
    """The approximate Bellman operator.
    Parameters: w is a LinInterp object (i.e., a
    callable object which acts pointwise on arrays).
    Returns: An instance of LinInterp that represents the optimal operator.
    w is a function defined on the state space.
    """
    vals = []
    for k in grid:
        kmax=F(k)
        h = lambda kp: U(kmax + (1-delta) * k - kp,sigma) + beta * w(kp)
        vals.append(maximum(h, 0, kmax))
    return LinInterp(grid, vals)

# Optimal policy
def policy(w):
    """
    For each function w, policy(w) returns the function that maximizes the
    RHS of the Bellman operator.
    Replace w for the Value function to get optimal policy.
    The approximate optimal policy operator w-greedy (See Stachurski (2009)).
    Parameters: w is a LinInterp object (i.e., a
    callable object which acts pointwise on arrays).
    Returns: An instance of LinInterp that captures the optimal policy.
    """
    vals = []
    for k in grid:
        kmax=F(k)
        h = lambda kp: U(kmax + (1-delta) * k - kp,sigma) + beta * w(kp)
        vals.append(maximizer(h, 0, kmax))
    return LinInterp(grid, vals)

V0=LinInterp(grid,U(grid))

fig, ax = plt.subplots()
ax.set_xlim(grid.min(), grid.max())
ax.plot(grid,Va(grid), label='Actual', color='k', lw=2, alpha=0.6);

count=0
maxiter=200
tol=1e-6
while count<maxiter:
    V1=bellman(V0)
    err=np.max(np.abs(np.array(V1(grid))-np.array(V0(grid))))
    if np.mod(count,10)==0:
        ax.plot(grid,V1(grid), color=plt.cm.jet(count / maxiter), lw=2, alpha=0.6);
        #print '%d %2.10f ' % (count,err)
    V0=V1
    count+=1
    if err<tol:
        print(count)
        break
ax.plot(grid,V1(grid), label='Estimated', color='r', lw=2, alpha=0.6);
ax.legend(loc='lower right')
plt.draw();

plt.show()





