from scipy.misc import derivative
import math


#first, we define our function, as well as initial values
# If we are to look at a different function, the definition needs to be changed
def func(x):
    return (x-10)* math.exp(-x**2)+5

#the initial values 0 and 10 suffice, as f(0) = -5, f(10)=5
#Admittedly this is a sloppy way of finding such values.
# A more general way may take more computational time though

#This is for bisection
def bisection(func, a, b, e, d):
    if func(a) * func(b) >= 0:
        raise BaseException("Wrong interval!")
    if func(a) > func(b):
        a, b = b, a

    m = (a + b) / 2
    while abs(b - a) > e:
        if abs(func(m)) < d:
            return m
        elif func(a) * func(m) < 0:
            b = m
        else:
            a = m
        m = (a + b) / 2

    return m

print(bisection(func, 0, 10, 0.00001, 0.00001))
print(func(0.7821178436279297))


d = derivative(func, 1.0, dx=1e-3)

#part 2: Newton's Method
def newton(func, x0, e, de):
    xt = x0
    d = derivative(func, xt, dx=1e-3)
    #note that we are using numerical estimates of deriv itself.
    # This allows numpy to handle less conventional cases, such as this one.

    xt1 = xt - func(xt)/d
    while abs(xt-xt1) > e * (1+abs(xt1)):
        if d ==0:
            return "Derivative reaches 0 before solution found"
        else:
            xt = xt1
            d = derivative(func, xt, dx = 1e-3)
            xt1 = xt - func(xt)/d

    if func(xt1) > de:
        return "result not found"
    else:
        return "value found: "+ str(xt1)


print(newton(func, 1, 1e-3, 1e-3))

print([0]*5)