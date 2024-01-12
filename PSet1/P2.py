from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()

# syntax for 3-D projection
ax = plt.axes(projection='3d')

r = 0.95**(-1)-1
# defining all 3 axis: x is sigma_a, y is sigma_b, z is b_a=b_b
x = np.linspace(0.1, 2, 100)
y = np.linspace(0.1, 2, 100)
z = (1-(0.95**(-(1/x)))*(4+2*r)*(1+r)**(-(1/x)))/(1-(0.95**(-(1/x)))*(1+r)*(1+r)**(-(1/x)))

# plotting
ax.plot3D(x, y, z, 'green')
ax.set_title('3D plot')
plt.show()