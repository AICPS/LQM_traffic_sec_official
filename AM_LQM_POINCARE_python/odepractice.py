import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint

def dy_dx(y,x):
	return x-y

# approximate plot; numerical solution
xs = np.linspace(0,5,100)
print(xs)
y0 = 1.0
ys = odeint(dy_dx, y0, xs)
ys = np.array(ys).flatten()

plt.rcParams.update({'font.size':14})
plt.xlabel("x")
plt.ylabel("y")
plt.plot(xs,ys);

# exact plot; analytical solution
y_exact = xs - 1 + 2*np.exp(-xs)
y_difference = ys - y_exact
plt.plot(xs, ys, xs, y_exact, "+");
plt.show()
# difference between two series
y_diff = np.abs(y_exact - ys)
plt.semilogy(xs, y_diff)
plt.ylabel("Error")
plt.xlabel("x")
plt.title("Error in numerical integration");

plt.show()