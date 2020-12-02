import matplotlib.pyplot as plt
import numpy as np


def function_tau(x, in_radius, out_radius):
    if abs(x) <= in_radius:
        return 1
    elif in_radius < abs(x) < out_radius:
        radius = (abs(x)-in_radius)/(out_radius-in_radius)
        return np.exp((2*np.exp(-1/radius))/(radius-1))
    else:
        return 0


def non_periodic_function(x):
    return np.exp(np.sin(5.4*np.pi*x-2.7*np.pi)-np.cos(2*np.pi*x))


tau_x = []  # Storing x values for plotting tau function
tau_y = []  # Storing y values for plotting tau function
np_y = []  # Storing the y values for plotting non-periodic function
first_kind_y = []  # Storing the y-values for plotting the graph of First kind Fourier continuation.

in_tau_radius = float(input("enter the inner radius for tau:"))  # Inputs the value of inner radius for tau function
out_tau_radius = float(input("enter the outer radius for tau:"))  # Inputs the value of outer radius for tau function

plot_left_end = -1.5  # Left end point of the graph
plot_right_end = 1.5  # Right end point of the graph
plot_points = 2000  # Number of points used for plotting the graph
physical_left = -1.0  # Left end point of the physical interval.
physical_left_ii = 0
physical_right = 1.0  # Right end point of the physical interval.
physical_right_ii = 0

slice_interval = (plot_right_end-plot_left_end)/plot_points  # The length the function value is unknown

for ii in range(plot_points+1):
    x_prime = plot_left_end + slice_interval*ii
    if abs(x_prime-physical_left) < slice_interval:
        physical_left_ii = ii
    elif abs(x_prime-physical_right) < slice_interval:
        physical_right_ii = ii
    tau_x += [x_prime]
    tau_y += [function_tau(x_prime, in_tau_radius, out_tau_radius)]
    np_y += [non_periodic_function(x_prime)]
    first_kind_y += [(tau_y[ii])*(np_y[ii])]


plt.plot(tau_x, tau_y, 'red', linewidth=1)
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Window function $\mathscr{T}$")
#plt.legend()
plt.show()

plt.plot(tau_x, np_y, 'orange', linestyle='dotted', linewidth=1)
plt.plot(tau_x[physical_left_ii:physical_right_ii+1], np_y[physical_left_ii:physical_right_ii+1], 'orange', linewidth=1)
plt.title("Lyon function")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
# plt.legend()
plt.show()

plt.plot(tau_x, first_kind_y, 'orange', linewidth=1, label='$\hat{f}=\mathscr{T}f$')
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Fourier Continuation of First kind")
plt.legend()
plt.show()
