import numpy as np
import math as ma
import scipy as sp
import scipy.linalg as lg
from matplotlib import pyplot as plt
from matplotlib import style

style.use('ggplot')


def function_values_array(n, delta, a):  # Creating the matrix f i.e. column vector of size (N+1)X(1)
    val = []  # The column vector stored as row vector
    for ii in range(n+1):
        val = val + [function_lyon(a+ii*delta)]
    val = np.array(val)
    return val.transpose()  # Making the row vector to column vector


def function_lyon(t):  # Function definition
    return np.exp(np.sin(5.4*np.pi*t-2.7*np.pi)-np.cos(2*np.pi*t))  # For now f(x) = x as an sample example.


def function(t):
    return t


def function_val_ploting(t, fourier_coefficients, extended_interval, interval, m, a):
    ld = interval+extended_interval
    period_number = t//ld
    t = t-ld*period_number
    g_z = fourier_coefficients[0]
    for ii in range(m):
        g_z += fourier_coefficients[2*ii+1]*np.cos(2*np.pi*(ii+1)*(t-a)/ld)
        g_z += fourier_coefficients[2*(ii+1)]*np.sin(2*np.pi*(ii+1)*(t-a)/ld)
    return g_z


def coefficient_matrix(m, n, delta, a, b, gamma):  # Creating the matrix A of size (2M+1)X(N+1)
    l = b-a  # Length of the interval
    d = gamma*delta  # Length of Fourier extension
    ld = l + d  # Length of the period after extension
    coefficient_matrix_array = []  # The coefficient Matrix A
    for ii in range(2*m+1):
        row_array = []  # Temporary matrix, a row vector
        for jj in range(n+1):
            row_array = row_array+[np.exp(2*np.pi*1j*(ii-m)*jj*delta/ld)]
        coefficient_matrix_array = coefficient_matrix_array+[row_array]
    # Arranges the coefficient_matrix_array in actual matrix form
    return (np.array(coefficient_matrix_array)).transpose()


def fourier_complex_to_real(f_cap, m):
    f_cap_real_dash = []
    if abs(f_cap[m].real) < 10**(-12):  # Setting degree of precision to be 10^(-11)
        f_cap_real_dash += [0.0]
    else:  # Truncating the number after 11 digits
        f_cap_real_dash += [ma.trunc(f_cap[m].real * (10 ** 11)) / (10 ** 11)]

    for alp in range(m):
        if abs(f_cap[m+1+alp].real) < 10**(-12):
            f_cap_real_dash += [0.0]
        else:
            f_cap_real_dash += [2*ma.trunc(f_cap[m+1+alp].real*(10**11))/(10**11)]
        if abs(f_cap[m+1+alp].imag) < 10**(-12):
            f_cap_real_dash += [0.0]
        else:
            f_cap_real_dash += [-2*ma.trunc(f_cap[m+1+alp].imag*(10**11))/(10**11)]
    f_cap_real_dash = np.array(f_cap_real_dash)
    return f_cap_real_dash.transpose()


# Input left end-point of the interval
a = float(input("Enter the left end-point of the interval:"))
# Input right end-point of the interval
b = float(input("Enter the right end-point of the interval:"))
# The number of terms in the truncated fourier series is 2*m+1
m = int(input("Enter the value of M:"))
# The number of points interpolation points chosen in the interval [a,b]. It should be at least 2m+1 i.e. n>=2m
n = int(input("Enter the value of N:"))
# The gamma is the length of the extension. A 30% of the interval length is sufficient. n>=gamma>=0.3*n
gamma = int(input("Enter the value of gamma:"))
# We calculate the length of the smallest piece of the interval after dividing it in equal n parts.
delta = (b-a)/n
"""This is the length of the extended interval i.e. the given non-periodic function becomes periodic after we
extended its interval to [a, b+extended_interval]. The period of the function becomes b+extended_interval-a."""
extended_interval = gamma*delta
# The length of the interval on which the non-periodic function values are exactly known.
interval = b-a

# It is the Fourier coefficient matrix of sin and cos formed by evaluation of the fourier series at n different points
A_cap = coefficient_matrix(m, n, delta, a, b, gamma)
# It is the data vector i.e. the values of the function at n+1 different points.
g_cap = function_values_array(n, delta, a)

# The svd decomposition of the Fourier coefficient matrix of sin and cos
# u , s, v = lg.svd(A_cap)

# It is the inverse of the Fourier coefficient matrix of sin and cos obtained using SVD.
A_cap_inv = lg.pinv(A_cap)
# It is the matrix containing the values of unknown complex Fourier coefficient obtained by solving linear equation AF=G
f_cap = np.matmul(A_cap_inv, g_cap)
# It is the matrix containing the values of unknown real Fourier coefficients obtained from the complex coefficients.
f_cap_real = fourier_complex_to_real(f_cap, m)

# This will contain x-axis points for plotting Global fourier continuation graph
x = []
# This will contain x-axis points for plotting function values
x_f = []
# This contains global fourier continuation function value corresponding to the point on x-axis
y = []
# The interval to interpolate i.e. from [a, a+alpha]
alpha = 4
# This contains the function values corresponding to the values of the point on the x-axis.
y_f = []

for ii in range(1000):
    # We are taking 1000 points to plot the graph accurately on the graph length 4.
    x_prime = a + ii*alpha/1000
    x += [x_prime]
    ld = interval + extended_interval  # Value of the new period.
    y_prime = function_val_ploting(x_prime, f_cap_real, extended_interval, interval, m, a)
    y += [y_prime]
    if x_prime <= (a+interval):
        x_f += [x_prime]
        y_f += [function_lyon(x_prime)]


# This will plot the Graph of the Global Fourier continuation
plt.plot(x, y, 'black', linestyle='dashed', label='Fourier continuation')
# This will plot the Graph of the function
plt.plot(x_f, y_f, 'red', label='f(x)')
# This will add title to the graph
plt.title("Global Fourier Extension of f(x)=x")
# These will add labels to the different axis in the graph
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
# This will create legends in the graph at appropriate place.
plt.legend()
# This will show the graph in a separate window.
plt.show()

print("The condition number for the coefficient matrix is ", lg.norm(A_cap_inv)*lg.norm(A_cap))
