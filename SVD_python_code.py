import numpy as np
import scipy as sp
import scipy.linalg as lg
from matplotlib import pyplot as plt
from matplotlib import style

style.use('ggplot')

def function_values_array( n, delta, a ): # Creating the matrix f i.e. column vector of size (N+1)X(1)
    val = [] # The column vector stored as row vector
    for ii in range(n+1):
        val = val +[function(a+ii*delta)]
    val = np.array(val)
    return val.transpose() # Making the row vector to column vector

def function_lyon( t ): # Function definition
    return np.exp(np.sin(5.4*np.pi*t-2.7*np.pi)-np.cos(2*np.pi*t)) # For now f(x) = x as an sample example.

def function(t):
    return t

def function_val_ploting(t, fourier_coefficients, extended_interval, interval, m, a):
    ld = interval+extended_interval
    period_number = t//ld
    t = t-ld*period_number
    g_z = 0
    for ii in range(2*m+1):
        g_z += fourier_coefficients[ii]*np.exp(2*np.pi*1j*(ii-m)*(t-a)/ld)
    return  g_z

def coefficient_matrix(m, n , delta, a, b, gamma): # Creating the matrix A of size (2M+1)X(N+1)
    l = b-a # Length of the interval
    d = gamma*delta # Length of Fourier extension
    ld = l + d # Length of the period after extension
    coefficient_matrix_array = [] # The coefficient Matrix A
    for ii in range(2*m+1):
        row_array = [] # Temporary matrix, a row vector
        for jj in range(n+1):
            row_array = row_array+[np.exp(2*np.pi*1j*(ii-m)*jj*delta/ld)]
        coefficient_matrix_array = coefficient_matrix_array+[row_array]
    return (np.array(coefficient_matrix_array)).transpose() # Arranges the coefficient_matrix_array in actual matrix form

a = float(input("Enter the left end-point of the interval:"))
b = float(input("Enter the right end-point of the interval:"))
m = int(input("Enter the value of M:"))
n = int(input("Enter the value of N:"))
gamma = int(input("Enter the value of gamma:"))
delta = (b-a)/n
extended_interval = gamma*delta
interval = b-a

A_cap = coefficient_matrix(m, n, delta, a, b, gamma)
g_cap = function_values_array(n,delta, a)


#u , s, v = lg.svd(A_cap)

A_cap_inv = lg.pinv(A_cap)
f_cap = np.matmul(A_cap_inv, g_cap)

x = []
x_f = []
y = []
alpha =4
y_f = []
for ii in range(1000):
    x_prime = a + ii*alpha/1000
    x += [x_prime]
    ld = interval + extended_interval
    y_prime = function_val_ploting(x_prime, f_cap, extended_interval, interval, m, a)
    y += [y_prime]
    if x_prime <= (a+interval):
        x_f += [x_prime]
        y_f += [function(x_prime)]


plt.plot(x,y,'g', label='Fourier continuation')
plt.plot(x_f,y_f, 'c', label = 'f(x)')
plt.title("Global Fourier Extension of f(x)=x")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.legend()
plt.show()

print("The condition number for the coefficient matrix is ", lg.norm(A_cap_inv)*lg.norm(A_cap))
