import numpy as np
import scipy as sp
import scipy.linalg as lg

def function_values_array( n, delta, a ): # Creating the matrix f i.e. column vector of size (N+1)X(1)
    val = [] # The column vector stored as row vector
    for ii in range(n+1):
        val = val +[function(a+ii*delta)]
    val = np.array(val)
    return val.transpose() # Making the row vector to column vector

def function( t ): # Function definition
    return t # For now f(x) = x as an sample example.

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
A_cap = coefficient_matrix(m, n, delta, a, b, gamma)
g_cap = function_values_array(n,delta, a)

#print(A_cap.shape)
#print(g_cap.shape)
#print(lg.svdvals(A_cap))

#u , s, v = lg.svd(A_cap)

#print("Coefficient matrix is\n", A_cap)
#print("matrix U is \n",u)
#print("matrix S is \n", s)
#print("matrix V is \n", v)

f_cap = lg.solve(A_cap, g_cap)
print("The solution is:\n", f_cap)
#print(np.allclose(np.dot(A_cap, f_cap), g_cap))
