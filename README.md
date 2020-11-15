# Global-Fourier-Continuation
We present a method for global fourier continuation based on O.P. Bruno Paper in 2002 given https://doi.org/10.1006/jcph.2002.7023

#We use singular value decomposition(SVD) library from numpy to solve for fourier Coefficeints.


The SVD_python_trial.py is the trail version of the code for the continuation which has been successfully corrected.
The corrected file name: SVD_python_code.py

# These are the following desired result:
The code takes the following inputs:
1. The left and right end points of the inetrval in which the function is known analytically.
2. The largest value of M in the truncated complex fourier series i.e. if the series is represented as <a href="https://www.codecogs.com/eqnedit.php?latex=g(z)=\sum_{k=-M}^M&space;c_je^{\frac{2\pi&space;ik(z-a)}{L&plus;d}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g(z)=\sum_{k=-M}^M&space;c_je^{\frac{2\pi&space;ik(z-a)}{L&plus;d}}" title="g(z)=\sum_{k=-M}^M c_je^{\frac{2\pi ik(z-a)}{L+d}}" /></a> .
