# Global-Fourier-Continuation
We present a method for global fourier continuation based on O.P. Bruno Paper in 2002 given https://doi.org/10.1006/jcph.2002.7023 <br/>
We use singular value decomposition(SVD) library from numpy to solve for fourier Coefficeints.

# Note:
The SVD_python_trial.py is the trail version of the code for the continuation which has been successfully corrected.<br/>
The corrected file name: SVD_python_code.py

# These are the following desired result:
The code takes the following inputs:
1. The left and right end points of the inetrval in which the function is known analytically.
2. The largest value of M in the truncated complex Fourier series i.e. if the series is represented as <br/>
<a href="https://www.codecogs.com/eqnedit.php?latex=g(z)=\sum_{k=-M}^M&space;c_je^{\frac{2\pi&space;ik(z-a)}{L&plus;d}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g(z)=\sum_{k=-M}^M&space;c_je^{\frac{2\pi&space;ik(z-a)}{L&plus;d}}" title="g(z)=\sum_{k=-M}^M c_je^{\frac{2\pi ik(z-a)}{L+d}}" /></a> .
3. The number of points of Evalution N of the truncated complex Fourier series. It should be greater than or equal to 2 * M.
4. The interval of extension i.e. the interval on which we want to extend.<br/>
<a href="https://www.codecogs.com/eqnedit.php?latex=\triangle\delta&space;=\frac{L}{N}=\frac{b-a}{N}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\triangle\delta&space;=\frac{L}{N}=\frac{b-a}{N}" title="\triangle\delta =\frac{L}{N}=\frac{b-a}{N}" /></a> , <a href="https://www.codecogs.com/eqnedit.php?latex=d=\triangle\delta\times\gamma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?d=\triangle\delta\times\gamma" title="d=\triangle\delta\times\gamma" /></a>

The results produced:
1. The graph of the Fourier extension
2. The conditional number for the linear system.
