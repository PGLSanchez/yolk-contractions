from __future__ import division,print_function
import os,sys
from pylab import *
from math import atan2
from os import path,walk
from scipy.optimize import leastsq,minimize
from scipy.signal import hilbert,cwt,ricker,lombscargle,welch,morlet
from matplotlib import rc
from xlrd import *

from numpy import linspace as lp


#=======================Autoregressive Process=============================================


def ar1_sim(alpha,sigma,N,x0 = None):

    N = int(N)
    sol = zeros(N)

    if x0 is None:
        x0 = randn()
        
    sol[0] = x0

    for i in range(1,N):
        sol[i] = alpha*sol[i-1] + sigma*randn()

    return sol


def ar1_log_likelihood(par,sol):

    alpha = par[0]
    sigma = par[1]

    N = len(sol)
    sum_ = sum((sol[1:] - alpha*sol[:-1])**2)

    logL = 0.5 * log( (1-alpha**2)/sigma**2 ) - (1-alpha**2)/(2*sigma**2)*sol[0]**2 - (N-1)/2.*log(sigma**2) - 1/(2*sigma**2) * sum_

    # print alpha,sigma,-logL
    return -logL


def ar1_MLE(sol, bounds =  ( (1e-5,1),(1e-5,10) ) ,full_output = False):

    sol = sol/std(sol)

    ''' Estimate AR1 parameters    '''
    res = minimize(ar1_log_likelihood, x0 = array((0.3,0.4)),args = (sol,),method ='SLSQP',bounds = bounds)

    if full_output:
        return res
    
    return res['x']

def ar1_powerspec(alpha,periods,dt):
    res = (1-alpha**2)/(1+alpha**2 - 2*alpha*cos(2*pi*dt/periods))

    return res

#=====================================================================================

def smooth(x,window_len=11,window='bartlett',data = None):
    """smooth the data using a window with requested size.

    input:
    x: the input signal
    window_len: the dimension of the smoothing window; should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    flat window will produce a moving average smoothing.
    data: if not None, will be used as evaluated window!

    """

    x = array(x)

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        raise ValueError("window must not be shorter than 3")

    if window_len%2 is 0:
        raise ValueError("window_len should be odd")

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman','triang']:
       raise ValueError("Window is none of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman','triang'")

    # use externally derieved window evaluation
    if data is not None:
        window_len = len(data)
        window = 'extern'
   
    s=r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
                                        #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')

    elif window == 'triang':
        w = triang(window_len)

    elif window == 'extern':
        w = data
        
    else:
        w=eval(window+'(window_len)')

    y=convolve(w/w.sum(),s,mode='valid')
    
    return y[(window_len-1)/2:len(y)-(window_len-1)/2]

def sinc_filter(M, f_c = 0.2):

    ''' 
    f_c in sampling frequency unit, max 0.5! 
    M is blackman window length and even!

    '''

    # not very effective, but should be get called only once per convolution

    assert M%2 == 0,'M must be even!'
    res = []

    for x in arange(0,M+1):
            
        if x == M/2:
            res.append(2*pi*f_c)
            continue
    
        r = sin(2*pi*f_c*(x - M/2))/( x - M/2 ) # the sinc filter unwindowed
        r = r * (0.42 - 0.5*cos(2*pi*x/M) + 0.08*cos(4*pi*x/M)) # blackman window
        res.append(r)

    res = array(res)
    res = res/sum(res)
            
    return res


def sinc_detrend(raw_signal,T_c,dt):

    signal = array(raw_signal)
    dt = float(dt)

    # relative cut_off frequency
    f_c = dt/T_c
    M = len(signal) # max for sharp roll-off

    # M needs to be even
    if M%2 != 0:
        M = M - 1

    w = sinc_filter(M, f_c)  # the evaluated windowed sinc filter
    sinc_smoothed = smooth(signal, data = w)
    sinc_detrended = signal - sinc_smoothed

    return sinc_detrended

def detrend(raw_signal,window = 'bartlett',winsize = 7, data = None):

    avsignal = smooth(raw_signal,winsize,window = 'flat', data = data) 
    dsignal = raw_signal - avsignal             # detrend by subtracting filter convolution

    return dsignal

# acorr still better
def acf(x,maxlag = 10):

    res = [1] + [corrcoef(x[:-i],x[i:])[0,1] for i in range(1,maxlag)]
    return array(res)


#================Generic Fitting========================

def FitFunc(func,yvec,ipars,xvec = None):

    if xvec is None:
        xvec = arange(len(yvec))
        
    pfit , success = leastsq(errfunc, ipars, args = (xvec,yvec,func) )
    
    return pfit

def errfunc(p,xvec,yvec,tfunc):
    return tfunc(p,xvec) - yvec

def expfunc(p,x):
    return p[0]*exp(-p[1]*x)
