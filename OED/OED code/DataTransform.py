# -*- coding: utf-8 -*-
"""
Created on Wed Feb 06 10:12:02 2019

@author: huckg
"""
#from __future__ import division
import numpy
import copy
import itertools

''' function that takes a list of dicts and turns it into Dict of lists'''
def LDtoDL(LD):
    nd={}
    for d in LD:
        for k,v in d.items():
            try:
                nd[k].append(v)
            except KeyError:
                nd[k]=[v] 
    return nd

def DLtoLD(DL):
    return [dict(zip(DL,t)) for t in zip(*DL.values())]
        
def unpack(parameters,scores,observables,criteria = [],binary = []):
    if type(parameters) == list:
        parameters = {i:parameters[i] for i in range(len(parameters))}
    p_unpack = {i:[] for i in list(parameters.values())[-1]}
    for i in range(len(parameters)):
        for parameter,value in parameters[i].items():
            p_unpack[parameter].append(value)
            
#    unpack the scores in a dictionary of ordered lists by state observed       
    s_unpack = {i:{} for i in observables}
    for k,v in s_unpack.items():
        s_unpack[k] = {i:[] for i in criteria}
    s_bin = copy.deepcopy(s_unpack)
    pdict = {i:[] for i in list(p_unpack.keys())}
    p_bin = {i:copy.deepcopy(pdict) for i in observables}
    for i in range(len(scores)):
        for state in observables:
            for criterion,value in scores[i][state].items():
                s_unpack[state][criterion].append(value)
    if binary:
        for state in observables:
            for bit in binary:
                for i in range(len(s_unpack[state][bit])):
                    if s_unpack[state][bit][i]:
                        for criterion in criteria:
                            s_bin[state][criterion].append(s_unpack[state][criterion][i])
                        for parameter,values in p_unpack.items():
                            p_bin[state][parameter].append(values[i])
                    if not s_unpack[state][bit][i]:
                        for criterion in criteria:
                            s_bin[state][criterion].append(s_unpack[state][criterion][i])
                        for parameter,values in p_unpack.items():
                            p_bin[state][parameter].append(values[i])                    
    return p_unpack,s_unpack,p_bin,s_bin
                

def quantiles(scrs,prmt,observables,criteria,binary,fraction = 0.1):
    '''get the quantiles which seperate the scores in their extremes to find
    parameter sets which are grouped together'''
    lower,upper = ({},{})
    for state in observables:
        lower[state] = {}
        upper[state] = {}
        for cr in criteria:
            '''rank the elements in place of the list by creating a list that ranks
            the elements from low to high''' 
            ranked = sorted(range(len(scrs[state][cr])),key=scrs[state][cr].__getitem__)
            '''list indices of the elements that were ranked in place and sort them 
            according to an upper and lower, its total number limited by
             the fraction ''' 
            ilow  = [i for i in range(len(ranked)) if ranked[i] < int(len(ranked)*fraction)]
            ihigh = [i for i in range(len(ranked)) if ranked[i] > int(len(ranked)*(1-fraction))]  

            '''the dictionaries with the parameter values needed for ''' 
            lower[state][cr] = {p:[] for p in prmt.keys()}
            upper[state][cr] = {p:[] for p in prmt.keys()}
            for k in prmt.keys():
                for i in ilow:
                    lower[state][cr][k].append(prmt[k][i])                
                for i in ihigh:
                    upper[state][cr][k].append(prmt[k][i])       
    return [lower,upper]


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    #taken from http://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html, not coded by @bobvansluijs
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data. It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    """
    import numpy as np
    from math import factorial
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    
    # pad the signal at the extremes with values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def interpolate(t_vector,s_time,ts):
    inrtpl_data = numpy.interp(t_vector, s_time,ts)
    sid = savitzky_golay(inrtpl_data,101,3)
    return sid

    
def interpolate_dataframe(tv,dv,tu,desired = "min"):
    units  = ["sec","min","hour","day"];factors = [1,60,3600,86400]    
    permutations = [i for i in itertools.product(units,repeat = 2)]
    
    for i in range(len(permutations)):
        start,end = permutations[i]
        istart = units.index(start)
        iend   = units.index(end)
        if istart == iend:
            permutations[i] += (1,)
        else:
            fct = factors[istart]/factors[iend]
            permutations[i] += (fct,)
    for i in permutations: 
        initial,des,factor = i
        if (initial,des) == (tu,desired):
            conversion = factor
        else:
            conversion = 1
    c = 0      
    for i in dv:
        if i > 10**100:
            dv[c] = dv[c-1]
        c += 1
        
    intp = numpy.interp(numpy.linspace(0,int(conversion*tv[-1]),int(conversion*tv[-1])),tv*conversion,dv)
    profile = savitzky_golay(intp,11,3)
    time = range(0,int(conversion*int(tv[-1])),1)
    return profile,time,conversion

def powerspectrum(data,dt = 1):
    #these functions are meant to transform datasets to something more manageable e.g. a fourier transform
    data       = numpy.array(data)
    normalized = data - numpy.mean(data)
    transform  = numpy.abs(numpy.fft.fft(normalized))**2
    freqs      = numpy.fft.fftfreq(int(data.size),dt)
    arg        = numpy.argsort(freqs)
    data       = transform[arg]
    data       = data[int(len(data)/2):-1]
    return data

