# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 14:25:47 2019

@author: huckg
"""
from __future__ import division
import numpy
import matplotlib.pylab as plt
from scipy.signal import argrelextrema as agr
from itertools import chain
from Operations import *
from DataTransform import *

def extremasearch(cps):
    mn,mx = (list(agr(cps, numpy.less)[0]),list(agr(cps, numpy.greater)[0]))
    return mn,mx,list(sorted(mn+mx))

evaluate = {}
task = lambda f: evaluate.setdefault(f.__name__, f)

""" this is a list of functions which either fit or extract information from a 1D
numpy array, to be used in generalized evaluation and or assessment functions"""

"""calculates the amplitude only works for stable oscillations use fourier for 
unstable oscillations"""
@task    
def amplitude(data,attr = []):
    if not attr:
        attr = extremasearch(data)
    mnm,mxm,ext = attr
    try:
        scr = abs(data[mnm[-1]] - data[mxm[-1]])
    except IndexError:
        scr = 0
    return  scr


"""calculates period, same conditions as amplitude apply"""
@task
def period(data,attr = []):
    if not attr:
        attr = extremasearch(data)
    mnm,mxm,ext = attr
    try:
        scr = abs(mxm[-1] - mxm[-2])
    except IndexError:
        scr = 0
    return scr


""" gets mean value i.e. around which mean does my system oscillate"""
@task
def mean(data,attr = []):
    mean = numpy.mean(data)
    if  10**-10 < mean < 10**10:
        mean = 0
    return numpy.mean(data)

""" Check if there is a minimum number of minima or maxima to classify oscillations exist"""
@task
def oscillations(data,attr = []):
    if not attr:
        attr = extremasearch(data)
    minima,maxima,extrema = attr   
    if len(maxima) <= 3 or len(minima) <= 2:
        boolean = False
    else:
        boolean = True
    return boolean

""" perform a fourier transform of the data for optimization of fitting of complex oscillations"""
@task
def fouriertransform(data,attr):
    if not attr:
        dt = 1
    else:
        dt = attr
    data = powerspectrum(data, dt=dt)
    return data

""" LSQ fit, note depending on the concentration differential a 
normalized fit would work better i.e. CHI squared,  the absolute difference 
between 2 and 200 is two orders but 2 could deviate 50% from the target whereas 200
only 10% if the concentrations fitted are 4 and 2200"""
@task
def leastsquares(data,attr = []):
    return numpy.sum((numpy.array(data) - numpy.array(attr))**2)


""" LSQ fit of the differential, same applies"""
@task
def leastsquaresdifferential(data,attr = []):
    return numpy.sum((numpy.gradient(numpy.array(data)) - numpy.gradient(numpy.array(attr)))**2)

"""bistability assessment of first, secnd, ration and conditional"""
@task
def firststate(data,attr):
    try:
        length = int(len(data)/10)
        first = data[0:length][-1]
    except:
        first = 0
    return first
@task
def secondstate(data,attr):
    try:
        second = data[-1]
    except:
        second = 0
    return second
@task
def stateratio(data,attr):
    try:
        firstlength  = int(len(data)/10)
        secondlength = data[-1]
        ratio  = float(data[0:firstlength][-1]) / float(data[secondlength:len(data)][-1]) 
    except:
        ratio = 1
    return ratio

@task
def bistable(data,attr):
    try:
        firstlength  = int(len(data)/10)
        ratio  = float(data[0:firstlength][-1]) / float(data[-1]) 
    except:
        ratio = 1
    boolean = True
    if 0.8 <= ratio <= 1.2:
        boolean = False   
    if len(data) != 10000:
        boolean = False
  
    plt.plot(data)
    plt.show()
    return boolean

@task
def steadystates(data,attr = []):
    return data[-1]
    

@task
def normalized_integral(data,attr = []):
    data = numpy.diff(data)
    data = [abs(i) for i in data]
    average = mean(data)
    new  = data/average
    integral = numpy.trapz(new)
    return integral 
    
@task
def extrema_score(data, attr = []):
        """obtain the exrema indices"""
        maxima = agr(data, numpy.greater) 
        minima = agr(data, numpy.less)    
        extrema_chain = list(chain((maxima,minima)))
        
        """obtain the exrema values"""
        extrema = []
        for i in extrema_chain:
            if len(i[0]) > 0:
                for j in i[0]:
                    extrema.append(j)
        extrema = sorted(extrema)
    
        """calculate the score"""
        score = 20
        if len(extrema) > 1:
            count = 0
            sumfactor = 0
            for i in range(len(extrema)):
                index_start = extrema[i]
                index_end   = extrema[i+1]
                numerator   = abs(data[index_start] - data[index_end])
                denominator = data[index_start] + data[index_end]
                factor      = min(1, abs(data[index_start] - data[index_end]))
                sumfactor   += (float(numerator)/float(denominator))*factor
                count += 1
                if count == len(extrema)-1:
                    break
            score -= 2*sumfactor
        return score

""" this is the generalized function that takes a list of tasks and performs the scoring if needed
we place it in task object function dictionaries to prevent use of exec and to feed string to scores
making it in essence dynamic"""

def data_evaluation(data,tasks,attr= []):
    evaluation = {}
    for i in tasks:
        evaluation[i] = evaluate[i](data,attr = attr) 
    return evaluation



