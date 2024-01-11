# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:26:31 2018

@author: bobva
"""

import math
import numpy

d_func = {}
task = lambda f: d_func.setdefault(f.__name__, f)

@task
def log_uniform_dist(lower,upper,n):
    set_dist = numpy.logspace(numpy.log10(lower),numpy.log10(upper),n) #put this in task format and keep it optional
    return {i:round(set_dist[i],6) for i in range(n)}  
 
@task
def log_normal_dist(mean,sd,n):
    set_dist = list(sorted(numpy.random.lognormal(mean, sd, n)))
    return {i:round(set_dist[i],6) for i in range(n)}

@task
def uniform_dist(lower,upper,n):   
    set_dist = numpy.linspace(lower,upper,n)	
    return {i:round(set_dist[i],6) for i in range(n)}
            
            
            

