# -*- coding: utf-8 -*-
"""
Created on Wed Feb 06 10:57:16 2019

@author: huckg
"""
from __future__ import division
import numpy
import matplotlib.pylab as plt
import copy
from scipy.signal import argrelextrema
from Operations import *
from DataEvaluation import *
from itertools import chain

class __score__:
    def __init__(self, number = '',score = '',observable = '',criterion = '',sensitivity = '',weight = 1,parameters = [],conditions = []):
        ''' the score is given'''
        self.score      = score
        ''' for an observalbe'''
        self.observable = observable
        ''' given a criterion function''' 
        self.criterion  = criterion
        '''within a series of measurements assigned a number
        in the __measurement__ object all with potentially
        and differnt scoring criteria for different observables'''
        self.number     =  number
        """weight given to the score biasing it to this scoring function
        could be usefull in light of the fitness landscape"""
        self.weight     = weight
        """"parameters are stored in this object to plot sensitivities afterwards"""
        self.parameters = parameters
        """"conditions are stored to plot sensitivity to conditions """
        self.conditions = conditions        
        """the sensitivtiy ODE that is scored"""
        self.sensitivity = sensitivity
        
        
assessment = {}
task = lambda f: assessment.setdefault(f.__name__, f)
###############################################################################
""""fitting function that the algorithm uses to optimize towards a measurement (data)"""
@task
def E_Fisher(simdata = [],measurement = [],model = [],criteria = False):
    assessment = {}
    tasklist = ['median parameter collinearity']
    if criteria:
        return (tasklist,[])
    """create the sensitivity vector"""
    sensitivity_vector = {i:[] for i in model.sensitivity_parameters}
    for o in model.observables:
        for p in model.sensitivity_parameters:
            sensitivity_vector[p].append(simdata.forward_sensitivities[p][o])

    """sensitivity vectors"""
    sensitivity_vector = {k:numpy.concatenate(v) for k,v in sensitivity_vector.items()}
    """the indices you take to calculate the eigenvalues with, more points = higher resolutions
    but is takes quite a bit longer to calculate the actual SVD values (larger polynomial)"""
    for parameter,values in sensitivity_vector.items():
        sensitivity_vector[parameter] = numpy.array([values[i] for i in range(0,len(values),5)])
    """"perform a single value decomposition"""
    collvalues,cumsense = single_value_decomposition(sensitivity_vector,model.observables)
    """extract collinearities"""
    collinearities = [i.collinearity for i in collvalues]
    """"we want this value to be as low as possible thus we maintain this"""
    assessment["E_Fisher"] = numpy.mean(collinearities)
    """return the assessment"""
    return assessment

@task
def D_Fisher(simdata = [],measurement = [],model = [],criteria = False):
    assessment = {}
    tasklist = ['median parameter collinearity']
    if criteria:
        return (tasklist,[])
    """create the sensitivity vector"""
    sensitivity_vector = {i:[] for i in model.sensitivity_parameters}
    for o in model.observables:
        for p in model.sensitivity_parameters:
            sensitivity_vector[p].append(simdata.forward_sensitivities[p][o])
    """sensitivity vectors"""
    sensitivity_vector = {k:numpy.concatenate(v) for k,v in sensitivity_vector.items()}
    """the indices you take to calculate the eigenvalues with, more points = higher resolutions
    but is takes quite a bit longer to calculate the actual SVD values (larger polynomial)"""
    for parameter,values in sensitivity_vector.items():
        sensitivity_vector[parameter] = numpy.array([values[i] for i in range(0,len(values),3)])
    """perform a single value decomposition"""
    array = numpy.array([sensitivity_vector[i] for i in sorted(sensitivity_vector.keys())])
    FIS   = numpy.matmul(array,array.T)
    sign, log = numpy.linalg.slogdet(FIS)
    if log < 0:    
        assessment["D_Fisher"] = abs(log)
    else:
        assessment["D_Fisher"] = 1./log 
    """return the assessment"""
    return assessment
