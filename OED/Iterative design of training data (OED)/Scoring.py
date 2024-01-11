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
def remove_indices(arr, indices):
    import numpy as np
    # Convert the list of indices to a set for faster lookup
    index_set = set(indices)

    # Get the number of rows and columns in the original array
    num_rows, num_cols = arr.shape

    # Create a mask to exclude the specified rows and columns
    row_mask = np.array([i not in index_set for i in range(num_rows)])
    col_mask = np.array([i not in index_set for i in range(num_cols)])

    # Use the masks to create a new array with rows and columns removed
    new_arr = arr[row_mask][:, col_mask]

    return new_arr


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
    print(assessment)
    """return the assessment"""
    return assessment

@task
def D_Fisher(simdata = [],measurement = [],model = [],criteria = False):
    assessment = {}
    tasklist = ['median parameter collinearity']
    if criteria:
        return (tasklist,[])
    """create the sensitivity vector"""
    sensitivity_vector = {}
    for o in model.observables:
        for p in model.sensitivity_parameters:
            if 'k_cat' in p:
                if p not in sensitivity_vector:
                    sensitivity_vector[p] = []
                sensitivity_vector[p].append(simdata.forward_sensitivities[p][o])
            if 'K_M' in p:
                if p not in sensitivity_vector:
                    sensitivity_vector[p] = []
                sensitivity_vector[p].append(simdata.forward_sensitivities[p][o])
    """sensitivity vectors"""
    sensitivity_vector = {k:numpy.concatenate(v) for k,v in sensitivity_vector.items()}
    """the indices you take to calculate the eigenvalues with, more points = higher resolutions
    but is takes  quite a bit longer to calculate the actual SVD values (larger polynomial)"""
    for parameter,values in sensitivity_vector.items():
        sensitivity_vector[parameter] = numpy.array([values[i] for i in range(0,len(values),3)])
    """perform a single value decomposition"""
    array = numpy.array([sensitivity_vector[i] for i in sorted(sensitivity_vector.keys())])
    FIS   = numpy.matmul(array,array.T)
    # print(sensitivity_vector.keys())
    """data and indices"""
    data    = []
    indices = []
    
    """check the boolean"""
    j = 0
    for i in FIS:
        data.append(i)
        testbool = False
        for n in i:
            if n != 0:
                testbool = True
        if testbool == False:
            indices.append(j)
        j += 1
    FIS = remove_indices(FIS,indices)
    sign, log = numpy.linalg.slogdet(FIS)
    print(sign,log)
    if log < 0:    
        assessment["D_Fisher"] = abs(log)
    else:
        assessment["D_Fisher"] = 1./log 
    """return the assessment"""
    return assessment
        
    
@task
def ofit(simdata = [],measurement = [],model = [],criteria = False):
    tasklist = ['mean','amplitude','period','oscillations']
    if criteria:
        return (tasklist,['oscillations'])
    assessment = {state:'' for state in measurement.observables}    
    for state in measurement.observables:
        y = simdata[:,model.map[state]]
        attr = extremasearch(y)
        evaluation = data_evaluation(y,tasklist,attr = attr)
        assessment[state] = {i:abs(evaluation[i] - measurement.evaluation[state][i]) for i in tasklist}
    return assessment

@task
def fourierfit(simdata = [],measurement = [],model = [],criteria = False):
    tasklist = ['fouriertransform']
    if criteria:
        return (tasklist,[])
    assessment = {state:'' for state in measurement.observables}    
    for state in measurement.observables:
        y = simdata[:,model.map[state]]
        attr = measurement.dt
        evaluation = data_evaluation(y,tasklist,attr = attr)
        try:
            assessment[state] = {i:numpy.sum((evaluation[i] - measurement.evaluation[state][i])**2) for i in tasklist}
        except ValueError: #operands could not be broadcast together, solver error
            assessment[state] = {i:float('inf') for i in tasklist}
    return assessment

@task        
def leastsquaresfit(simdata = [],measurement = [],model = [],criteria = False):
    tasklist = ['leastsquares']
    if criteria:
        return (tasklist,[])    
    assessment = {state:'' for state in measurement.observables if state in model.states}  
    for state in assessment.keys():
        try:
            assessment[state] = data_evaluation(simdata[state],tasklist,attr = measurement.profile[state])
        except ValueError: #operands could not be broadcast together, solver error
            raise ValueError("The LSQ score is adding and or substracting wrong NaN values or arrays or mismatched in shape")
            assessment[state] = {i:float('inf') for i in tasklist}
    return assessment

@task
def standardeviation(simdata = [],measurement = [],model = [],criteria = False):   
    tasklist = ['standard_deviation']
    if criteria:
        return (tasklist,[])    
    assessment = {state:'' for state in measurement.observables if state in model.states}     
    average_SD,n = 0,0
    for state in measurement.observables:
         """update the counter"""
         n += 1
         """take difference"""
         factor = ( numpy.array(simdata[state]) / numpy.array(measurement.profile[state]))   
         """take absolute values"""
         deviation = []
         for i in factor:
             if i <= 1:
                 deviation.append(i**-1)
             else:
                 deviation.append(i)
         SD = numpy.mean([abs(i - 1) for i in deviation])
         """update average SD"""
         average_SD += SD
    average_SD = average_SD/float(n)
    for observable,v in assessment.items():
        assessment[observable] = {"standard_deviation":average_SD}
    return assessment

@task
def LSQmeancorrected(simdata = [],measurement = [],model = [],criteria = False):   
    tasklist = ['LSQ mean corrected']
    if criteria:
        return (tasklist,[])    
    assessment = {state:'' for state in measurement.observables if state in model.states}     
    average_SD,n = 0,0
    for state in measurement.observables:
         """update the counter"""
         n += 1
         """take difference"""
         factor = sum(( numpy.array(simdata[state]) - numpy.array(measurement.profile[state]))**2)/ numpy.average(measurement.profile[state])
         """update average SD"""
         average_SD += SD
    average_SD = average_SD/float(n)
    for observable,v in assessment.items():
        assessment[observable] = {"LSQ mean corrected":average_SD}
    return assessment


@task
def sugarfit(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    assessment = {state:None for state,index in observables}
    tasklist = ['leastsquares','leastsquaresdifferential']
    if criteria:
        return (tasklist,[])  
    for state,indices in observables:
        y = numpy.zeros((len(simdata[-1])))
        for i in indices:
           y += simdata[model.map[i]]
        assessment[state] = numpy.sum((y - measurement[state])**2)
    return assessment

###############################################################################
""""evaluation function that the algorithm uses to assess general behavior ( nodata)"""
@task
def ESS(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    assessment = {state:None for state,index in observables}
    tasklist = ['mean','steadystates']
    if criteria:
        return (tasklist,[])
    for state,index in observables:
        y = simdata[:,index]
        evaluation = data_evaluation(y,tasklist)
        for k,v in evaluation.items():
            if numpy.isnan(v) or v < 0:
                evaluation[k] = 0
        assessment[state] = copy.copy(evaluation)
    return assessment

@task
def oscillations(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    assessment = {state:None for state in observables}
    print('bababaababababaq')
    tasklist = ['mean','amplitude','period','oscillations']
    if criteria:
        return (tasklist,['oscillations'])
    for state in observables:
        y = simdata[state]
        base = {i:False for i in tasklist}
        minima,maxima,extrema = extremasearch(y)
        values  = [y[i] for i in extrema]  
        evaluation = data_evaluation(y,tasklist,attr = (minima,maxima,extrema))
        for k,v in evaluation.items():
            if numpy.isnan(v) or v < 0:
                evaluation[k] = 0
        assessment[state] = evaluation
    return assessment

@task
def pHoscillations(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    import math
#    for k,v in simdata.items():
#        plt.plot(v,label = k)
#    plt.legend(fancybox = True)
#    plt.show()
#    for o in observables:
#        plt.plot(simdata[o],label = o)
#    plt.yscale("log")
#    plt.legend(fancybox = True)
#    plt.show()
#    
#    plt.plot(- numpy.log10(simdata["OH"]) - numpy.log10(simdata["H"]))
#    plt.show()
    assessment = {state:None for state in observables}
    tasklist = ['mean','amplitude','period','oscillations']
    if criteria:
        return (tasklist,['oscillations'])
    for state in observables:
        y = simdata[state]
        base = {i:False for i in tasklist}
        minima,maxima,extrema = extremasearch(y)
        values  = [y[i] for i in extrema]  
        evaluation = data_evaluation(y,tasklist,attr = (minima,maxima,extrema))
        for k,v in evaluation.items():
            if numpy.isnan(v) or v < 0:
                evaluation[k] = 0
        assessment[state] = evaluation
    return assessment


@task
def oeval(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    assessment = {state:None for state,index in observables}
    tasklist = ['mean','amplitude','period','oscillations']
    if criteria:
        return (tasklist,['oscillations'])
    for state,index in observables:
        y = simdata[:,index]
        base = {i:False for i in tasklist}
        minima,maxima,extrema = extremasearch(y)
        values  = [y[i] for i in extrema]  
        evaluation = data_evaluation(y,tasklist,attr = (minima,maxima,extrema))
        for k,v in evaluation.items():
            if numpy.isnan(v) or v < 0:
                evaluation[k] = 0
        assessment[state] = evaluation
    return assessment

@task
def biseval(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    assessment = {state:None for state,index in observables}
    tasklist = ["firststate","secondstate","stateratio",'bistable']
    if criteria:
        return (tasklist,["bistable"])
    for state,index in observables:
        y = simdata[:,index]
        base = {i:False for i in tasklist} 
        evaluation = data_evaluation(y,tasklist)
        for k,v in evaluation.items():
            if numpy.isnan(v) or v < 0:
                evaluation[k] = 0
            assessment[state] = evaluation            
    return assessment

@task
def retroaldolflux(simdata = [],measurement = [],model = [],criteria = False,observables = []):
    assessment = {state:None for state,index in observables}
    tasklist = ['aldolflux','retroaldolflux',"sugarbreakdown","sugarcycle"]
    if criteria:
        return (tasklist,[]) 
    y = numpy.zeros(len(simdata[-1]))
    for state,indices in observables:
        for i in indices:
            pass

###############################################################################
""""scoring function that the algorithm uses as optmization criteria (no data)"""
@task
def score_extrema(simdata = [],observables = [],model = [],criteria = False): 
    tasklist = ["extrema_score"]
    if criteria:
        return (tasklist,[])
    assessment = {state:'' for state in model.observables}     
    for state in model.observables:
        y = simdata[state]
        evaluation = data_evaluation(y,tasklist)
        for k,v in evaluation.items():
            if numpy.isnan(v) or v < 0:
                evaluation[k] = float('inf')
            assessment[state] = copy.copy(evaluation)  
    return assessment


@task
def score_ft(simdata = [],measurement = [],model = [],criteria = False): 
    tasklist = ["extrema_score"]
    if criteria:
        return (tasklist,[])
    assessment = {state:'' for state in model.observables}   

    
    for i in m_target:
        j = cp[:,i]
        j = j[25*60:]
        k = j - numpy.mean(j)
        trf = numpy.abs(numpy.fft.fft(k))**2
        ts = 1./3600.
        freqs = numpy.fft.fftfreq(k.size,ts)
        arg   = numpy.argsort(freqs)
        x = freqs[arg]
        cp = trf[arg]
    #this scoring function is tailored to the multiobjective optimization of systems (large amplitude vs high frequencies)
    mdl = int(len(x)/2)
    srf = numpy.trapz(cp[mdl:],x[mdl:])

    f_max = numpy.max(cp)
    idx_max    = (list(cp)).index(f_max)
    ls  = cp[0:idx_max]
    rs =  cp[idx_max:-1]   
    try:
        fmax = (abs(1./float(f_max)))
    except ZeroDivisionError:
        fmax = abs(1*10**10)     
    try:
        smax = (abs(float(srf)/float(f_max)))
    except ZeroDivisionError:
        smax = abs(1*10**10)  
    score = (fmax,smax)    
    return score
 
@task
def score_pw(simdata = [],measurement = [],model = [],criteria = False): 
    tasklist = ["extrema_score"]
    if criteria:
        return (tasklist,[])
    assessment = {state:'' for state in model.observables}   
    
    for i in m_target:
        j = cp[:,i]
        j = j[25*60:]
        k = j - numpy.mean(j)  
        trf = numpy.abs(numpy.fft.fft(k))**2
        ts = 1./3600.
        freqs = numpy.fft.fftfreq(k.size,ts)
        arg   = numpy.argsort(freqs)
        x = freqs[arg]
        cp = trf[arg]
    #this scoring function is tailored to the multiobjective optimization of systems (large amplitude vs high frequencies)
    mdl = int(len(x)/2)
    srf = numpy.trapz(cp[mdl:],x[mdl:])

    f_max = numpy.max(cp)
    if math.isnan(f_max):
        return (10*10,10*10,10*10)
    idx_max    = (list(cp)).index(f_max)
    ls  = cp[0:idx_max]
    rs =  cp[idx_max:-1]   
    
    lb_idx = numpy.searchsorted(ls,float(f_max)/3.)
    rb_idx = numpy.searchsorted(rs,float(f_max)/3.)  
    
    xl = x[lb_idx]
    xr=  x[rb_idx] 
    width = abs(xl - xr)
    if width == 0.0:
        width = 10
    ff = x[idx_max]
    try:
        fmax = (abs(1./float(f_max)))
    except ZeroDivisionError:
        fmax = abs(1*10**10)   
    try:
        smax = (abs(float(srf)/float(f_max)))
    except ZeroDivisionError:
        smax = abs(1*10**10)
    try:
        frq= abs(1./float(ff))
    except ZeroDivisionError:
        frq = abs(1*10**10)
    score = (fmax,smax,frq)    
    return score
         
#@task
#def parameter_identification(simdata = [],measurement = [],model = [],time = [],criteria = False):
#    """task list of the bayes factor comparing least squares and steady states"""
#    tasklist = ["steady_states"]
#    if criteria:
#        return (tasklist,[]) 
#
#    assessment = {state:'' for state in model.observables}
#    assessment.update({state:'' for state in model.lumped.keys()})
#    
#    """this is for observables, specifically individual components"""
#    for state in model.observables:
#        y = simdata[:,model.map[state]]
#        evaluation = data_evaluation(y,tasklist,attr = time)
#        assessment[state] = evaluation
#
#    """this is for observables lumped, specifically lumped sugars"""
#    lo = {i:numpy.zeros((len(simdata[:,-1]))) for i in model.lumped.keys()}
#    for l_obs,indices in model.lumped.items():
#        for i in indices:
#            """y vector with data corresponding to a """
#            y = simdata[:,i]
#            lo[l_obs] += y
#            
#    for state,y in lo.items():
#        evaluation = data_evaluation(y,tasklist,attr = time)
#        assessment[state] = evaluation
#    return assessment