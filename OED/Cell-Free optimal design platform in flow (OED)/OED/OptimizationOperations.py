    # -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 12:40:18 2019

@author: vansl
"""
import numpy
import random
import math
import copy
import matplotlib.pylab as plt 
import statsmodels.api as sm

from Scoring import *
from Scoring import __score__
from SolveSystem import ModelSolver
from Model import ModelVariables

from Operations import *
from Distributions import *
from DataTransform import *


#for leas squares fit uncomment this one and comment the function below
# def calculate_deviation(simdata,measurement):
#     m,SD = 0.0,[]
#     """loop through the measurements"""
#     for n in measurement.keys():
#         for state in measurement[n].observables:
             
#               """update the counter"""
#               m += 1
#               """take difference"""
#               maximum   = float(numpy.amax(measurement[n].profile[state]))
#               deveation = sum((simdata[n][state] - measurement[n].profile[state])**2)
#               """the SD in the list"""
#               SD.append(copy.copy(deveation))
#     """average deveation"""
#     average = sum(SD)/float(m)     
#     return average

def calculate_deviation(simdata,measurement):
    m,SD = 0.0,[]
    """loop through the measurements"""
    for n in measurement.keys():
        for state in measurement[n].observables:
             
              smean = numpy.mean(simdata[n][state])
              mmean = numpy.mean(measurement[n].profile[state])
              diff = smean/mmean
              if diff < 1:
                  diff = 1./diff
              """update the counter"""
              m += 1
              """take difference"""
              maximum   = float(numpy.amax(measurement[n].profile[state]))
              N = len(simdata[n][state])
              deveation = sum((2.0/N * np.abs(numpy.fft.fft(simdata[n][state])[:N//2]) - (2.0/N * np.abs(numpy.fft.fft( measurement[n].profile[state]))[:N//2]))**2)
              """the SD in the list"""
              SD.append(copy.copy(deveation)*diff)

    """average deveation"""
    average = sum(SD)/float(m)    
    return average

def simulate_measurement(model,measurements,parameters,forward = False,sf = ""):
    scores = []
    """store the simulation data of the different measuremenst"""
    simdata,forward_sensitivities = {},{}
    """loop through the measurements and simulate the measurements"""
     
    for nmbr, measure in measurements.items(): 
        """reset the variables of the model"""
        variables = ModelVariables(model,conditions = measure.conditions,modification = parameters)  
        time_dependent_parameters = measure.time_dependent_parameters
        """solve the system using """
        solution = ModelSolver(model,variables = variables,measurement_time = measure.time,manual_TDI = time_dependent_parameters,forward_sensitivity = forward)
        """the data of the simulation"""
        data = solution.__getData__()
        """simdata nmbr"""
        simdata[nmbr] = data.simdata
        """forward sensitivities"""
        forward_sensitivities[nmbr] = data.forward_sensitivities
        
        """insim scoring to fit to data"""
        if sf == "leastsquaresfit":
            lsq = leastsquaresfit(simdata = simdata[nmbr],measurement = measure,model = model)
            for state,information in lsq.items():
                for criterion,value in information.items():
                    scores.append(__score__(score = value,observable = state,criterion = criterion))

        
    """score the simulation"""
    if sf == "standardeviation":
        average_SD = calculate_deviation(simdata,measurements)
        """update scores"""
        scores.append(__score__(score = average_SD,observable = model.observables,criterion = 'average'))
    return scores,simdata,forward_sensitivities

"""incorporate the simulation in the model"""
def simulate_measurements(measurements,parameters,forward = False,sf = ""):
    scores = []
    """store the simulation data of the different measuremenst"""
    simdata,forward_sensitivities = {},{}
    """loop through the measurements and simulate the measurements"""
     
    for nmbr, measure in measurements.items(): 
        model = measure.model
        """the model has explicit variables that it allows thus we need to prune those that do not exist"""
        pset = {}
        for p,value in parameters.items():
            if p in model.parameters:
                pset[p] = value
            else:
                pass

        """reset the variables of the model"""
        variables = ModelVariables(model,conditions = measure.conditions,modification = pset)  
        time_dependent_parameters = measure.time_dependent_parameters
        """solve the system using """
        solution = ModelSolver(model,variables = variables,measurement_time = measure.time,manual_TDI = time_dependent_parameters,forward_sensitivity = forward)
        """the data of the simulation"""
        data = solution.__getData__()
        """simdata nmbr"""
        simdata[nmbr] = data.simdata
        data.rdata
        
        """forward sensitivities"""
        forward_sensitivities[nmbr] = data.forward_sensitivities

        """insim scoring to fit to data"""
        if sf == "leastsquaresfit":
            lsq = leastsquaresfit(simdata = simdata[nmbr],measurement = measure,model = model)
            for state,information in lsq.items():
                for criterion,value in information.items():
                    scores.append(__score__(score = value,observable = state,criterion = criterion,model = model.name))

    """score the simulation"""
    if sf == "standardeviation":
        average_SD = calculate_deviation(simdata,measurements)
        """update scores"""
        scores.append(__score__(score = average_SD,observable = model.observables,criterion = 'average'))
    return scores,simdata,forward_sensitivities

         
""" these functions simulate measured data and score it according to the fit function
stored in the measuremnt"""
def sort_by_attr(seq, attr):
    import operator
    intermed = map(None, map(getattr, seq, (attr,)*len(seq)),
        range(len(seq)), seq)
    intermed.sort( )
    return map(operator.getitem, intermed, (-1,)*len(intermed))

def sort_by_attr_inplace(lst, attr):
    lst[:] = sort_by_attr(lst, attr)   

""" these are two multustart functions which canb be used to scour the fitness landscape prior to an optimization"""
def multistart(multistart,model,include,fitdata = [],agents = 2,start = True,sf = ''):
    scores = {i:[] for i in range(len(multistart))}
    for i in range(len(multistart)):
        scores[i],simdata,fwd = simulate_measurement(model,fitdata,multistart[i], sf = sf)
    return scores

def multimodalstart(multistart,measurements,start = True,sf = ''):
    scores = {i:[] for i in range(len(multistart))}
    for i in range(len(multistart)):
        scores[i],simdata,fwd = simulate_measurements(measurements,multistart[i], sf = sf)
    return scores


""" the following functions are to be used as operations within and evolutioanry algorithm"""
def recombine(parents):
    '''Get 1 randomly selected parents for all parameters'''
    sequence = [random.choice(range(len(parents))) for i in range(len(parents[-1]))]
    p = sorted(parents[-1].keys())
    '''fill up a new child with the right sequence''' 
    child = {}
    for i in range(len(sequence)):
        child[p[i]] = parents[sequence[i]][p[i]]
    return child  

def weighted_random_choice(choices):
    max = sum(choices.values())
    pick = random.uniform(0, max)
    current = 0
    for key, value in choices.items():
        current += value
        if current > pick:
            return key
        
def mutationrange(space):
    number = 2
    size   = random.choice(range(1,int(space/2),1))
    rnd    = random.choice([0,1])
    return number,size,rnd

def directed_mutation_range(space, percentage = (0.05,0.15),m_nmbr = 1):
    number = random.choice(range(1,m_nmbr,1))
    number = 1
    """set the size of the mutation"""
    size = int(space*random.choice(numpy.linspace(percentage[0],percentage[1],100)))
    return number,size
        
def select(scores,select = 2): 
    sort   = {i:[] for i in range(len(list(scores.values())[-1]))}
    weight = {i:list(scores.values())[-1][i].weight for i in range(len(list(scores.values())[-1]))} 
    for i in range(len(sort)):
        for j in range(len(scores)):
            sort[i].append(scores[j][i].score) 
            
    ranks = numpy.array(([0]*len(sort[0])))
    for i in range(len(sort)):    
        """"include previous generation generations fittest structures! to make sure you do not lose fitness"""
        scr = sort[i] #loop through scores
        output = [0] * len(scr) #create output file
        for j, x in enumerate(sorted(range(len(scr)), key=lambda y: scr[y])): #enumerate the different scores in the lists and rank them accoriding to their enumerated number
            output[x] = j #update output
        
        ranks += numpy.array(output) * weight[i] #final ranks    
    indices = numpy.argsort(ranks)
    selection = [indices[i] for i in range(select)]
    return selection  

def set_selection_treshold(measurements):
    treshold = random.choice(numpy.linspace(0.51,1,100))
    """the treshold cannot be lower than 50 percent!"""
    msr = 0
    for nmbr,measure in measurements.items():
        msr += len(measure.profile)
    number =int(treshold*msr)
    return number

def set_static_selection_treshold(criteria):
    treshold = random.choice(numpy.linspace(0.51,1,100))
    """the treshold cannot be lower than 50 percent!"""
    msr = len(criteria)
    """set criteria"""
    number =int(treshold*msr)
    """return the number"""
    return number
