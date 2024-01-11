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


#
#def calculate_deviation(simdata,measurement):
#    m,SD = 0.0,[]
#    """loop through the measurements"""
#    for n in measurement.keys():
#        for state in measurement[n].observables:
#              smean = numpy.mean(simdata[n][state])
#              mmean = numpy.mean(measurement[n].profile[state])
#              diff = smean/mmean
#              if diff < 1:
#                  diff = 1./diff
#              """update the counter"""
#              m += 1
#              """take difference"""
#              maximum   = float(numpy.amax(measurement[n].profile[state]))
#              deveation = sum((simdata[n][state] - measurement[n].profile[state])**2)#/(numpy.average(measurement[n].profile[state])* len(measurement[n].profile[state]))
#              """the SD in the list"""
#              SD.append(copy.copy(deveation))
#              SD.append(copy.copy(deveation)*diff)
#              # absolute = [abs(i) for i in deveation]
#              # summed = numpy.sum(absolute)
#              # SD.append(copy.copy(deveation))
#             
#    """average deveation"""
#    average = sum(SD)/float(m)  
#    if average < 0:
#        average = 10**15
#    return average


# def calculate_deviation(simdata,measurement):
#     m,SD = 0.0,[]
#     """loop through the measurements"""
#     for n in measurement.keys():
#         add = []
#         for state in measurement[n].observables:
             
#               """update the counter"""
#               m += 1
#               """take difference"""
#               for i in range(len(simdata[n][state])):
#                   if simdata[n][state][i] < 0.5:
#                       add.append(simdata[n][state][i])
#                       simdata[n][state][i] = 0.1
#                   else:
#                       add.append(simdata[n][state][i])
                
              
#               add =  numpy.array(add)[200:1200]
#               add = [i for i in add if i < 0.5]
             
#               # print(numpy.amin(measurement[n].profile[state][200:1200]),"MIIIIIIIIIIIIIN")
#               # print(numpy.amax(measurement[n].profile[state][200:1200]),"MAAAAAAAAAAAAAX")
#               d1 = [simdata[n][state][i]/measurement[n].profile[state][i] for i in range(len(simdata[n][state]))]
#               d2 = [measurement[n].profile[state][i]/simdata[n][state][i] for i in range(len(simdata[n][state]))]                                                                        
#               deveation = numpy.sum((numpy.array(d1) + numpy.array(d2))/(2*len(simdata[n][state])))                                                             
#               print(deveation)
#               # deveation = sum((simdata[n][state] - measurement[n].profile[state])**2)#/(numpy.average(measurement[n].profile[state])* len(measurement[n].profile[state]))
#               """the SD in the list"""
#               if len(add) > 10:
#                   factor = 10
#               else:
#                   factor = 1
#               SD.append(copy.copy(deveation*factor))
#               # absolute = [abs(i) for i in deveation]
#               # summed = numpy.sum(absolute)

#     """average deveation"""
#     average = sum(SD)/float(m)   
#     print(average)
#     return average


#def calculate_deviation(simdata,measurement):
#    m,SD = 0.0,[]
#    
#    """loop through the measurements"""
#    
#    for n in measurement.keys():
#        for state in measurement[n].observables:
#             
#              smean = numpy.mean(simdata[n][state])
#              mmean = numpy.mean(measurement[n].profile[state])
#              diff = smean/mmean
#              if diff < 1:
#                  diff = 1./diff
#                  
#                  
#           
#              """update the counter"""
#              m += 1
#              """take difference"""
#              maximum   = float(numpy.amax(measurement[n].profile[state]))
#              N = len(simdata[n][state])
#              deveation = sum((2.0/N * np.abs(numpy.fft.fft(simdata[n][state])[:N//2]) - (2.0/N * np.abs(numpy.fft.fft( measurement[n].profile[state]))[:N//2]))**2)
#                # deveation = sum((simdata[n][state] - measurement[n].profile[state])**2)#/(numpy.average(measurement[n].profile[state])* len(measurement[n].profile[state]))
#              """the SD in the list"""
#              SD.append(copy.copy(deveation)*diff)
#              # absolute = [abs(i) for i in deveation]
#              # summed = numpy.sum(absolute)
#              # SD.append(copy.copy(deveation))
#             
#    """average deveation"""
#    average = sum(SD)/float(m)  
#    if average < 0:
#        average = 10**15
#    return average



#get the inter and intra set 

# def calculate_deviation(simdata,measurement):
#     m,SD = 0.0,[]    
#     for n in measurement.keys():
#         for state in measurement[n].observables:
#               """update the counter"""
#               m += 1
#               """take difference"""
#               d,df = [],[]
#               for i in range(len(simdata[n][state])):
#                   ds = (simdata[n][state][i]+1) / (measurement[n].profile[state][i]+1)
#                   if ds < 1:
#                       ds = 1./ds
#                   ds -= 1
#                   d.append(copy.copy(ds/len(measurement[n].profile[state])))
                  
#               for i in range(len(simdata[n][state])-1):
#                   diff_sim = numpy.gradient(simdata[n][state]) + 1.0
#                   diff_m   = numpy.gradient(measurement[n].profile[state]) + 1.0
#                   ds = abs(diff_sim[i] - diff_m[i])
#                   df.append(copy.copy(ds/len(measurement[n].profile[state])-1))
                  
#               """the SD in the list"""
#               df = numpy.array(df)/max(df)
#               SD.append(copy.copy(sum(d) + copy.copy(numpy.sum(df))))
              
#     average = sum(SD)/float(m)    
#     return average

# def calculate_deviation(simdata,measurement):
#         m,SD = 0.0,[]
#         """Loop through the measurements"""
#         for n in measurement.keys():
#             for state in measurement[n].observables:
#                   smean = numpy.mean(simdata[n][state])
#                   mmean = numpy.mean(measurement[n].profile[state])
#                   diff = smean/mmean
#                   if diff < 1:
#                       diff = 1./diff

#                   """Update the counter"""
#                   m += 1
                  
#                   """Take difference"""
#                   maximum   = float(numpy.amax(measurement[n].profile[state]))
#                   N = len(simdata[n][state])
#                   deveation = sum((simdata[n][state]/float(maximum) - measurement[n].profile[state]/float(maximum))**2)
                  
#                   """The SD in the list"""
#                   SD.append(copy.copy(deveation)*diff)

#         """average deveation"""
#         average = sum(SD)/float(m) 
#         # print(average)
#         if average < 0:
#             average = 10**15
#         return average
   
def calculate_error(a, b):
    # Find the smaller of the two numbers
    smaller = min(a, b)
    
    # Calculate the percentage difference
    if smaller == 0:
        return 0

    diff = abs(abs(a) - abs(b))
    percentage_diff = abs((diff / smaller) * 100)
    return percentage_diff

# def calculate_deviation(simdata,measurement):
#     """absolute difference"""
#     average = 0
    
#     """Loop through the measurements"""
#     n = 0
#     for number in measurement.keys():
#         for state in measurement[number].observables:
#             simulated_state = simdata[number][state]
#             """Error per state"""
#             true_error = []
#             """loop through the dataset"""
#             for time,data in measurement[number].profile[state]:
#                 try:
#                     error = calculate_error(simulated_state[time],data)
#                     true_error.append(error)
#                 except IndexError:
#                     pass
#             n+= 1

#             average += numpy.average(true_error)
#     return float(average/float(n))

def calculate_deviation(simdata,measurement):
    import math
    """absolute difference"""
    average = 0
    
    """Loop through the measurements"""
    n = 0
    for number in measurement.keys():
        for state in measurement[number].observables:
            simulated_state = simdata[number][state]
            """Error per state"""
            true_error = []
            """loop through the dataset"""
            for time,data in measurement[number].profile[state]:
                try:
                    error = numpy.sqrt((math.log2(simulated_state[time]) - math.log2(data))**2)
                    true_error.append(error)
                    # print(error,'boo')
                except:
                    pass
            n+= 1
            average += numpy.average(true_error)
    # print(float(average/float(n)))
    return float(average/float(n))
     

def simulate_measurement(model,measurements,parameters,forward = False,sf = ""):
    """This function takes the model and measurements, solved the model and 
    checks the fit of the data to the measurements"""
    scores = []
    
    """Store the simulation data of the different measuremenst"""
    simdata               = {}
    forward_sensitivities = {}
    
    """Loop through the measurements and simulate the measurements"""
    for nmbr, measure in measurements.items():
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
        
    """this is the fit error"""
    average_FE = calculate_deviation(simdata,measurements)
    """update score"""
    scores.append(__score__(score = average_FE,observable = model.observables,criterion = 'average'))
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
        """the time dependent inputs"""
        time_dependent_parameters = measure.time_dependent_parameters
        """solve the system using """
        solution = ModelSolver(model,variables = variables,measurement_time = measure.time,manual_TDI = time_dependent_parameters,forward_sensitivity = forward)
        """the data of the simulation"""
        data = solution.__getData__()
        """simdata nmbr"""
        simdata[nmbr] = data.simdata
        """forward sensitivities"""
        forward_sensitivities[nmbr] = data.forward_sensitivities
        
    """this is the fit error"""
    average_FE = calculate_deviation(simdata,measurements)
    """update scores"""
    scores.append(__score__(score = average_FE,observable = model.observables,criterion = 'average'))
    return scores,simdata,forward_sensitivities


"""incorporate the simulation in the model"""
def simulate_optimal_parametersets(measurements,parameters,forward = False,sf = ""):
    scores = []
    """store the simulation data of the different measuremenst"""
    simdata = {}
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
        """the time dependent inputs"""
        time_dependent_parameters = measure.time_dependent_parameters
        """solve the system using """
        solution = ModelSolver(model,variables = variables,measurement_time = measure.time,manual_TDI = time_dependent_parameters,forward_sensitivity = forward)
        """the data of the simulation"""
        data = solution.__getData__()
        """simdata nmbr"""
        simdata[nmbr] = data.simdata
    return simdata

         
""" these functions simulate measured data and score it according to the fit function
stored in the measuremnt"""
def sort_by_attr(seq, attr):
    import operator
    intermed = map(None, map(getattr, seq, (attr,)*len(seq)),
        range(len(seq)), seq)
    intermed.sort()
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
    import time as trackseconds
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


#def multistart(multistart,model,include,agents = 2,sf = "ESS",start = True,time = 100):
#    scores = {i:[] for i in range(len(multistart))}
#    for i in range(len(multistart)):
#        control = {k:v for k,v in multistart[i].items() if k in model.control}
#        """"obtain the control parameters or the experimental conditions"""
#        scores[i],simdata,fwd = simulate_system(model,multistart[i],conditions,sf = sf,time = time,dt = 1)
#    return scores
#
#
#
#    
#def simulate(model,parameters = {},time_dependent_parameters = {},time = (0,100),dt = 1,sf = 'ESS',forward = False,meta = False):
#    scores = []
#    """reset the variables of the model"""
#    variables = ModelVariables(model,modification = parameters)  
#    """"variables, update it with a pulse sequence"""
#    if len(time_dependent_parameters) != 0:
#        variables.set_time_dependent_variables(time_dependent_parameters,model.pulsespace)
#    """solve the system using """
#    solution = ModelSolver(model,variables = variables,simtime = time,dt = 1,forward_sensitivity = forward)
#    data = solution.__getData__()
#    """assess the fitness"""
#    metric = assessment[sf](simdata = data, model = model)
#    simdata = data.simdata
#    """get forward sensititivites"""
#    forward_sensitivities = data.forward_sensitivities
#    if not meta:
#        for state,scr in metric.items():
#            for fnc,value in scr.items():
#                scores.append(__score__(score = value,observable = state,criterion = fnc,parameters = parameters))
#    return

#def simulate_system(model,parameters,sequence,time = (0,100),dt = 1,sf = 'score'):
#    scores = []
#    """reset the variables of the model"""
#    variables = ModelVariables(model,modification = parameters)  
#    """"variables, update it with a pulse sequence"""
#    variables.generate_pulse(sequence)
#    """solve the system using """
#    solution = ModelSolver(model,variables = variables,time = time,dt = 1,forward_sensitivity = True)
#    data = solution.__getData__()
#    """assess the fitness"""
#    metric = assessment[sf](simdata = data, model = model)
#    simdata = data.simdata
#    """get forward sensititivites"""
#    forward_sensitivities = data.forward_sensitivities
#    for state,scr in metric.items():
#        for fnc,value in scr.items():
#            scores.append(__score__(score = value,observable = state,criterion = fnc,parameters = parameters))
#    return scores,simdata,forward_sensitivities
#  
#def simulate_likelihood(model,conditions,parameters,time,dt,forward = False):
#    """store the scores for the simulation"""
#    scores = []
#    """store the simulation data of the different measuremenst"""
#    simdata = {}
#    
#    for nmbr,condition in enumerate(conditions):
#        """reset the variables of the model"""
#        variables = ModelVariables(model,conditions = condition,modification = parameters)  
#        '''Time vectors which need to be similar in order to score the fit''' 
#        solution = ModelSolver(model,variables = variables,simtime = time,dt = dt,forward_sensitivity = forward)
#        data     = solution.__getData__()
#        """assess the fitness"""
#        metric = assessment[fitfunc](simdata = data.rdata, model = model)
#        """simulate the state and number"""
#        simdata[nmbr] = data.rdata
#        for state,scr in metric.items():
#            for fnc,value in scr.items():
#                scores.append(__score__(number = nmbr,score = value,observable = state,criterion = fnc,parameters = parameters,conditions = cnd,simdata = simdata))
#    return scores,simdata,forward
#     


# def calculate_deviation(simdata,measurement):
#     m,SD = 0.0,[]
#     """loop through the measurements"""
#     for n in measurement.keys():
#         for state in measurement[n].observables:
#              """update the counter"""
#              m += 1
#              """take difference"""
#              maximum   = float(numpy.amax(measurement[n].profile[state]))
#              deveation = (simdata[n][state]/maximum - measurement[n].profile[state]/maximum)
#              """take absolute values"""
#              absolute = [abs(i) for i in deveation]
#              summed = numpy.sum(absolute)/float(len((simdata[n][state])))
#              SD.append(copy.copy(summed))
#     """average deveation"""
#     average = sum(SD)/float(m)    
#     return average