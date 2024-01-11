# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:18:09 2019

@author: huckg
"""

""""turn yur script in an autorun module that way you can combine python 2 and 3 in 1 program
Note that your boundaries cannot be zero or fixed, you can set your fixed to zero if needed and not add them
to the "include" variable"""
import random
import numpy

from Parse import *
from Measurements import *
from OptimizationOperations import *
from Scoring import __score__

def simulate_experiment(model,time_dependent_parameters = {},formatted = False,conditions = {},time = (0,100),dt = 1,sf = 'ESS',forward = True):
    scores = []
    """reset the variables of the model"""
    variables = ModelVariables(model,conditions = conditions)  
    """"variables, update it with a pulse sequence"""
    if len(time_dependent_parameters) != 0:
        if formatted == False:
            variables.set_time_dependent_variables(time_dependent_parameters,model.pulsespace)            
            """solve the system using """
            solution = ModelSolver(model,variables = variables,simtime = time,dt = dt,forward_sensitivity = forward)
        else:
            solution = ModelSolver(model,variables = variables,manual_TDI = time_dependent_parameters,simtime = time,dt = dt,forward_sensitivity = forward)

    data = solution.__getData__()
    """assess the fitness"""
    metric = assessment[sf](simdata = data, model = model)
    simdata = data.simdata
    """get forward sensititivites"""
    forward_sensitivities = data.forward_sensitivities
    for fnc,value in metric.items():
        scores.append(__score__(score = value,observable = model.observables,criterion = fnc))     
    return scores,simdata,forward_sensitivities

def find_nearest_value(space,value):
    number = min(space.values(), key=lambda x:abs(x-value)) 
    inverted = {v:k for k,v in space.items()}
    return inverted[number]

def recombine_pulsepattern(inputvectors):
    '''Get 1 randomly selected parents for all parameters'''
    sequence = [random.choice(range(len(inputvectors))) for i in range(len(inputvectors[-1]))]
    p = sorted(inputvectors[-1].keys())
    '''fill up a new child with the right sequence''' 
    child = {}
    for i in range(len(sequence)):
        child[p[i]] = inputvectors[sequence[i]][p[i]]
    return child   

def gerenate_mutated_pattern(coordinates,
                             pindex  = 10,
                             mutation_size = 10,
                             mutation_number = 1):
    
        times = list(sorted(list(coordinates.values())[-1].keys()))
        """mutatory list"""
        mutationlist = []
        """loop through the mutations"""
        for i in range(mutation_number):
            """randomly chosen size"""
            length = random.choice(range(mutation_size))
            """parameter selected for mutation"""
            parameter = random.choice(list(coordinates.keys()))
            """time of the pulse"""
            pulse = random.choice(range(len(times)))
            """pulse go left or right"""
            sequence = []
            if random.random() <= 0.5:
                """go left along the number line [(),(),()...<- (chosen),(),(),()]"""
                start_index = pulse-length
                if start_index < 0:
                    start_index = 0
                for t in range(start_index,pulse,1):
                    sequence.append(times[t])

            else:
                """go right along the number line [(),(),(),(chosen) ->...,(),(),()]""" 
                end_index = pulse+length
                if end_index > len(times)-1:
                    end_index = len(times)
                for t in range(pulse,end_index,1):
                    sequence.append(times[t])

            """size of the index space"""
            size = random.choice(range(pindex))
            """append mutation"""    
            if len(sequence) != 0:
                mutationlist.append((parameter,sequence,size))

        """tracking mutations"""
        for parameter,sequence,index in mutationlist:
            for i in sequence:
                """update the coordinate vector"""
                coordinates[parameter][i] = index 
        return coordinates

def generate_random_input(pulsevector,largest_coordinate):
    """"rewrite the mutaiton options as boolean"""
    parameters = pulsevector.keys()
    """"parameter vector"""
    conditional_vector = []
    """mutation vector"""
    mutation_vector = numpy.linspace(0,1,4)
   
    """cooridnate of parameter"""
    coordinate_control  = {}
    for parameter,vector in pulsevector.items():
        coordinate_array = {}
            
        for time,coordinate in pulsevector[parameter].items():
            """Probability of mutation vector"""
            probabilities = sorted(mutation_vector)
            """include the random prob"""
            p = random.random()
            """smallest indices"""
            difference = [abs(probabilities[n]-p) for n in range(len(probabilities))]
            """selection of vector"""
            select = difference.index(min(difference))
            """select the mutation"""    
            if select == 0:
                """move down 1 or to indices"""
                if pulsevector[parameter][time] < 2:
                    coordinate_array[time] = pulsevector[parameter][time] + random.choice([1,2])                         
                elif pulsevector[parameter][time] >= 2:
                    coordinate_array[time] = pulsevector[parameter][time] - random.choice([1,2]) 
                                   
            elif select == 1:
                """move up 1 or to indices"""
                if pulsevector[parameter][time] > largest_coordinate - 3:
                    coordinate_array[time] = pulsevector[parameter][time] - random.choice([1,2])                         
                elif pulsevector[parameter][time] <= largest_coordinate - 3:
                    coordinate_array[time] = pulsevector[parameter][time] + random.choice([1,2]) 

            elif select == 2:
                """randomly select coordinate"""
                coordinate_array[time] = random.choice(range(largest_coordinate))
                    
            elif select == 3:
                coordinate_array[time] = pulsevector[parameter][time]

        """update overall conditions with parameters"""
        coordinate_control[parameter] = coordinate_array
    return coordinate_control

def smooth_start_pulse_conditions(pulsevector,largest_coordinate):
    """"rewrite the mutaiton options as boolean"""
    parameters = pulsevector.keys()
    """"parameter vector"""
    conditional_vector = []
    """mutation vector"""
    mutation_vector = numpy.array([0.05,0.1])
   
    """parameters in vector"""
    boolean_control  = {}
    """cooridnate of parameter"""
    coordinate_control  = {}
    for parameter,vector in pulsevector.items():
        coordinate_array = {}

        for time,coordinate in pulsevector[parameter].items():
            """Probability of mutation vector"""
            probabilities = sorted(mutation_vector)
            """include the random prob"""
            p = random.random()
            """smallest indices"""
            difference = [abs(probabilities[n]-p) for n in range(len(probabilities))]
            """selection of vector"""
            select = difference.index(min(difference))
            
            if select == 0:
                """randomly select coordinate"""
                coordinate_array[time] = random.choice(range(largest_coordinate))
                    
            elif select == 1:
                coordinate_array[time] = pulsevector[parameter][time]
                
        """update overall conditions with parameters"""
        coordinate_control[parameter] = coordinate_array
            
        """update overall vector"""
        conditional_vector.append((coordinate_control,pulsevector))   
    return conditional_vector

def sort_list(reference,output):
    indices = numpy.argsort(reference)
    output  = [output[out] for out in indices]
    reference = [reference[out] for out in indices]
    return reference,output

def time_dependent_mutationmeta(start,end,plenght):
    return random.choice(range(smallest,largest,1)), random.choice([0.51+(i*0.01) for i in range(len(0,50,1))])