# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 12:09:59 2022

@author: Bob van Sluijs
"""


import os 
import sys

"""the path of the dependencies i.e. where you store the code, add it to current path"""
path = os.path.join(os.path.join(os.path.expanduser('~')), 'Onedrive\\Desktop\\OED\\OED\\')
path = sys.path.append(path)

"""the libaries that are needed"""
import math
import copy
import numpy  
import time as timing
import matplotlib.pylab as plt

"""the dependencies that are needed"""
from SolveSystem import ModelSolver
from Model import ModelVariables

"""Operations and parse"""
from Operations import *
from Parse import *
    

"""import the model"""
import ToyModel as model
models = model.main() 

"""create th SBML model and Amici models, needed because packacges conflict"""
for i in range(len(models)):
    model = models[i]
    model.SBMLconversion()
    model.PytoCompile()
    
"""import modules for experiment generation and identifiability analysis"""
from Measurements import GenerateExperiments
from IdentifiabilityAnalysis import IdentifiabilityAnalysis

"""name of folder where identifiability of experiment will be stored"""
name = model.name + 'First Test'

"""create a 'fake' measurement objects in the GenerateExperiment class 
(you can do batch, random pulse etc. just check measurement folder)
this returns a dict with a number of analyses"""
m = GenerateExperiments(model)
measurement_1 = m.return_random_pulse_measurement(store = True,name = name,time = (0,100))
measurement_2 = m.return_random_pulse_measurement(store = True,name = name,time = (0,100))

"""of you want to add an optimized pulse"""
# measurement_3 = m.optimized_pulse_measurement(store = True,name = name,time = (0,100))
#...
#...
#...

"""its a dict {0:Measurement object}, lets combine 2  experiments together and 
give them to the identifiabilitiy analysis (i.e. it combines them and does the idenfiability analysis as if its 1)"""
dataset = {}

"""build a dataset, Note we only have 1 model but 
the dataset can consists of measurements with multiple models"""
dataset[len(dataset)] = measurement_1[0]
dataset[len(dataset)] = measurement_2[0]


"""ANALYSE: Collinearity index, correlation and PCA"""
"""plots will be stored in a folder called identification on your desktop!"""
EXP = IdentifiabilityAnalysis(dataset)
EXP.Correlation(name = name)
EXP.PCA(name = name)
EXP.SVD(name = name) 


from ParameterOptimization import BoxplotOptimization
"""FIT: fit the measurements and retrieve parameter estimates"""
fit = BoxplotOptimization(dataset,
                     optimization_number = 5,
                     generations = 100,
                     agents = 10,
                     storedata = '[LOCATION]',
                     startsamples = 1000)

# fit.analyse()