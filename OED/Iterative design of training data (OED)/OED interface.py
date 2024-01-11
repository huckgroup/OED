
"""
Created on Mon Jan 10 12:09:59 2022

@author: Bob van Sluijs
"""

"""the libaries that are needed"""
import math
import copy
import numpy  
import time as timing
import matplotlib.pylab as plt

"""import the model"""
import __NucleotideSalvagePathway__ as model
models,control = model.main() 

"""compile the imported models to sbml, stored on desktop in __SBML__ folder and compile to AMICI stored in __AMICI__
using the name of the model (as foldername)"""

"""please insure you have installed AMICI, tellurium and libSBML for a comprehensive overview see READ me"""
for i in range(len(models)):
    model = models[i]
    model.SBMLconversion()
    model.PytoCompile()

"""you can load multiple models simultenously, we only have one"""    
NSP_model = models[0]

"""import modules for experiment generation and identifiability analysis"""
from Measurements import GenerateExperiments
m = GenerateExperiments(NSP_model)
                        
"""build an in silico measurement object and use this to train the model"""
measurement = m.return_random_pulse_measurement(store = False,name = 'RATE_ID_TEST',time = (0,12*60),plength = 20)
measurement.show(show=True)

"""the scripts that are needed"""
from ParameterOptimizationModelSets import BoxplotOptimization

"""define the model that belongs to this measurement"""
measurement.model = NSP_model

"""fit the data to the model and store parmeter estimates (individually in folder)"""
directory = " [PATH], folder will be created at this path if not it will use this name in the current folder "  #folder where estimtates are stored
BoxplotOptimization({0:measurement},
                      optimization_number = 50,
                      generations = 150,
                      agents = 5,
                      storedata = directory,
                      startsamples = 100)

