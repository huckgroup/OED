
"""
@author: Bob van Sluijs
"""
"""
NOTE: THERE IS A CONFLICT BETWEEN AMICI AND TELLURIUM
THE LATTER TRANSLATES OUR TEXT BASED MODEL TO SMBML, THIS
SBML FILE IS THEN USED TO CREATE AMICI C++ OBJECT, AFTER CONVERSION
TO SBML MODEL.PY WILL RESTART THE SCRIPT AND LOAD THE DIRECTORIES IN
A DIFFERENT ORDER SO AMICI HAS ACCESS TO THE RIGHT PACKAGES
(If you have not read the READme installing amici can take some effort please check their website)
For any questions: contact bob.vansluijs@gmail.com
"""
"""the libaries that are needed"""
import math
import copy
import numpy  
import time as timing
import matplotlib.pylab as plt

"""import the model"""
import DEMO_model as model
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

