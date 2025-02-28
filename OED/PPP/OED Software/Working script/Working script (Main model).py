# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 13:20:20 2023

@author: huckg
"""

import copy
import math

"""the libaries that are needed"""
import math
import copy
import numpy  
import time as timing
import matplotlib.pylab as plt
  
"""the scripts that are needed"""
from SolveSystem import ModelSolver
from Model import ModelVariables
from ParameterOptimizationModelSets import BoxplotOptimization
from Measurements import MeasurementObject

"""the modelbuilder scripts"""
from ModelBuilderOperations import *
from ParseMeasurementData import   *
from Angewandte_Figures import *


"""
This is a working script, it takes the two OED experiments, both the inputs given and the raw data
and converts that into a measurementObject, this is subsequently used to train a model, the model
is defined here by the reaction equations and constructed by the model builder i.e. it literally
takes the list with reaction equations and builds ODEs as shown in the demo model, and compiles with amici!
"""



def translate_TDI(TDI,marker = 'kflow(channel_{})',step = 12):
    for time,channel in TDI.items():
        alteration = {}
        for c,flow in channel.items():
            alteration[marker.format(str(c))] = flow/60.
        TDI[time] = alteration
    start = {list(sorted(TDI.keys()))[0]: TDI[list(sorted(TDI.keys()))[0]]}
    
    new_TDI = {}
    for time,vector in TDI.items():
        a,b = time
        new_TDI[(a,b)] = vector
    start.update(new_TDI)
    return start

# ################################################### --- OED ONE --- ################################################################################
"""Obtain the experimental data and time dependent inputs"""
path = "[Insert path of OED mearument folder]"
observables,time_dependent_inputs,known_parameters_M1 = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load,stock = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI,step = 12)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')

"""Prune the measurements"""
prune = (0,3,'ATP')
measurement.prune(prune)
prune = (0,3,'ADP')
measurement.prune(prune)
prune = (0,3,'AMP')
measurement.prune(prune)
prune = (0,3,'NADP')
measurement.prune(prune)

"""delete NADPH"""
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
del measurement.profile['NADPH']
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
total = sum(list(TDI[list(sorted(list(TDI.keys())))[0]].values()))
inverted_syringe_load = {i:'kflow(channel_{})'.format(k) for k,v in syringe_load.items() for i in v}

"""conditions updated from inflow stock"""
measurement.conditions = {state:value * TDI[list(sorted(list(TDI.keys())))[0]][inverted_syringe_load[state]]/total for state,value in stock.items()}
for state,data in measurement.profile.items():
        measurement.conditions[state] = data[1][1]

M = {}
M[(len(M))] = copy.copy(measurement)
################################################### --- OED TWO --- ##################################################################################
"""Obtain the experimental data and time dependent inputs"""
path = "[Insert path of OED mearument folder]"
observables,time_dependent_inputs,known_parameters_M2 = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load,stock = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI, step = 15)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')

"""Prune the first 3 minutes of measurement"""
prune = (0,3,'ATP')
measurement.prune(prune)
prune = (0,3,'ADP')
measurement.prune(prune)
prune = (0,3,'AMP')
measurement.prune(prune)
prune = (0,3,'NADP')
measurement.prune(prune)

"""delete NADPH"""
measurement.observables = [i for i in measurement.observables if i not in ['NADPH']]
try:
    del measurement.profile['NADPH']
except:
    pass

"""show the measurement"""
M[(len(M))] = copy.copy(measurement)
plot_raw_data(M[0],exp ='exp 1')
plot_raw_data(M[1],exp ='exp 2')


###############################################################################
"""Define the model of the system"""
#multiple ways to define a network in this instance we can make a list of lists
#where the entries are the defining terminology
###############################################################################
from random import shuffle            
modelnames = [
                'PPP_Literature_GH_rev',
                'PPP_Literature_inh_GH_rev',
              ]

for iterations in range(10):
    for Q in ['M1','M2']:
        for m in reversed(modelnames):   
                if m =='PPP_Literature_GH_rev': 
                    modelname = 'PPP_Literature_GH_rev'            
                    network  =    [
                                  ['HK'   , True, 'GH', ['Glucose','ATP'],['ADP','G6P'],[],[]],
                                  ['G6PDH', True, 'GH', ['NADP','G6P'],['6PGL','NADPH'],[],[]],
                                  [False  , True,False, ['6PGL'],['6PGA'],[],[]],
                                  ['6PGDH', True, 'GH', ['NADP','6PGA'],['RU5P','NADPH'],[],[]],  #NADPH, ATP,
                                  ['PRI'  , True, 'MM', ['RU5P'],['R5P'],[],[]],
                                  ['PRPPS', True, 'GH', ['ATP','R5P'],['PRPP','AMP'],[],[]]
                                  ]
                    
                if m =='PPP_Literature_inh_GH_rev': 
                    modelname = 'PPP_Literature_inh_GH_rev'            
                    network  =    [
                                  ['HK'   , True, 'GH', ['Glucose','ATP'],['ADP','G6P'],[],[]],
                                  ['G6PDH', True, 'GH', ['NADP','G6P'],['6PGL','NADPH'],[],[]],
                                  [False  , True,False, ['6PGL'],['6PGA'],[],[]],
                                  ['6PGDH', True, 'GH', ['NADP','6PGA'],['RU5P','NADPH'],[],[]],  #NADPH, ATP,
                                  ['PRI'  , True, 'MM', ['RU5P'],['R5P'],[],[]],
                                  ['PRPPS', True, 'GH', ['ATP','R5P'],['PRPP','AMP'],[],['ADP']]
                                  ]
                    
                ###############################################################################
                ###############################################################################
                """Get the reactions and categorize them"""
                reactions = categorize_reactions(network)
                """Build the equations for the reaction and reactor kinetics"""
                reactor_kinetics,reaction_kinetics = construct_model_fluxterms(reactions,reactormix = syringe_load)    
                """Construct the model"""
                basemodel = construct_kinetic_combinations(reactor_kinetics,reaction_kinetics,name = modelname)
                """directory where we store the results"""
                directory = "C:\\Users\\huckg\\OneDrive\\Desktop\\PPP Angewandte Adjusted\\{0} {1}\\".format(modelname,Q) 
          
                """add model to measurement"""
                if Q == 'M1':
                    measurement = M[0]
                if Q == 'M2':
                    measurement = M[1]
                
                if Q == 'M1':
                    model = basemodel.return_model(fixed = known_parameters_M1)
                if Q == 'M2':
                    model = basemodel.return_model(fixed = known_parameters_M2)
        
                """Prune the measurements"""
                measurement.model = model #assign model to measurement 
        
                BoxplotOptimization({0:measurement},
                                  include             = model.include,
                                  optimization_number = 10,
                                  generations         = 400,
                                  agents              = 25,
                                  storedata           = directory,
                                  startsamples        = 1000)
                
        
                