# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 13:41:23 2023

@author: huckg
"""

enzymes_TOY_1 = [
              ['HK'   , True, 'GH', ['Glucose','ATP'],['ADP','G6P'],[],[]],
              ['G6PDH', True, 'GH', ['G6P','NADP'],['6PGL','NADPH'],[],[]],
              [False  , True,False, ['6PGL'],['6PGA'],[],[]],
              ['6PGDH', True, 'GH', ['6PGA','NADP'],['RU5P','NADPH'],[],[]], 
              ['PRI'  , True, 'MM', ['RU5P'],['R5P'],[],[]],
              ['PRPPS', True, 'GH', ['R5P','ATP'],['PRPP','AMP'],[],[]],
              
              ['HK'   , True, 'GH', ['Fructose','ATP'],['ADP','F6P'],[],[]],
              ['GDH',   True, 'GH', ['GluAct','NADPH'],['NADP','Glucose'],[],[]],
              ['GPI',   True, 'MM', ['G6P'],['F6P'],[],[]],
              ['F6PK',  True, 'GH', ['F6P','ATP'],['F16BP','ADP'],[],[]],
              ['Ald',   True, 'MM', ['F16BP'],['DHAP','GlyP'],[],[]],
              ['Gc3PDH',True, 'GH', ['Gly3P','NADP'],['DHAP','NADP'],[],[]],
              ['TPI',   True, 'MM', ['DAP'],['Gly3P'],[],[]],
              ['GAPDH', True, 'GH', ['Gly3P','NADP'],['NADPH','BPgly'],[],[]],
              ['3PGK',  True, 'GH', ['BPgly','ADP'],['P3gly','ATP'],[],[]],
              ['PGM',   True, 'MM', ['P3gly'],['2Pgly'],[],[]],
              ['ENO',   True, 'MM', ['2Pgly'],['PEP'],[],[]],
              ['LDH',   True, 'GH', ['Pyruvate','NADPH'],['Lactate','NADP'],[],[]],
              
              ['PK',   True, 'GH', ['ADP','PEP'],['Pyruvate','ATP'],[],[]],
              ['PK',   True, 'GH', ['UDP','PEP'],['Pyruvate','UTP'],[],[]],
              ['PK',   True, 'GH', ['GDP','PEP'],['Pyruvate','GTP'],[],[]],
              ['UMPK', True, 'GH', ['UMP','ATP'],['UDP','ADP'],[],[]],
              ['GMPK', True, 'GH', ['GMP','ATP'],['GDP','ADP'],[],[]],
              ['AK',   True, 'GH', ['AMP','ATP'],['ADP','ADP'],[],[]],
              ['APRT', True, 'GH', ['PRPP','Adenine'],['AMP'],[],[]],
              ['UPRT', True, 'GH', ['PRPP','Uracil'],['UMP'],[],[]]
              ]

enzymes_TOY_2 = [   
                  ['Ex_glc' , False, 'MM', ['Glucose'],['Glucose_D'],[],[]],
                  ['HEX1'   , False, 'GH', ['Glucose_D','ATP'],['ADP','G6P'],[],[]],
                  ['pts'    , False, 'GH', ['PEP','Glucose_D'],['Pyruvate','G6P'],[],[]],
                  ['PGMT'   , True, 'MM',  ['G6P'],['G1P'],[],[]],
                  ['GLGC'   , False, 'GH', ['G1P','ATP'],['adpglc'],[],[]],
                  ['G6PDH2r', True, 'GH',  ['G6P','NADP'],['NADPH','6PGL'],[],[]],
                  ['PGL'    , False, 'MM',  ['6PGL',],['6PGC'],[],[]],
                  ['GND'    , False, 'GH',  ['6PGC','NADP'],['NADPH','ru5p'],[],[]],
                  ['EDD'    , False, 'MM', ['2ddg6p'],['6PGC'],[],[]],
                  ['EDA'    , False, 'GH',  ['g3p','Pyruvate'],['2ddg6p'],[],[]],
                  ['PGI'    , False, 'MM', ['G6P'],['f6p'],[],[]],

                  ['PFKA'   , False, 'GH', ['f6p','ATP'],['ADP','fdp'],[],[]],
                  ['PFKB'   , False, 'GH', ['f6p','ATP'],['ADP','fdp'],[],[]],
                  ['FBP'    , False, 'GH',  ['ADP','fdp'],['f6p','ATP'],[],[]],
                  ['FBA'    , True, 'MM',  ['fdp'],['dahp','g3p'],[],[]],
                  ['TPI'    , True, 'MM', ['g3p'],['dahp'],[],[]],
                  ['TKT2'   , True, 'GH', ['g3p','f6p'],['x5p','e4p'],[],[]],
                  
                  ['GAPD'   , True, 'GH', ['NAD','g3p'],['NADH','13dpg'],[],[]],
                  ['PGK'    , True, 'GH', ['13dpg','ADP'],['ATP','g3p'],[],[]],
                  ['PGM'    , True, 'MM', ['g3p'],['2pg'],[],[]],
                  ['ENO'    , True, 'MM', ['2pg'],['PEP'],[],[]],
    
                  ['DDPA'   , False, 'GH', ['PEP','e4p'],['2dda7p'],[],[]],
                  ['2DDA7Pt', False, 'MM', ['2dda7p'],[],[],[]],
                  ['TALA'   , True, 'GH', ['g3p','s7p'],['e4p','f6p'],[],[]],
                  ['TKT1'   , True, 'GH', ['r5p','x5p'],['g3p','s7p'],[],[]],
                  ['R5PP'   , False, 'MM', ['r5p'],['rip_D'],[],[]],
                  ['RPI'    , True, 'MM', ['ru5p'],['r5p'],[],[]],
                  ['RPE'    , True, 'MM', ['x5p'],['ru5p'],[],[]],
                  
                  ['PPS'    , False, 'GH', ['Pyruvate','ATP'],['ADP','PEP'],[],[]],
                  ['PYKA'   , False, 'GH', ['ADP','PEP'],['ATP','Pyruvate'],[],[]],
                  ['PYKf'   , False, 'GH', ['PEP','ADP'],['ATP','Pyruvate'],[],[]],
                  ['PDH'    , False, 'GH', ['Pyruvate','NAD'],['NADH','accoa'],[],[]],
                  ['LDH_D'  , True, 'GH', ['NADH','Pyruvate'],['NAD','Lac_D'],[],[]],
                  ['PFL'    , False, 'GH', ['Pyruvate','coa'],['FOR','accoa'],[],[]],
  
                  ['PTAr'     , True, 'MM', ['accoa'],['coa','actp'],[],[]],
                  ['ACALD'    , True, 'GH', ['NADH','accoa'],['NAD','acald'],[],[]],
                  ['ACS'      , False, 'MM', ['AC'],['accoa'],[],[]],
                  ['ACKr'     , True, 'GH', ['actp','ADP'],['ATP','AC'],[],[]],
                  ['ACALD2x'  , True, 'GH', ['acald','NADH'],['NAD','etoh'],[],[]],
                  ['CS'       , False, 'GH', ['accoa','oaa'],['coa','cit'],[],[]],
                
                  ['ACONTa'   , True, 'MM', ['cit'],['acon_C'],[],[]],
                  ['ACONTb'   , True, 'MM', ['acon_C'],['icit'],[],[]],
                  ['ICL'      , False, 'MM', ['icit'],['glx','succ'],[],[]],
                  ['MALS'     , False, 'GH', ['glx','accoa'],['mal_L','coa'],[],[]],
                  ['ICDhyr'   , True, 'GH', ['icit','NADP'],['akg','NADPH'],[],[]],
                  ['AKGDH'    , False, 'GH', ['akg','coa'],['succoa'],[],[]],
                  ['SACOAS'   , False, 'GH', ['succoa','ADP'],['succ','ATP','coa'],[],[]],
                  ['SUCDi'    , False, 'GH', ['succ','Q'],['fum','QH2'],[],[]],
                  ['FUM'      , True, 'MM', ['fum'],['mal_L'],[],[]],
                  ['MDH'      , True, 'GH', ['mal_L','NAD'],['NADH','oaa'],[],[]],
                  ['MDH2'     , False, 'GH', ['mal_L','Q'],['oaa','QH2'],[],[]],
                  ['FRD2'     , False, 'GH', ['fum','QH2'],['succ','Q'],[],[]],
                  ['PPC'      , True, 'MM', ['PEP'],['oaa'],[],[]],
                  ['PPCK'     , False, 'GH', ['oaa','ATP'],['ADP','PEP'],[],[]],

                  ['GLUDy'    , True, 'GH', ['akg','NADPH'],['glu_L','NADP'],[],[]],
                  ['GLUSy'    , False,'GH', ['akg','gln_L'],['glu_L'],[],[]],
                  ['GLNS'     , False, 'GH', ['glu_L','ATP'],['ADP','gln_L'],[],[]],
                  ['GLUN'     , False, 'MM', ['gln_L'],['glu_L'],[],[]],
    
                  ['ME2'      , False,'GH', ['mal_L','NADP'],['Pyruvate','NADPH'],[],[]],
                  ['ME1'      , False, 'GH', ['mal_L','NAD'],['Pyruvate','NADH'],[],[]],
    
                  ['ADK1'     , True,'GH', ['ATP','AMP'],['ADP','ADP'],[],[]],
                  ['THD2'     , False, 'GH',  ['NADH','NADP'],['NADPH'],[],[]],
                  ['ATPSrpp'  , True,'MM', ['ADP'],['ATP'],[],[]],
                  ['CYTBD'    , False, 'MM',  ['QH2'],['Q'],[],[]],
                  ['NADTRHD'  , False, 'GH',  ['NADPH','NAD'],['NADP','NADH'],[],[]],
                  ['NADH16'   , False,'GH', ['QH2','NADH'],['Q','NAD'],[],[]],
                  ['ATPM'     , False, 'MM',  ['ATP'],['ADP'],[],[]]
                  ]  

"""the libaries that are needed"""
import copy
import math
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

"""Enzyme All Enzymes model"""
network      = enzymes_TOY_2 
modelname    = 'ModelBuilderTest'   

"""Syringe load"""
syringe_load = {}

"""Get the reactions and categorize them"""
reactions = categorize_reactions(network)

"""Build the equations for the reaction and reactor kinetics"""
syringe_load = {}

"""indata and enzymes"""
indata       = []
enzymes      = []
for i in reactions:
    if i['E']:
        if i['E'] not in indata:
            indata.append(i['E'])
            syringe_load[len(syringe_load)] = [i['E']]
            enzymes.append(i['E'])
            
"""LOAD THIS PUPPY UP WITH ENERGY"""
syringe_load[len(syringe_load)] = ['PEP']
syringe_load[len(syringe_load)] = ['ATP']
syringe_load[len(syringe_load)] = ['GMP']
syringe_load[len(syringe_load)] = ['Adenine']
syringe_load[len(syringe_load)] = ['Uracil']
syringe_load[len(syringe_load)] = ['PRPP']

"""construct the models"""
reactor_kinetics,reaction_kinetics = construct_model_fluxterms(reactions,reactormix = syringe_load)   
"""Construct the model"""
basemodel = construct_kinetic_combinations(reactor_kinetics,reaction_kinetics,name = modelname,reactor_type = 'reactor_predefined_control_single_reactor')
print(basemodel)
"""Get the base model to work"""
model = basemodel.return_model()
