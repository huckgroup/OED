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

def plot_optimal_prediction(directory,measurements,name = '',mode = '',basepath = '',r2 = False):
    import pickle
    import os
    
    """the basepath """
    basepath = 'C:\\Users\\huckg\\OneDrive\\Desktop\PPP Angewandte Figures Insilico\\'
    
    """import simulation package and biuld a store folder"""
    from OptimizationOperations import simulate_optimal_parametersets
    try:
        os.mkdir(basepath + name + '\\')
    except:
        pass
    try:
        os.mkdir(basepath + name + '\\' + mode + '\\')
    except:
        pass
    path = basepath + name + '\\' + mode + '\\'
    
    """Get the listdir and the directory to import the parameters"""
    parameters = []
    """import the parameters and sort by score"""
    for i in os.listdir(directory + '\\New Parameters\\'):
        with open(directory + '\\New Parameters\\' + i, 'rb') as f:
            parameters.append(pickle.load(f))
            
    print(directory)
    """Sort parameters by score"""
    indices    = numpy.argsort([i['score'] for i in parameters])
    """Build a dataset"""
    parameters = [parameters[indices[i]] for i in range(len(parameters))] 
    print(parameters)
    # import pandas
    # import seaborn as sns
    # plt.figure(figsize = (4,10))
    # data = pandas.DataFrame(parameters)
    # sns.boxplot(data = data,orient="h")
    # plt.xscale('log')
    # plt.tight_layout()
    # plt.savefig(path + 'boxplot.png',dpi =600)  
    # plt.close()

    """Simulate the optimal parameters and store simdata"""
    data = []
    for p in parameters:
        data.append(simulate_optimal_parametersets(measurements,p))

        
    state_filter = ['AMP','ADP','NADP','ATP']
    for i in data:
        for k,v in i[0].items():
            if k in state_filter:
                plt.plot(v,label = k)
    plt.savefig(path + 'All Prediction.png',dpi =600)   
    plt.close()
    
    """scores for the parametersets"""
    scores     = []
    """tracks the scores per species i.e. goodness of fit analysis"""
    sf_species = {}

    R2 = {}
    # data = [data[10]]
    for simdata in data:
        sf = 0
        """extract the measurement (single measurements only for this analysis"""
        measurement = measurements[0]
        for state in measurement.profile.keys():
            if state in state_filter:
                if state not in R2:
                    R2[state] = []
                if state not in sf_species:
                    sf_species[state] = []
                time,concentration = zip(*measurement.profile[state])
                m_max = numpy.amax(concentration)
                m_min = numpy.amin(concentration)
                
                """iterate through timepoints"""
                k,s = 0,0
                for i in time:
                    try:
                        # splus = ( ( ( simdata[0][state][i] - m_min) / (m_max - m_min) ) - ( ( concentration[k] - m_min ) / (m_max - m_min) ) ) **2
                        splus =  numpy.sqrt((concentration[k] - simdata[0][state][i])**2)
                        if splus != math.inf:
                            R2[state].append((concentration[k],simdata[0][state][i]))
                            s += splus
                    except:
                        print('species not accounted for')
                    k += 1
                s = s/len(time)
                factor = 1
                if state == 'AMP':
                    factor = 2
                if state == 'ADP':
                    factor = 2
                s = s * factor
                sf_species[state].append(s)
                sf += copy.copy(s)
        scores.append(sf)
        
    """calculate mean_sf"""
    mean_sf_species = {k:numpy.mean(v) for k,v in sf_species.items()}
    order  = sorted(mean_sf_species.items(), key=lambda item: item[1], reverse=True)
    """get some colors to species"""
    c = ['orange', 'tab:blue', 'green', 'red', 'pink', 'turquoise', 'purple', 'y']
    """assign colors"""
    color = {order[i][0]:c[i] for i in range(len(order))}

    
    import pandas
    import seaborn as sns
    df = pandas.DataFrame(sf_species)
    plt.figure(figsize = (3,6))
    sns.barplot(data = df,edgecolor = 'k',orient = 'h',palette = color,order = [n for n,y in order])
    plt.gca().set_facecolor('whitesmoke')
    plt.xlabel('Goodness of fit',size = 10)
    plt.tight_layout()
    plt.savefig(path + 'Goodness of fit.png',dpi = 600)
    plt.savefig(path + 'Goodness of fit.svg',dpi = 600)
    plt.close()
    
    """get the order of the simulations"""
    rank = [numpy.argsort(scores)[i] for i in range(3)]
    highest_rank = print(rank[0],'THIS IS THE HIGHEST RANKED BLAALBLBALBLABLBALBA')
    data = [data[i] for i in range(5)]

    # """Plot the interpolated raw data Figure"""
    # plt.figure(figsize = (5,5))
    # plt.gca().set_facecolor('whitesmoke')
    # for number, measurement in measurements.items():
    #     for state in sorted(measurement.profile.keys()):
    #         if state in state_filter:
    #             time,concentration = zip(*measurement.profile[state])
    #             df = pandas.DataFrame({state:concentration,'Time (min)':time})
    #             df['Concentration (uM)'] = df[state].rolling(window=10, min_periods=1).mean()
    #             sns.lineplot(data = df,x = 'Time (min)',y = 'Concentration (uM)',errorbar= ('sd',1),color = color[state],label = state)
    #     plt.legend(fancybox = True)
    #     plt.savefig(path + 'Original data.png',dpi =600)   
    #     plt.savefig(path + 'Original data.svg',dpi =600)   
    #     plt.close()
        
    fullplotdata = {}
    for number, measurement in measurements.items():
        n = 1
        for state in sorted(measurement.profile.keys()):
            if state in state_filter:
                if state not in fullplotdata:
                    fullplotdata[state] = []
                plt.figure(figsize = (3,3))
                time,concentration = zip(*measurement.profile[state])
                plt.scatter([time[i] for i in range(0,len(time),1)],[concentration[i] for i in range(0,len(concentration),1)],s = 4,marker = '^',c = 'k',label = state)
                
                'import the libraries'
                import seaborn as sns
                import pandas
                
                """the plotdata and the time nad concentration"""
                plotdata = {'Time (min)':[],'Concentration (uM)':[]}
                for simdata in data:
                    statedata = simdata[number]
                    plotdata['Concentration (uM)'].extend([k for k in statedata[state]])
                    plotdata['Time (min)'].extend(range(len([k for k in statedata[state]])))

                df = pandas.DataFrame(plotdata)
                fullplotdata[state] = df
                sns.lineplot(data = df,x = 'Time (min)',y = 'Concentration (uM)',errorbar= ('sd',2),color = color[state])
                sns.lineplot(data = df,x = 'Time (min)',y = 'Concentration (uM)',errorbar= ('sd',1),color = color[state])
                plt.legend(fancybox = True,loc = 'upper right')
                plt.tight_layout()
                lk = ''
                if len(measurements) > 1:
                    lk = 'measurement'
                """Save the figure"""
                plt.savefig(path + 'Prediction {}.png'.format(state),dpi =600)   
                # plt.savefig(path + 'Prediction {}.svg'.format(state),dpi =600)       
                plt.close()

    r_species = {}
    for k,values in R2.items():
        x,y = [],[]
        d,p = zip(*values)
        for i in range(len(d)):
                if d[i] != math.isnan and p[i] != math.isnan:
                    if d[i] != math.inf and p[i] != math.inf:
                        x.append(d[i])
                        y.append(p[i])      
        plt.scatter(d,p,color = color[k],alpha = 0.25,label = k)
        r2,line1 = calculate_R(numpy.array(x),numpy.array(y))
        plt.plot(numpy.array(sorted(x)),line1,color = 'DarkRed')
        plt.gca().set_facecolor('whitesmoke')
        plt.ylabel('Model Prediction (Lognorm uM)',size = 15)
        plt.xlabel('data (Lognorm uM)',size = 15)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title(k + '| ' + str(round(r2,3)))
        plt.tight_layout()
        plt.legend(fancybox = True)
        plt.savefig(path + 'R_test_{}.png'.format(k),dpi = 600)
        plt.close()
        
        """get the R-species"""
        r_species[k] = copy.copy(r2)

    c = {'ATP':'slategray','NADP':'peachpuff','ADP':'royalblue','AMP':'gold'}
    x,y = [],[]
    for k,values in R2.items():
        d,p = zip(*values)
        for i in range(len(d)):
                if d[i] != math.isnan and p[i] != math.isnan:
                    if d[i] != math.inf and p[i] != math.inf:
                        x.append(d[i])
                        y.append(p[i])      
        plt.scatter(d,p,color = c[k],alpha = 0.25,label = k)
    r2,line1 = calculate_R(numpy.array(x),numpy.array(y))
    plt.plot(numpy.array(sorted(x)),line1,color = 'DarkRed')
    plt.ylabel('Model Prediction',size = 15)
    plt.xlabel('Data',size = 15)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(str(round(r2,3)))
    plt.tight_layout()
    plt.legend(fancybox = True)
    plt.savefig(path + 'R_test.png',dpi = 600)
    plt.close()    
    
    """scores for the parametersets"""
    scores     = []
    """tracks the scores per species i.e. goodness of fit analysis"""
    error = {}

    R2 = {}
    for simdata in data:
        sf = 0
        """extract the measurement (single measurements only for this analysis"""
        measurement = measurements[0]
        for state in measurement.profile.keys():
            if state in state_filter:
                if state not in R2:
                    R2[state] = []
                if state not in error:
                    error[state] = []
                time,concentration = zip(*measurement.profile[state])
                m_max = numpy.amax(concentration)
                m_min = numpy.amin(concentration)
                
                """iterate through timepoints"""
                k,s = 0,0
                for i in time:
                    try:
                        deviation = concentration[k] / simdata[0][state][i]
                        if deviation > 1:
                            if deviation < 8:
                                splus = deviation 
                        else:
                            ratio = 1./deviation
                            if ratio < 8:
                                splus = (1./deviation)
                        if splus != math.inf:
                            R2[state].append((concentration[k],simdata[0][state][i]))
                            s += splus
                            print(splus)
                    
                    except:
                        print('species not accounted for')
                    k += 1
                s = s/len(time)
                factor = 1
                if state == 'AMP':
                    factor = 9
                if state == 'ADP':
                    factor = 4
                # s = s * factor
                error[state].append((s-1)*100)
                sf += copy.copy((s-1)*100)
        scores.append(sf)    
    
    return r_species,fullplotdata, error

def plot_optimal_r2(directory,measurements,name = '',mode = '',basepath = '',r2 = False):
    import pickle
    import os
    
    """the basepath """
    basepath = 'C:\\Users\\huckg\\OneDrive\\Desktop\PPP Angewandte Figures Insilico\\'
    
    """import simulation package and biuld a store folder"""
    from OptimizationOperations import simulate_optimal_parametersets
    try:
        os.mkdir(basepath + name + '\\')
    except:
        pass
    try:
        os.mkdir(basepath + name + '\\' + mode + '\\')
    except:
        pass
    path = basepath + name + '\\' + mode + '\\'
    
    """Get the listdir and the directory to import the parameters"""
    parameters = []
    """import the parameters and sort by score"""
    for i in os.listdir(directory + '\\New Parameters\\'):
        with open(directory + '\\New Parameters\\' + i, 'rb') as f:
            parameters.append(pickle.load(f))
            
    """Sort parameters by score"""
    indices    = numpy.argsort([i['score'] for i in parameters])
    """Build a dataset"""
    parameters = [parameters[indices[i]] for i in range(len(parameters))] 

    """Simulate the optimal parameters and store simdata"""
    data = []
    for p in parameters:
        data.append(simulate_optimal_parametersets(measurements,p))

    """scores for the parametersets"""
    scores     = []
    """tracks the scores per species i.e. goodness of fit analysis"""
    sf_species = {}
    state_filter = ['AMP','ADP',"ATP","NADP"]
    R2 = {}
    sda = [data[10]]
    print(sda)
    for simdata in sda:
        sf = 0
        """extract the measurement (single measurements only for this analysis"""
        measurement = measurements[0]
        for state in measurement.profile.keys():
            if state in state_filter:
                if state not in R2:
                    R2[state] = []
                if state not in sf_species:
                    sf_species[state] = []
                time,concentration = zip(*measurement.profile[state])
                m_max = numpy.amax(concentration)
                m_min = numpy.amin(concentration)
                
                """iterate through timepoints"""
                k,s = 0,0
                for i in time:
                    try:
                        # splus = ( ( ( simdata[0][state][i] - m_min) / (m_max - m_min) ) - ( ( concentration[k] - m_min ) / (m_max - m_min) ) ) **2
                        splus =  numpy.sqrt((concentration[k] - simdata[0][state][i])**2)
                        if splus != math.inf:
                            R2[state].append((concentration[k],simdata[0][state][i]))
                            s += splus
                    except:
                        print('species not accounted for')
                    k += 1
                s = s/len(time)
                factor = 1
                if state == 'AMP':
                    factor = 2
                if state == 'ADP':
                    factor = 2
                s = s * 1
                sf_species[state].append(s)
                sf += copy.copy(s)
        scores.append(sf)
        
    """calculate mean_sf"""
    mean_sf_species = {k:numpy.mean(v) for k,v in sf_species.items()}
    order  = sorted(mean_sf_species.items(), key=lambda item: item[1], reverse=True)
    """get some colors to species"""

    print('dsablakjbdfslaksbdflkasjdbflkajbfdlakdsbfalksdjbflaskjdbflakdjbf')
    c = {'ATP':'slategray','NADP':'peachpuff','ADP':'royalblue','AMP':'gold'}
    x,y = [],[]
    
    plt.figure(figsize = (3,3))
    ax = plt.subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for k,values in R2.items():
        d,p = zip(*values)
        for i in range(len(d)):
                if d[i] != math.isnan and p[i] != math.isnan:
                    if d[i] != math.inf and p[i] != math.inf:
                        x.append(d[i])
                        y.append(p[i])      
        plt.scatter(d,p,color = c[k],alpha = 0.25,label = k)
    r2,line1 = calculate_R(numpy.array(x),numpy.array(y))
    plt.plot(numpy.array(sorted(x)),line1,color = 'DarkRed')
    plt.ylabel('Model Prediction (uM)',size = 12)
    plt.xlabel('Data (uM)',size = 12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # plt.title(str(round(r2,3)))
    plt.tight_layout()
    # plt.legend(fancybox = True)
    plt.savefig('C:\\Users\\huckg\\OneDrive\\Desktop\\Angewandte Paper Figures\\' + 'R_test.png',dpi = 600)
    plt.close()    
    


 ############################ --- TEST ONE --- ######################### 
"""Obtain the experimental data and time dependent inputs"""
path = "C:\\Users\\huckg\\OneDrive\\Desktop\\Projects\\-active- Free Enzyme OED\\PPP pathway\\Simulated Results\\Test data\\"
observables,time_dependent_inputs,known_parameters_TEST = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load_TEST,stock_test = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')

"""delete NADPH"""
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
if 'NADPH' in measurement.profile:
    del measurement.profile['NADPH']
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
total = sum(list(TDI[list(sorted(list(TDI.keys())))[0]].values()))
inverted_syringe_load_test = {i:'kflow(channel_{})'.format(k) for k,v in syringe_load_TEST.items() for i in v}

T = {}
T[(len(T))] = copy.copy(measurement)

 ############################ --- MEASUREMENT ONE --- ######################### 
"""Obtain the experimental data and time dependent inputs"""
path = "C:\\Users\\huckg\\OneDrive\\Desktop\\Projects\\-active- Free Enzyme OED\\PPP pathway\\Simulated Results\\Free OED\\"
observables,time_dependent_inputs,known_parameters_M1 = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load,stock = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')


"""delete NADPH"""
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
if 'NADPH' in measurement.profile:
    del measurement.profile['NADPH']
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
total = sum(list(TDI[list(sorted(list(TDI.keys())))[0]].values()))
inverted_syringe_load = {i:'kflow(channel_{})'.format(k) for k,v in syringe_load.items() for i in v}

# measurement.show(show = True)
M = {}
M[(len(M))] = copy.copy(measurement)

"""Define the model of the system"""
 ############################ --- MEASUREMENT ONE --- #################### 
"""Obtain the experimental data and time dependent inputs"""
path = "C:\\Users\\huckg\\OneDrive\\Desktop\\Projects\\-active- Free Enzyme OED\\PPP pathway\\Simulated Results\\Free SS\\"
observables,time_dependent_inputs,known_parameters_M1 = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load,stock = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')


"""delete NADPH"""
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
if 'NADPH' in measurement.profile:
    del measurement.profile['NADPH']
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
total = sum(list(TDI[list(sorted(list(TDI.keys())))[0]].values()))
inverted_syringe_load = {i:'kflow(channel_{})'.format(k) for k,v in syringe_load.items() for i in v}
# measurement.show(show = True)
M[(len(M))] = copy.copy(measurement)

"""Define the model of the system"""
 ########################## --- MEASUREMENT ONE --- ########################### 
"""Obtain the experimental data and time dependent inputs"""
path = "C:\\Users\\huckg\\OneDrive\\Desktop\\Projects\\-active- Free Enzyme OED\\PPP pathway\\Simulated Results\\Beads OED\\"
observables,time_dependent_inputs,known_parameters_M1 = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load,stock = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')

"""delete NADPH"""
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
if 'NADPH' in measurement.profile:
    del measurement.profile['NADPH']
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
total = sum(list(TDI[list(sorted(list(TDI.keys())))[0]].values()))
inverted_syringe_load = {i:'kflow(channel_{})'.format(k) for k,v in syringe_load.items() for i in v}
M[(len(M))] = copy.copy(measurement)

######################### --- MEASUREMENT ONE --- ############################# 
"""Obtain the experimental data and time dependent inputs"""
path = "C:\\Users\\huckg\\OneDrive\\Desktop\\Projects\\-active- Free Enzyme OED\\PPP pathway\\Simulated Results\\Beads SS\\"
observables,time_dependent_inputs,known_parameters_M1 = parse_flow_data(path)
"""Get the data and plug it in a database"""
TDI,syringe_load,stock = translate_time_dependent_inputs(time_dependent_inputs)
"""Plug these vectors into a measurement module"""
TDI = translate_TDI(TDI)
"""Measurement and measurement object"""
measurement = MeasurementObject(
                  observables,
                  time_dependent_parameters = TDI,
                  time_unit = 'min')

"""Delete NADPH"""
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
if 'NADPH' in measurement.profile:
    del measurement.profile['NADPH']
measurement.observables = [i for i in measurement.observables if i != 'NADPH']
total = sum(list(TDI[list(sorted(list(TDI.keys())))[0]].values()))
inverted_syringe_load = {i:'kflow(channel_{})'.format(k) for k,v in syringe_load.items() for i in v}
M[(len(M))] = copy.copy(measurement)

"""Define the model of the system"""
###############################################################################
from random import shuffle    
import os
        
modelname = 'PPP_Literature_inh_EO_rev'            
network  =    [
              ['HK'   , True, 'EO', ['Glucose','ATP'],['ADP','G6P'],[],[]],
              ['G6PDH', True, 'EO', ['NADP','G6P'],['6PGL','NADPH'],[],[]],
              [False  , True,False, ['6PGL'],['6PGA'],[],[]],
              ['6PGDH', True, 'EO', ['NADP','6PGA'],['RU5P','NADPH'],[],[]],  #NADPH, ATP,
              ['PRI'  , True, 'MM', ['RU5P'],['R5P'],[],[]],
              ['PRPPS', True, 'EO', ['ATP','R5P'],['PRPP','AMP'],[],['ADP']]
              ]

"""Get the reactions and categorize them"""
reactions = categorize_reactions(network)
"""Build the equations for the reaction and reactor kinetics"""
reactor_kinetics,reaction_kinetics = construct_model_fluxterms(reactions,reactormix = syringe_load)    
"""Construct the model"""
basemodel = construct_kinetic_combinations(reactor_kinetics,reaction_kinetics,name = modelname)
"""directory where we store the results"""
directory = "C:\\Users\\huckg\\OneDrive\\Desktop\\PPP Angewandte Final\\{0} {1}\\".format(modelname,'M1') 
"""return the base model"""
model = basemodel.return_model(fixed = known_parameters_M1)

from random import shuffle    
import os
        
modelname = 'PPP_Literature_inh_EO_rev'            
network  =    [
              ['HK'   , True, 'EO', ['Glucose','ATP'],['ADP','G6P'],[],[]],
              ['G6PDH', True, 'EO', ['NADP','G6P'],['6PGL','NADPH'],[],[]],
              [False  , True,False, ['6PGL'],['6PGA'],[],[]],
              ['6PGDH', True, 'EO', ['NADP','6PGA'],['RU5P','NADPH'],[],[]],  #NADPH, ATP,
              ['PRI'  , True, 'MM', ['RU5P'],['R5P'],[],[]],  
              ['PRPPS', True, 'EO', ['ATP','R5P'],['PRPP','AMP'],[],['ADP']]
              ]

"""Get the reactions and categorize them"""
reactions = categorize_reactions(network)
"""Build the equations for the reaction and reactor kinetics"""
reactor_kinetics,reaction_kinetics = construct_model_fluxterms(reactions,reactormix = syringe_load_TEST)    
"""Construct the model"""
basemodel = construct_kinetic_combinations(reactor_kinetics,reaction_kinetics,name = modelname)
"""directory where we store the results"""
directory = "C:\\Users\\huckg\\OneDrive\\Desktop\\PPP Angewandte\\{0} {1}\\".format(modelname,'M1') 
"""return the base model"""
model_test = basemodel.return_model(fixed = known_parameters_TEST)

import os
import pickle
import random
"""Get the listdir and the directory to import the parameters"""
parameters = []
"""import the parameters and sort by score"""
for i in os.listdir(directory + '\\New Parameters\\'):
    with open(directory + '\\New Parameters\\' + i, 'rb') as f:
        parameters.append(pickle.load(f))
scores = numpy.argsort([i['score'] for i in parameters])
chosen = parameters[scores[5]]  
# for k,v in chosen.items():
#     if 'k_cat' in k:
#         chosen[k] = v*7
#     if 'k_N_inh(PRPPS:ADP)' in k:
#         chosen[k]  = v*7
#         print(chosen[k])
# enzymes = ['PRPPS(in)','PRI(in)']
# for i in model.fixed.keys():
#     if i in enzymes:
#         model.fixed[i] = model.fixed[i] * 10
# model.fixed['HK(in)'] = model.fixed['HK(in)']*2
for k,v in chosen.items():
    if 'k_cat' in k:
        chosen[k] = v*3000
    if 'k_N_inh(PRPPS:ADP)' == k:
        chosen[k]  = 8
print(chosen ,'sadfasdfasdfasdfasdfasfd')
enzymes = ['PRPPS(in)','PRI(in)']
for i in model.fixed.keys():
    if i in enzymes:
        model.fixed[i] = model.fixed[i] * 25
model.fixed['HK(in)'] = model.fixed['HK(in)']*25

"""delete the score of the fit"""
del chosen['score']   
model_test.observables = ['ATP','ADP','AMP','NADP']
model.observables = ['ATP','ADP','AMP','NADP']

"""Prune the measurements"""
###############################################################################
from Measurements import simulate_measurement
measurement.model_test = model_test 
OED = T[0]
TEST = simulate_measurement(model_test,time_dependent_parameters = OED.time_dependent_parameters,parameters = chosen)
"""conditions updated from inflow stock"""
TEST.conditions = {state:value * TEST.time_dependent_parameters[list(sorted(list(TEST.time_dependent_parameters.keys())))[0]][inverted_syringe_load_test[state]]/total for state,value in stock_test.items()}
for state,data in TEST.profile.items():
    TEST.conditions[state] = data[0][1]
TEST.model = model_test

"""Prune the measurements"""
###############################################################################
from Measurements import simulate_measurement
measurement.model = model 
OED = M[0]
m1 = simulate_measurement(model,time_dependent_parameters = OED.time_dependent_parameters,parameters = chosen)
"""conditions updated from inflow stock"""
m1.conditions = {state:value * m1.time_dependent_parameters[list(sorted(list(m1.time_dependent_parameters.keys())))[0]][inverted_syringe_load[state]]/total for state,value in stock.items()}
for state,data in m1.profile.items():
        m1.conditions[state] = data[0][1]
        
###################################################################
OED_SS = M[1] 
m2 = simulate_measurement(model,time_dependent_parameters = OED_SS.time_dependent_parameters,parameters = chosen)
for state,vector in m2.profile.items():
    new_vector = []
    for time,value in vector:
        if int(time)%150 == 0:
            new_vector.append((time,value))
    m2.profile[state] = copy.copy(new_vector)
"""conditions updated from inflow stock"""
m2.conditions = {state:value * m2.time_dependent_parameters[list(sorted(list(m2.time_dependent_parameters.keys())))[0]][inverted_syringe_load[state]]/total for state,value in stock.items()}
for state,data in m2.profile.items():
        m2.conditions[state] = data[0][1]
    ###############################################################################       
Beads_OED = M[2]
m3 = simulate_measurement(model,time_dependent_parameters = Beads_OED.time_dependent_parameters,parameters = chosen)

prune = (0,300,'ATP')
m3.prune(prune)
prune = (0,300,'ADP')
m3.prune(prune)
prune = (0,300,'AMP')
m3.prune(prune)
prune = (0,300,'NADP')
m3.prune(prune)   
"""conditions updated from inflow stock"""
m3.conditions = {state:value * m3.time_dependent_parameters[list(sorted(list(m3.time_dependent_parameters.keys())))[0]][inverted_syringe_load[state]]/total for state,value in stock.items()}
for state,data in m3.profile.items():
        m3.conditions[state] = data[0][1]
#print(m3.time_dependent_parameters)
#m3.show(show = True)

###############################################################################
Beads_OED_SS = M[3] 
m4 = simulate_measurement(model,time_dependent_parameters = Beads_OED_SS.time_dependent_parameters,parameters = chosen)
for state,vector in m4.profile.items():
    new_vector = []
    for time,value in vector:
        if int(time)%150 == 0:
            new_vector.append((time,value))
    m4.profile[state] = copy.copy(new_vector)
    
prune = (0,300,'ATP')
m4.prune(prune)
prune = (0,300,'ADP')
m4.prune(prune)
prune = (0,300,'AMP')
m4.prune(prune)
prune = (0,300,'NADP')
m4.prune(prune)  

"""conditions updated from inflow stock"""
m4.conditions = {state:value * m4.time_dependent_parameters[list(sorted(list(m4.time_dependent_parameters.keys())))[0]][inverted_syringe_load[state]]/total for state,value in stock.items()}
for state,data in m4.profile.items():
        m4.conditions[state] = data[0][1]  
        
"""show the measurements"""
# m1.show(show = True)   
# m2.show(show = True)  
# m3.show(show = True)   
# m4.show(show = True)  

model.observables = ['ATP','ADP','AMP','NADP']
"""measurements"""
m1.model = model
m2.model = model
m3.model = model
m4.model = model

def plot_insilico_data(profile, step = False,name = ''):
    colors = {'ATP':(38/255.,70/255.,83/255.),
            'AMP':(233/255.,196/255.,106/255.),
            'NADP':(244/255.,162/255.,97/255.),
            'ADP':(42/255.,157/255.,83/255.)
            }
    
    
    plt.figure(figsize = (3.9,3.9))
    ax = plt.subplot(111)
    # Hide the right and top spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for k in ['AMP','ADP','ATP','NADP']:
        v = profile[k]

        time,values  = zip(*v)
        if name == 'OED Bead':
            print('asdnh;lakjhdsflakjshdflajhdsflajhdsf')
            time = [i-300 for i in time]
            print(time)
        if step == False:
            ax.plot(time,values,label = k,color = colors[k],linewidth=3)
        else:
            ax.step(time,values,label = k,color = colors[k],linewidth = 2)
            ax.scatter(time,values,color = 'k', marker = '^')

    # plt.gca().set_facecolor((1,253/255.,239/255.))
    # plt.legend(fancybox = True,loc = 'upper left',fontsize  = 15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=15)
    plt.xlabel('Time (min)',size = 16)
    plt.ylabel('Concentration (uM)',size = 16)
    plt.tight_layout()
    plt.savefig('C:\\Users\\huckg\\OneDrive\\Desktop\\Angewandte Paper Figures\\{}.png'.format(name),dpi = 800)
    plt.close()
    
plot_insilico_data(m1.profile,name = 'OED')
plot_insilico_data(m2.profile,step = True,name = 'OED SS')
plot_insilico_data(m3.profile,name = 'OED Bead')
plot_insilico_data(m4.profile,step = True,name = 'OED Bead SS')
    
from IdentifiabilityAnalysisModelSets import IdentifiabilityAnalysis
IDF = IdentifiabilityAnalysis({0:m1},split = 1,nm = 'OED',observables = model.observables)
IDF.SVD()   
IDF = IdentifiabilityAnalysis({0:m2},split = 150,nm = 'SS',observables = model.observables)
IDF.SVD()   
IDF = IdentifiabilityAnalysis({0:m3},split = 1,nm = 'OED Beads',observables = model.observables)
IDF.SVD()   
IDF = IdentifiabilityAnalysis({0:m4},split = 150,nm = 'Beads SS',observables = model.observables)
IDF.SVD()  

import random
for k,v in model.fixed.items():
    if k not in known_parameters_M1.keys():
        model.fixed[k] = v*random.random()
model.include = [i for i in model.include if i not in known_parameters_M1.keys()]
m1.model = model
m2.model = model
m3.model = model
m4.model = model
for optimize in range(10):
    i = 1
    for prune_time in reversed([15,30,60,120,180,240,420]):
        prune = (prune_time,1000,'ATP')
        m1.prune(prune)
        prune = (prune_time,1000,'ADP')
        m1.prune(prune)
        prune = (prune_time,1000,'AMP')
        m1.prune(prune)
        prune = (prune_time,1000,'NADP')
        m1.prune(prune) 
        BoxplotOptimization({0:m1},
                          include             = model.include,
                          optimization_number = 1,
                          generations         = 300,
                          agents              = 35,
                          storedata           = 'C:/Users/huckg/OneDrive/Desktop/PPP Angewandte In Silico Final Paper/' + 'm' + str(1) + ' ' + '{0}'.format(prune_time),
                          startsamples        = 1000)
        
        
    i += 1
    for prune_time in reversed([151,301,601,1201,1801,2401,4201]):
        prune = (prune_time,4000,'ATP')
        m2.prune(prune)
        prune = (prune_time,4000,'ADP')
        m2.prune(prune)
        prune = (prune_time,4000,'AMP')
        m2.prune(prune)
        prune = (prune_time,4000,'NADP')
        m2.prune(prune)  
        BoxplotOptimization({0:m2},
                          include             = model.include,
                          optimization_number = 1,
                          generations         = 300,
                          agents              = 35,
                          storedata           = 'C:/Users/huckg/OneDrive/Desktop/PPP Angewandte In Silico Final Paper/' + 'm' + str(2) + ' ' + '{0}'.format(prune_time),
                          startsamples        = 1000)
        
    i += 1
    for prune_time in reversed([15,30,60,120,180,240,360,420]):
        prune = (300 + prune_time,1000,'ATP')
        m3.prune(prune)
        prune = (300 + prune_time,1000,'ADP')
        m3.prune(prune)
        prune = (300 + prune_time,1000,'AMP')
        m3.prune(prune)
        prune = (300 + prune_time,1000,'NADP')
        m3.prune(prune)    
        BoxplotOptimization({0:m3},
                          include             = model.include,
                          optimization_number = 1,
                          generations         = 400,
                          agents              = 35,
                          storedata           = 'C:/Users/huckg/OneDrive/Desktop/PPP Angewandte In Silico Final Paper/' + 'm' + str(3) + ' ' + '{0}'.format(prune_time),
                          startsamples        = 1000)
    i += 1
    for prune_time in reversed([151,301,601,1201,1801,2401,3001,4201]):
        prune = (300 + prune_time,4000,'ATP')
        m4.prune(prune)
        prune = (300 + prune_time,4000,'ADP')
        m4.prune(prune)
        prune = (300 + prune_time,4000,'AMP')
        m4.prune(prune)
        prune = (300 + prune_time,4000,'NADP')
        m4.prune(prune)    
        BoxplotOptimization({0:m4},
                          include             = model.include,
                          optimization_number = 1,
                          generations         = 400,
                          agents              = 35,
                          storedata           = 'C:/Users/huckg/OneDrive/Desktop/PPP Angewandte In Silico Final Paper/' + 'm' + str(4) + ' ' + '{0}'.format(prune_time),
                          startsamples        = 1000)


##############################################################################
###############################################################################
###############################################################################
###############################################################################
import os
def process_and_plot_data(name,timelist,TEST,r2 = False):
    path              = 'C:\\Users\\huckg\\OneDrive\\Desktop\\PPP Angewandte In Silico Final Paper\\'
    store             = 'C:\\Users\\huckg\\OneDrive\\Desktop\\\\Angewandte Paper Figures\\'
    folder            = os.listdir(path)
    dataset           = []
    R_prediction_data = {}
    r_squared         = []
    error_data       = []
    
    import numpy as np
    import seaborn as sns
    
    for i in reversed(folder):
        try:
            experiment = i.split(' ')[0]
            time = i.split(' ')[1]
            if name == experiment:
                if int(time) in timelist:
                    print(path + i)
                    r, data, error = plot_optimal_prediction(path + i + '\\', {0: TEST}, name=i, mode='TEST',r2 = r2)
                    # if str(time) == str(420):
                    #    plot_optimal_r2(path + i + '\\', {0: TEST}, name=i, mode='TEST',r2 = r2)
                    r_squared.append((time, r))
                    error_data.append((time,error))
                    dataset.append((int(time), data))
                    if experiment not in R_prediction_data:
                        R_prediction_data[experiment] = []
                    R_prediction_data[experiment].append((time, r))
        except:
            pass

    plotdata = []
    vardata  = []
    plotime  = []
    data = sorted(error_data)
    plt.figure(figsize = (2.8,2.8))
    ax = plt.subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for time, r in data:
        print(r['AMP'])
        rfinal = numpy.mean(r['AMP'])
        plotdata.append(rfinal)
        plotime.append(int(time))
    # plt.gca().set_facecolor('whitesmoke')
    arg = np.argsort(plotime)
    plt.scatter([plotime[i] for i in arg], [plotdata[i] for i in arg])
    plt.plot([plotime[i] for i in arg], [plotdata[i] for i in arg], color='navy')
    # plt.errorbar([plotime[i] for i in arg], [plotdata[i] for i in arg],yerr = vardata, color='navy')
    plt.xlabel('Training time (min)', size=12)
    plt.ylabel('Prediction accuracy (%)', size=12)
    plt.tight_layout()
    plt.savefig(f'{store}ERROR_improvement_{name}.png', dpi=800)
    plt.close()
    
    c = {'ATP':'slategray','NADP':'peachpuff','ADP':'royalblue','AMP':'gold'}
    plt.figure(figsize = (3,1.6))
    ax = plt.subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    colors =  {180:'green',
              30:(250/255.,229/255.,134/255.),
              240:'turquoise',
              420:(175/255.,134/255.,44/255.) }
    
    # colors = {'ATP':(38/255.,70/255.,83/255.),
    #         'AMP':(233/255.,196/255.,106/255.),
    #         'NADP':(244/255.,162/255.,97/255.),
    #         'ADP':(42/255.,157/255.,83/255.)
    #         }
        
    
    
    for time, data in sorted(dataset):
        df = data['AMP']
        sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
        # sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
    time, values = zip(*TEST.profile['AMP'])
    plt.scatter(time, values, marker='^', color='k',label = 'Data (test)',s = 4)
    plt.xlabel('Time (min)', size=8)
    plt.ylabel('Concentration (uM)', size=8)
    plt.legend('',frameon=False)
    plt.ylim([0,100])
    plt.tight_layout()
    plt.savefig(f'{store}Training_time_AMP_improvement_{name}.png', dpi=800)
    plt.close()
    
    plt.figure(figsize = (3,1.6))
    ax = plt.subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    colors =  {180:'green',
              420:(32/255.,132/255.,63/255.),
              240:'turquoise',
              30:(90/255.,215/255.,135/255.) }  
    
    # colors =  {180:'green',
    #           30:(250/255.,229/255.,134/255.),
    #           240:'turquoise',
    #           420:(175/255.,134/255.,44/255.) }
    
    # colors = {'ATP':(38/255.,70/255.,83/255.),
    #         'AMP':(233/255.,196/255.,106/255.),
    #         'NADP':(244/255.,162/255.,97/255.),
    #         'ADP':(42/255.,157/255.,83/255.)
    #         }
        
    for time, data in sorted(dataset):
        df = data['ADP']
        sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
        # sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
    time, values = zip(*TEST.profile['ADP'])
    plt.scatter(time, values, marker='^', color='k',label = 'Data (test)',s = 3)
    plt.xlabel('Time (min)', size=8)
    plt.ylabel('Concentration (uM)', size=8)
    plt.legend('',frameon=False)
    plt.ylim([0,400])
    plt.tight_layout()
    plt.savefig(f'{store}Training_time_ADP_improvement_{name}.png', dpi=800)
    plt.close()

    plt.figure(figsize = (3,1.6))
    ax = plt.subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    colors =  {180:'green',
              30:'yellow',
              240:'turquoise',
              420:'royalblue' } 

    colors =  {180:'green',
              30:(75/255.,105/255.,130/255.),
              240:'turquoise',
              420:(24/255.,55/255.,73/255.) }
    
    # colors = {'ATP':(38/255.,70/255.,83/255.),
    #         'AMP':(233/255.,196/255.,106/255.),
    #         'NADP':(244/255.,162/255.,97/255.),
    #         'ADP':(42/255.,157/255.,83/255.)
    #         }
            
    for time, data in sorted(dataset):
        df = data['ATP']
        sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
        # sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
    time, values = zip(*TEST.profile['ATP'])
    plt.scatter(time, values, marker='^', color='k',label = 'Data (test)',s = 3)
    plt.xlabel('Time (min)', size=8)
    plt.ylabel('Concentration (uM)', size=8)
    plt.legend('',frameon=False)
    plt.ylim([0,450])
    plt.tight_layout()
    plt.savefig(f'{store}Training_time_ATP_improvement_{name}.png', dpi=800)
    plt.close()  

    plt.figure(figsize = (3,1.6))
    ax = plt.subplot(111)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    colors = {180:'green',
              30:'yellow',
              240:'turquoise',
              420:'royalblue' }  
    

    colors =  {180:'green',
              30:(180/255.,102/255.,37/255.),
              240:'turquoise',
              420:(253/255.,175/255.,113/255.) }
    
    # colors = {'ATP':(38/255.,70/255.,83/255.),
    #         'AMP':(233/255.,196/255.,106/255.),
    #         'NADP':(244/255.,162/255.,97/255.),
    #         'ADP':(42/255.,157/255.,83/255.)
    #         }    
    for time, data in sorted(dataset):
        df = data['NADP']
        sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
        # sns.lineplot(data=df, x='Time (min)', y='Concentration (uM)', errorbar=('sd', 2),color = colors[time])
    time, values = zip(*TEST.profile['NADP'])
    plt.scatter(time, values, marker='^', color='k',label = 'Data (test)',s = 3)
    plt.xlabel('Time (min)', size=8)
    plt.ylabel('Concentration (uM)', size=8)
    plt.ylim([0,320])
    plt.legend('',frameon=False)
    plt.tight_layout()
    plt.savefig(f'{store}Training_time_NADP_improvement_{name}.png', dpi=800)
    plt.close()
    
    
    
    
    plotdata = []
    plotime  = []
    data = sorted(r_squared)
    for time, r in data:
        rfinal = r['AMP']
        plotdata.append(rfinal)
        plotime.append(int(time))
    # plt.gca().set_facecolor('whitesmoke')
    arg = np.argsort(plotime)
    plt.scatter([plotime[i] for i in arg], [plotdata[i] for i in arg])
    plt.plot([plotime[i] for i in arg], [plotdata[i] for i in arg], color='navy')
    plt.xlabel('Training time (min)', size=16)
    plt.ylabel('$R^2$', size=16)
    plt.savefig(f'{store}R2_improvement_{name}.png', dpi=800)
    plt.close()
    
    plotdata = []
    vardata  = []
    plotime  = []
    data = sorted(error_data)
    for time, r in data:
        rfinal = numpy.mean(r['AMP'])
        vfinal = numpy.std(r['AMP'])
        vardata.append(vfinal)
        plotdata.append(rfinal)
        plotime.append(int(time))
    # plt.gca().set_facecolor('whitesmoke')
    arg = np.argsort(plotime)
    plt.scatter([plotime[i] for i in arg], [plotdata[i] for i in arg])
    plt.plot([plotime[i] for i in arg], [plotdata[i] for i in arg], color='navy')
    # plt.errorbar([plotime[i] for i in arg], [plotdata[i] for i in arg],yerr = vardata, color='navy')
    plt.xlabel('Training time (min)', size=16)
    plt.ylabel('Prediction accuracy (%)', size=16)
    plt.savefig(f'{store}No_ERROR_improvement_{name}.png', dpi=800)
    plt.close()
    
    plotdata = []
    vardata  = []
    plotime  = []
    data = sorted(error_data)
    for time, r in data:
        rfinal = (numpy.mean(r['AMP'])+numpy.mean(r['ADP']))/2.
        vardata.append(vfinal)
        plotdata.append(rfinal)
        plotime.append(int(time))
    # plt.gca().set_facecolor('whitesmoke')
    arg = np.argsort(plotime)
    plt.scatter([plotime[i] for i in arg], [plotdata[i] for i in arg])
    plt.plot([plotime[i] for i in arg], [plotdata[i] for i in arg], color='navy')
    # plt.errorbar([plotime[i] for i in arg], [plotdata[i] for i in arg],yerr = vardata, color='navy')
    plt.xlabel('Training time (min)', size=16)
    plt.ylabel('Prediction accuracy (%)', size=16)
    plt.savefig(f'{store}No_ALL_ERROR_improvement_{name}.png', dpi=800)
    plt.close()
    

process_and_plot_data('m1',[420],TEST)
try:
    process_and_plot_data('m2',[301,4201],TEST)
except:
    pass
try:
    process_and_plot_data('m3',[30,420],TEST)
except:
    pass
try:
    process_and_plot_data('m4',[301,4201],TEST)
except:
    pass

def process_and_plot_data_r2(name,timelist,TEST,r2 = False):
    path              = 'C:\\Users\\huckg\\OneDrive\\Desktop\\PPP Angewandte In Silico Final Paper\\'
    store             = 'C:\\Users\\huckg\\OneDrive\\Desktop\\\\Angewandte Paper Figures\\'
    folder            = os.listdir(path)
    dataset           = []
    R_prediction_data = {}
    r_squared         = []
    error_data       = []
    
    import numpy as np
    import seaborn as sns
    
    for i in reversed(folder):

            experiment = i.split(' ')[0]
            time = i.split(' ')[1]
            if name == experiment:
                if int(time) in timelist:
                    plot_optimal_r2(path + i + '\\', {0: TEST}, name=i, mode='TEST',r2 = r2)
  