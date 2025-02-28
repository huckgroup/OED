# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 11:49:58 2023

@author: huckg
"""
import copy 
import math

class Observable:
    def __init__(self,species,time,concentration,unit = 'min'):
        """name of species"""
        self.name = species
        
        """add the time"""
        self.time = [int(i) for i in time]
        self.unit = unit
        
        """the concentrations to match the times"""
        self.data = concentration

class TimeDependentInputs:
    def __init__(self,species,time,concentration,unit = 'min'):
        """Name of species"""
        self.name  = species
        
        """Add the time"""
        self.time  = time['Time']
        self.unit  = unit

        """The stock concentrations"""
        self.stock = species
        self.stock_concentration = {}
        
        """The concentrations to match the times"""
        self.data = concentration

def translate_time_dependent_inputs(time_dependent_inputs):
    """Time dependent inputs"""
    TDI = {}
    
    """Add the concentrations"""
    for control, obj in time_dependent_inputs.items():
        """The statevector with command inputs"""
        vector = {}
        
        """The statevector and time (start, end)"""
        for i in range(len(obj.time)-1):
            time      = (obj.time[i],obj.time[i+1])
            if time not in TDI:
                TDI[time] = copy.deepcopy(vector)
            TDI[time][control] = obj.data[i]
      
    """Stock species within each syringe"""
    syringe_load = {}
    for control, obj in time_dependent_inputs.items():
        syringe_load[control] = obj.stock
     
    """Stock concentration solution"""
    stock = {}
    for control, obj in time_dependent_inputs.items():
        for species in obj.stock:
            stock[species] = obj.stock_concentration[species]
    return TDI,syringe_load,stock

def datetime_list_to_time_vector(datetime_list):
    """import datetime"""
    from datetime import time, timedelta 
    if not datetime_list:
        return []

    """the time vector of the experiment"""
    time_vector = []
    
    """the start time and datetime is converted to a list with flowts"""
    start_time  = datetime_list[0]
    
    for dt in datetime_list:
        try:
            start_time_seconds = (start_time.hour * 60 + start_time.minute) * 60 + start_time.second
            dt_seconds = (dt.hour * 60 + dt.minute) * 60 + dt.second
            time_difference = (dt_seconds - start_time_seconds)/60.
            time_vector.append(time_difference)
        except:
            pass
    return time_vector

def parse_flow_data(folder,extentions = []):
    import os
    import pandas
    
    """Get the folders in the list """
    folder_items = os.listdir(folder)
    
    """Input in the folder"""
    for i in folder_items:
        """the input folder which contains flowrates and input output map for syringe"""
        if 'channel composition' in i:
            input_map = pandas.read_excel(folder + '\\' + i)        
        if 'channel flows' in i:
            control_inputs = pandas.read_excel(folder + '\\' + i)       
        """the data in the xlsx folder """
        if 'online data' in i:
            online_data = pandas.read_excel(folder + '\\' + i)
        if 'offline data' in i:
            offline_data = pandas.read_excel(folder + '\\' + i)
        if 'stock concentration' in i:
            stock = pandas.read_excel(folder + '\\' + i)
        if 'channel volume' in i:
            volume = pandas.read_excel(folder + '\\' + i)
                  
    """stock concentrations"""
    stock_concentration = {name:concentration[0] for name,concentration in stock.to_dict().items()}  
    
    """define the parameters that we know from the datafiles"""
    known_parameters    = {name:value[0] for name,value in volume.to_dict().items()}
    for name,value in stock_concentration.items():
        known_parameters[name+"(in)"] = value

    """Get the concentrations of species"""        
    offline_concentrations    = {}
    offline_experimental_time = {}
    
    for i in offline_data.columns:
        if i == 'Time':
            if ':' in str(offline_data[i][1]):
                offline_experimental_time[i] = datetime_list_to_time_vector(list(offline_data[i]))
            else:
                offline_experimental_time[i] = list(offline_data[i])
        elif i == 'Sample':
            pass
        else:
            offline_concentrations[i] = list(offline_data[i])
   
    """Live measurements of species"""
    online_concentrations    = {}
    online_experimental_time = {}
    
    for i in online_data.columns:
        if i == 'Time':
            if ':' in str(online_data[i][1]):
                online_experimental_time[i] = datetime_list_to_time_vector(list(online_data[i]))
            else:
                online_experimental_time[i] = list(online_data[i])
        elif  i == 'Sample':
            pass
        else:
            online_concentrations[i]    = list(online_data[i])
 
    
    """Get the inputs given to the reactor"""
    input_times = {}
    input_flows = {}
    
    for i in control_inputs.columns:
        if i == 'Time':
            input_times[i] = list(control_inputs[i]) 
        else:
            input_flows[i] = list(control_inputs[i])
   
    """online offline measurements"""
    online_time  = online_experimental_time['Time']
    offline_time = offline_experimental_time['Time']
    
    """Get and sort the observables"""
    observables,time_dependent_inputs = [],{}
    for name,concentrations in offline_concentrations.items():
        observables.append(Observable(name,offline_time,concentrations))
    for name,concentrations in online_concentrations.items():
        observables.append(Observable(name,online_time,concentrations))

        
    """Get and sort the time dependent inputs"""
    for name,concentrations in input_flows.items():
        time_dependent_inputs[name] = TimeDependentInputs(name,input_times,concentrations)    

    """channels in the dataset"""
    channels = {} 
    
    """Get the input map"""
    for i in input_map.columns:
        if i not in channels:
            channels[i] = []
        for n in input_map[i][0].split(','):
            channels[i].append(n.replace(' ',''))
   
    for c,stock in channels.items():
        time_dependent_inputs[c].stock = stock 
        for species in stock:
            time_dependent_inputs[c].stock_concentration.update({species:stock_concentration[species]})

    """get the observables and the time dependent inputs"""
    return observables,time_dependent_inputs,known_parameters