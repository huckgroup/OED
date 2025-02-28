# -*- coding: utf-8 -*-
"""
Created on Mon May 06 14:11:45 2019

@author: huckg
"""
import numpy
import scipy.integrate  
from pyDOE import *
import matplotlib.pyplot as plt 
import math
import random
import seaborn as sns
import os

from collections import defaultdict
from DataTransform import *
from Scoring import *

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
        
#this is a simple class which stores the information you are trying to extract from the model i.e. the parameters you are simulating and the scores you give 
#the simulation. This means that 
def simulate_measurement(model,
                           name = '',
                           #conditions
                           conditions = {},
                           #parameters
                           parameters = {},
                           #lumped
                           lumped = {},
                           #see if there are time dependent parameters
                           time_dependent_parameters = {},
                           #coordinates instead of TDI
                           coordinates = (),
                           #fitfunction of the measurement
                           fitfunc = 'leastsquaresfit',
                           #flowrate
                           flowrate = 0, 
                           time = (0,200),
                           dt = 1,
                           
                           store = '',
                           #show the result
                           show = False):
    
    from Model import ModelVariables
    from SolveSystem import ModelSolver
    
    start,end = time
    """create the parameter vectors and initial conditions to solve the system"""
    variables = ModelVariables(model,conditions = conditions,modification = parameters)
    if len(coordinates) != 0:
        """extract the coordinates"""
        coordinates,space = coordinates
        """set the time dependent parameters"""
        variables.set_time_dependent_variables(coordinates,space)
        """extract the manual TDI"""
        time_dependent_parameters = variables.amici_TDI
        """"solve the system"""

    solution = ModelSolver(model,variables = variables,simtime = time,dt = dt,manual_TDI = time_dependent_parameters)
    """add the sensitivities to the vector"""
    data = solution.__getData__()            
    """"data of the sim"""
    experimental_data = [(data.rdata,data.time)]
    """store measurements"""  
    simdata = []
    """get the observables"""
    if model.observables == []:
        model.observables = model.states
    for species,sd in data.simdata.items():
        if species in model.observables:
            simdata.append(Observable(species,data.time,sd))
    measurement       = MeasurementObject(simdata,
                                          time = data.time,
                                          name = name,
                                          conditions = conditions,
                                          lumping = lumped,
                                          time_unit = 'min',
                                          time_dependent_parameters=time_dependent_parameters, 
                                          store = store,
                                          show = show)
    measurement.model = model
    return measurement


""""this is the M template used in the pipeline to fit data and or calculate likelihoods 
for (after some minor adjustments are made to your version)"""
class MeasurementObject:
    def __init__(self,
                 #name,rawdata and time
                 rawdata,
                 
                 #supply a time vector to the dataste
                 time = False,
                 
                 #name of the measurement
                 name = '',
                
                 #Initial conditions
                 initial_conditions        = [],               
                
                 # there can be conditions and initial conditions
                 conditions                = {},
                 
                 #Time dependent parameters i.e. if this exist solvers are broken appart
                 time_dependent_parameters = {},    
                 
                 #Define which observables are lumped into a single observalbe
                 lumping                   = {},
                 
                 #time units within the data
                 desired_time    = 'min',
                 time_unit       = 'min',

                 interpolation   = False,
                 store           = False,
                 show            = False):
        
        
        self.name = name
        """raw data as presented in the excell sheet"""
        self.rawdata           = rawdata
        """time over which this was measured"""
        self.measurement_time  = time
        """time dependent parameters i.e. inputs as experiment goes"""
        self.time_dependent_parameters = time_dependent_parameters   

        """should there be an interpolation of the data"""
        self.interpolation = interpolation
        
        """set the timescales"""
        self.desired_time    = desired_time
        self.time_unit       = time_unit

        """controlparameters"""
        if self.time_dependent_parameters != {}:
            """plug the measurements and models in to the same folder
            and run asssign the right model to the right measurement"""
            dct = list(self.time_dependent_parameters.values())[0]
            self.control = list(sorted(list(dct.keys())))
            
        """deal with the open and or closed system we are working in flow reactors or
        under batch conditions compounds that flow into the reactor as initiators are denoted 
        in the measurement file, however if there is no flow this becomes an initial_condition"""
        self.conditions         = conditions
        """ïf there is no flow i.e. just batch the initial conditions become the conditions"""
        self.initial_conditions = initial_conditions
        if not self.initial_conditions:
            self.initial_conditions = self.conditions
            
        """set the inputs vector"""
        self.inputs = self.conditions.keys()

        """interpolate data vector to perform least squares fitting"""
        self.profile = {}

        if not self.interpolation:
            pass
        else:
            if type(self.rawdata) == dict:
                self.interpolation = True
            else:
                self.interpolation = False  
        if self.interpolation:
            """extract the information on the time vector """
            self.t_start = self.measurement_time[0]
            self.t_end   = self.measurement_time[-1]

            """get the datavector and interpolate the dataframe, if interpolation == False then we do NOT!"""
            for species,data in self.rawdata.items():      
                self.profile[species],self.time,self.conversion = data,self.measurement_time,1
                if len(data) != len(range(int(self.t_end))):
                    self.profile[species],self.time,self.conversion = interpolate_dataframe(self.measurement_time,data,self.time_unit,desired = desired_time)  
           
            """timevector"""                    
            self.time     = list(self.time)
            
            """get the timestep by taking time vector and """
            self.dt       = self.time[1] - self.time[0]
            self.t_start *= self.conversion
            self.t_end   *= self.conversion
        else:
            for obj in self.rawdata:
                self.profile[obj.name] = [(obj.time[i],obj.data[i]) for i in range(len(obj.time))]
            self.time = 0
            for name,profile in self.profile.items():
                t,d = profile[-1]
                if t > self.time:
                    self.time = copy.deepcopy(t)
            self.time = list(range(self.time))

            
        """get the observables"""
        self.observables = list(self.profile.keys())
        
        """see which observables in the profile are lumped this needs to be automated in
        the script from which you either load your model or manually set"""
        self.lumped      = lumping

        """the path of the desktop folder"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')    
        """directory where the files are stored"""
        self.directory = desktop + "\\__SimulatedExperiments__\\"
    
    def prune(self,prune):
        """this removes part of the dataset with obvious calibration errors"""
        start,end,species = prune
        """the integer time values deleted"""
        time = range(start,end,1)
        """update profile"""
        self.profile[species] = [(t,d) for t,d in self.profile[species] if t not in time] 
        
    def show(self,directory = '',name = '',scale = '',store = True,show = False,temperature = False,flowrate = False,observables = []):
        sns.set()
        """"overwrite the dictionaries if needed"""
        if directory != '':
            self.directory = directory
        if name != '':
            self.name = name
            
        """plot data of experiment"""
        fig = plt.figure(figsize=(10,5))
        ax = plt.subplot(1,2,1)
        if observables == []:
            observables = [i for i in self.profile.keys()]
        
        print(self.observables)
        for k,v in self.profile.items():
            if k in observables:
                if self.interpolation:
                    ax.plot(v,label = k)
                    ax.scatter(self.measurement_time*self.conversion,self.rawdata[k])
                else:
                    time,concentration = zip(*self.profile[k])
                    ax.scatter(time,concentration,s = 5,marker = '^',c = 'k')
                    ax.plot(time,concentration,label = k)
        ax.set_title("Measurement Data",size = 14)
        ax.set_xlabel("Time ({})".format(self.desired_time),size = 14)
        ax.set_ylabel("Concentration",size = 14)
        # plt.yscale('log')
        ax.legend(fancybox = True)

        if len(self.profile) < 10:
            ax.legend(fancybox = True)

        ax = plt.subplot(1,2,2)
        if len(self.time_dependent_parameters) == 0:
            print("no tdi",self.time_dependent_parameters)
            """subplot with the conditions of the experiment"""
            if temperature and flowrate:
                    textstr = '\n'.join(('Temperature: {}'.format(temperature,),
                                         'Flow rate: {}'.format(flowrate,),))
                    props = dict(boxstyle='round', facecolor='grey')
                    # place a text box in upper left in axes coords
                    ax.text(0.08, 0.97, textstr, transform=ax.transAxes, 
                            fontsize=14,verticalalignment='top', bbox=props)     
            """plot bar plot to show condition with which experiment was performed"""
            xrn    = range(len(self.conditions))
            labels = []
            values = []
            for k,v in self.conditions.items():
                if k != "O":
                    labels.append(k)
                    values.append(v) 
            """logarithmic values of the conditions that were used"""
            logvalues = []
            for i in values:
                if type(i) == list or type(i) == tuple:
                    print ("""A list or tuple has been past as a number
                    to the array There is probably a delimiter issue at the parsing section of the data
                    adjust it to ; or , and check if it solves the issue""")
                if i != 0:
                    logvalues.append(math.log10(float(i)))
                else:
                    logvalues.append(0)
                    
            xrn = numpy.argsort(numpy.argsort(logvalues))
            """optional scale to plot the different conditions"""
            if scale == "log":
                ax.bar(xrn, logvalues, align='center',edgecolor = 'k',alpha = 0.7)
            else:
                ax.bar(xrn, values, align='center',edgecolor = 'k',alpha = 0.7)
            """set the xticks of the experiment"""
            ax.set_xticks(xrn)
            ax.set_xticklabels(labels,rotation=90)
            ax.set_title("Control Parameters")
            plt.legend(fancybox = True)
            
            if store:
                plt.savefig(self.directory + self.name + "_EXP.png",dpi = 600)
                plt.savefig(self.directory + self.name + "_EXP.svg",dpi = 600)
            if show:
                plt.show()
            # plt.close()
            
        else:
            time,control = [],[]
            t = sorted(self.time_dependent_parameters.keys())
            for i in t:
                start,end = i
                time.append(start)
                parameters = self.time_dependent_parameters[i]
                control.append(parameters)
            control = LDtoDL(control)
            for parameter,value in control.items():
                plt.plot(time,value,label = parameter)
            if len(parameters) < 20:
                plt.legend(fancybox = True)
            
            plt.xlabel("Time")
            plt.ylabel("uL")
            plt.tight_layout()
            
            if store:
                plt.savefig(self.directory + self.name + "_EXP.png",dpi = 600)
                plt.savefig(self.directory + self.name + "_EXP.svg",dpi = 600)
            if show:
                plt.show()
                


class GenerateExperiments:
    def __init__(self,model,
                 #number of control sets sampled
                 samples = 100,
                 experiments = 5,
                 #simulation time
                 simulation_time = (0,500),
                 dt = 1,
                 #superset created i.e. 100 supersets of 5 experimental conditions
                 superset = 100):
        
        """the optimization and and information package"""
        from OptimizeExperiment import OptimizeInformation
        """the model"""
        self.model = model
        """set measurement dicts"""
        self.measurements = {}
        """"generate a set of experiments"""
        conditions = build_global_search(model.p_space,include = model.control,samples = samples, order = int(math.log10(len(list(model.p_space.values())[-1]))))

        for i in range(len(conditions)):
            self.measurements[i] = simulate_measurement(model,parameters = model.fixed,conditions = conditions[i],time = simulation_time,fitfunc = "leastsquaresfit",show = False)
        """simulation data of the system"""
        data = {i:[] for i in model.observables}
        """lsq values from the mean devation"""
        lsq  = {i:[] for i in model.observables}
        
        if len(conditions) < 10:
            self.measurement_superset = self.measurements
        else:
            """"return both the most divergent and least divergent sets of the mean"""
            for o in model.observables:
                """ the mean of the data for a given observable"""
                mean =  numpy.zeros((len(self.measurements[0].profile[o])))
                for i in self.measurements.values():  
                    t,y = zip(*i.profile[o])
                    mean += y
                mean   /= len(self.measurements)
                for i in range(len(self.measurements)):
                    t,y = zip(*self.measurements[i].profile[o])
                    data[o].append(y)
                    lsq[o].append(sum(y-mean))
    
            """metric per observable"""
            self.distance_set = {i:[] for i in model.observables}
            """sort the data according to their lsq distance from the mean"""
            for observable,error in lsq.items():    
                ranks = list(numpy.argsort(error))
                """add min max from mean"""
                self.distance_set[observable].append(ranks.index(min(ranks)))
                self.distance_set[observable].append(ranks.index(max(ranks)))
                """"distance from the mean"""
                for i in range(0,len(error),int(len(error)/experiments-2)):
                    """find index of a this ranked simulation"""
                    idx = ranks.index(i)
                    """add to set"""
                    self.distance_set[observable].append(idx)
    
            """append the conditions to be used"""
            self.conditions = []
            for o,indices in self.distance_set.items():
                for i in indices:
                    self.conditions.append(conditions[i])
                           
            """this selection is rank based however these are equistant based on index
            ergo it does not take into acccount the underlying distribution of ranks into account
            therefore we attempt to flatten the distribution by inverting the probability of the 
            densities and subserquently randomly sample conditional sets (to be used for
            collinearity and senstivity experiments), therefore we cut it into pieces (number of conditions) and
            select from windows"""                                            
            self.superset = []
            for n in range(superset):
                for observable,error in lsq.items():
                    lsq_selection = []
                    """get window size"""
                    factor = int(len(error)/experiments)
                    """get indices"""
                    window_size = [0 + factor for i in range(0,len(error),factor)]
                    """loop through windows"""
                    index = 0
                    for window in range(len(window_size)):
                        arr_slice = numpy.array(sorted(error))[index:window_size[window]]
                        """window index"""
                        index += window
                        """select value in the window and add to list"""
                        lsq_selection.append(random.choice(arr_slice))
                    """add the list to the superset"""
                    self.superset.append([error.index(i) for i in lsq_selection])
                    
            """set of measurements"""
            self.measurement_superset = []
            for i in self.superset:
                self.measurement_superset.append([self.measurements[j] for j in i])
            self.conditions = []
            for i in self.measurement_superset:
                self.conditions.append([j.conditions for j in i])
    
        """the path of the desktop folder"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')    
        """directory where the files are stored"""
        self.directory = desktop + "\\__SimulatedExperiments__\\"
        if not os.path.exists(self.directory):            
            os.makedirs(self.directory)    
        filenames = []
          
            
    def return_single_measurement(self,scale = '',store = False,desired_time = "min",name = ''):
        """return the single measurement"""
        measurement = random.choice(random.choice(self.measurement_superset))
        if store == True:
            sns.set()
            """plot data"""
            fig = plt.figure(figsize=(10,5))
            ax = plt.subplot(1,2,1)
            for k,v in measurement.profile.items():
                ax.plot(v,label = k)
                ax.scatter(measurement.measurement_time,measurement.rawdata[k])
    
            ax.set_title("Measurement Data",size = 14)
            ax.set_xlabel("Time ({})".format(desired_time),size = 14)
            ax.set_ylabel("Concentration (M)",size = 14)
            ax.legend(fancybox = True)
    
            ax = plt.subplot(1,2,2)
            """subplot with the conditions of the experiment"""
            if measurement.temperature and measurement.flowrate:
                    textstr = '\n'.join(('Temperature: {}'.format(measurement.temperature,),
                                         'Flow rate: {}'.format(measurement.flowrate,),))
                    props = dict(boxstyle='round', facecolor='grey')
                    # place a text box in upper left in axes coords
                    ax.text(0.08, 0.97, textstr, transform=ax.transAxes, 
                            fontsize=14,verticalalignment='top', bbox=props) 
                    
            """plot bar plot to show condition with which experiment was performed"""
            labels,values = [],[]
            
            """the measurement conditions"""
            for k,v in measurement.conditions.items():
                labels.append(k)
                values.append(v) 
            xrn = range(len(values))
            
            """logarithmic values of the conditions that were used"""
            logvalues = []
            for i in values:
                if type(i) == list or type(i) == tuple:
                    print ("""A list or tuple has been past as a number
                    to the array There is probably a delimiter issue at the parsing section of the data
                    adjust it to ; or , and check if it solves the issue""")
                if i != 0:
                    logvalues.append(math.log10(float(i)))
                else:
                    logvalues.append(0)
                        
                """optional scale to plot the different conditions"""
                if scale == "log":
                    ax.bar(xrn, logvalues, align='center',edgecolor = 'k',alpha = 0.7)
                else:
                    ax.bar(xrn, values, align='center',edgecolor = 'k',alpha = 0.7)
                """set the xticks of the experiment"""
                ax.set_xticks(xrn)
                ax.set_xticklabels(labels,rotation=90)
                ax.set_title("Control Parameters")
                plt.tight_layout()
                plt.savefig(self.directory + "SM " + name + ".png",dpi = 600) 
        return {0:measurement} 
    
    def return_random_measurents(self,number = 3,store = False,desired_time = "min",name = ''):
        """measurements"""
        measurements = {}
        """potential"""
        potential = list(range(len(self.measurements)))
        for i in range(number):
            c = random.choice(potential)
            """"delete potential"""
            potential.remove(c)
            """update measurement list"""
            measurements[i] = self.measurements[c]

        if store == True:           
            sns.set()
            """get identical colorset"""
            colorset = [numpy.random.rand(3,) for i in range(len(measurements[0].profile))]
            """plot data"""
            fig = plt.figure(figsize=(12,5))
            ax = plt.subplot(1,2,1)
            """plot the measurement raw data"""
            for m in range(len(measurements)):
                cnt = 0
                for k,v in measurements[m].profile.items():
                    if m == 0:
                        ax.plot(v,label = k,c = colorset[cnt])
                    else:
                       ax.plot(v,c = colorset[cnt])                       
                    cnt += 1
                
            """set data labels and measurements"""
            ax.set_title("Measurement Data",size = 14)
            ax.set_xlabel("Time ({})".format(desired_time),size = 14)
            ax.set_ylabel("Concentration (nM)",size = 14)
            ax.legend(fancybox = True)
    
            """second subplot with the conditions"""
            ax = plt.subplot(1,2,2)
            conditions = [measurements[m].conditions for m in range(len(measurements))]
            conditions = LDtoDL(conditions)
            
            """get the labels and the values"""
            labels = conditions.keys()
            values = [conditions[i] for i in labels]
            xrn = range(len(labels))
            
            """plot the barplot of the conditions"""
            for i in range(len(values)):
                for value in values[i]:
                    ax.bar(xrn[i], value, align='center',edgecolor = 'k',alpha = 0.2)
                    
            """set the xticks of the experiment"""
            ax.set_xticks(xrn)
            ax.set_xticklabels(labels,rotation=90)
            ax.set_title("Control Parameters")
            plt.tight_layout()
            plt.savefig(self.directory + "MRS " + name + ".png",dpi = 600) 
            plt.close()
        return measurements
    
    def return_multiplexed_measurements(self,store = False,desired_time = "min",name = ''):
        """return a set of measurements instead of the single measurement"""
        measurements = random.choice(self.measurement_superset)
        if store == True:           
            sns.set()
            """get identical colorset"""
            colorset = [numpy.random.rand(3,) for i in range(len(measurements[0].profile))]
            """plot data"""
            fig = plt.figure(figsize=(12,5))
            ax = plt.subplot(1,2,1)
            """plot the measurement raw data"""

            for m in range(len(measurements)):
                cnt = 0
                for k,v in measurements[m].profile.items():
                    if m == 0:
                        ax.plot(v,label = k,c = colorset[cnt])
                    else:
                       ax.plot(v,c = colorset[cnt])   
                    cnt += 1
                
            """set data labels and measurements"""
            ax.set_title("Measurement Data",size = 14)
            ax.set_xlabel("Time ({})".format(desired_time),size = 14)
            ax.set_ylabel("Concentration (M)",size = 14)
            ax.legend(fancybox = True)
    
            """second subplot with the conditions"""
            ax = plt.subplot(1,2,2)
            conditions = [measurements[m].conditions for m in range(len(measurements))]
            conditions = LDtoDL(conditions)
            
            """get the labels and the values"""
            labels = conditions.keys()
            values = [conditions[i] for i in labels]
            xrn = range(len(labels))
            
            """plot the barplot of the conditions"""
            for i in range(len(values)):
                for value in values[i]:
                    ax.bar(xrn[i], value, align='center',edgecolor = 'k',alpha = 0.2)
            """set the xticks of the experiment"""
            ax.set_xticks(xrn)
            ax.set_xticklabels(labels,rotation=90)
            ax.set_title("Control Parameters")
            plt.tight_layout()
            plt.savefig(self.directory  + name +"Multiplex.png",dpi = 600)
            plt.close()
        return {i:measurements[i] for i in range(len(measurements))}
    
    def return_oscillating_measurement(self,features = {},time = (0,200),dt = 2,store = False,name = ''):
        from scipy import signal
        from UliEngineering.SignalProcessing.Simulation import sine_wave
        from UliEngineering.SignalProcessing.Simulation import square_wave
        from UliEngineering.SignalProcessing.Simulation import sawtooth
        
        """set times"""
        start,end = time
        """"the arrange generates a timespan from start -> end -dt"""
        end += dt
        """toi vector"""
        ctimes = range(start,end,dt)
     
        """length of the features"""
        if len(features) == 0:
            features = {i:'' for i in self.model.control}
            for i in self.model.control:
                """"set amplitude at halfpoint parameters"""
                amplitude = self.model.p_space[i][int(self.model.spacing/2)]
                """"random period"""
                period    = random.choice(range(15,30))
                """oscillate around base"""
                base      = self.model.p_space[i][int(self.model.spacing/2)]
                """update features"""
                features[i] = copy.copy((period,amplitude,base))
                
        """"set the pattern you want and generate it"""
        swt,osc,bnr = {},{},{}
        for parameter,variables in features.items():
            """ampltidudes are hardset"""
            period,amplitude,base = variables
            
            """binary"""
            bnr[parameter] = {}                    
            x = numpy.arange(start,end,1)
            y = square_wave(frequency=(end/period), samplerate = end, amplitude=amplitude, offset=base) 
            y += base
            """update binary pattern"""
            for i in range(len(ctimes)-1):
                bnr[parameter][(ctimes[i],ctimes[i+1])] = y[ctimes[i]]   
            """sawtooth"""
            swt[parameter] = {}
            x = numpy.arange(start,end,1)
            y = sawtooth(frequency=(end/period), samplerate = end, amplitude=amplitude, offset=base)
            y += base
            """update sawtooth pattern"""
            for i in range(len(ctimes)-1):
                swt[parameter][(ctimes[i],ctimes[i+1])] = y[ctimes[i]]   
            """oscillations"""         
            osc[parameter] = {} 
            x = numpy.arange(start,end,1)
            y = sine_wave(frequency=(end/period), samplerate = end, amplitude=amplitude, offset=base)
            y += base
            """update oscillations"""
            for i in range(len(ctimes)-1):
                osc[parameter][(ctimes[i],ctimes[i+1])] = y[ctimes[i]]   
        """reorient values"""
        bn,sw,os = {},{},{}
        """count"""
        cnt = 0
        for i in [bnr,swt,osc]:
            for time in list(i.values())[-1]:
                """"create the amici and python tdi input vectors"""
                dct = {}
                """loop throug paramteers to create dicts"""
                for parameter in i.keys():
                    dct[parameter] = i[parameter][time]
                if cnt == 0:
                    bn[time] = dct
                if cnt == 1:
                    sw[time] = dct
                if cnt == 2:
                    os[time] = dct
            cnt += 1
        """simulate this stuff and create a measurement object for each pattern
        by doing so we can start fitting against this input"""
        binary = simulate_measurement(self.model,dt = 1,name = name + "Binary",time_dependent_parameters = bn,store = store)  
        sawtooth = simulate_measurement(self.model,dt = 1,name = name + "Sawtooth",time_dependent_parameters = sw,store = store)  
        oscillations = simulate_measurement(self.model,dt = 1,name = name +  "Oscillation",time_dependent_parameters = os,store = store)    
        return {0:oscillations}
    
    def return_random_pulse_measurement(self,
                            #index and mutation size and mutation number
                             pindex          = 30,
                             mutation_size   = 50,
                             mutation_number = 100,
                             #lenth of pulses
                             plength = 15,
                             #pulse start
                             pulsestart = 3,
                             #Forced mutation
                             forced_mutation = {},
                             conditions      = {},
                             #time
                             time = (0,48),
                             name = '',
                             store = False,
                             desired_time = 'Hour'):
        
        coordinates,space = {},{}
        """"here we will generate a number of 
        spaces where whe can sample and create
        a input space in a MIMO system"""   
        start,end = time
        for i in self.model.control:
            controlbounds = {i:self.model.boundaries[i]}
            """define boundaries of the pulse space"""
            cs,exclude  = define_parameter_boundaries(controlbounds,spacing = pindex) 
            """pulse vector and plength"""
            coordinates[i] = {(j,j+plength):int(pindex/2) for j in range(pulsestart,end,plength)}
            """space of the system"""
            space.update(cs)
        """pulse space transferred to the model"""
        self.model.pulsespace = space
        """pulse sequnece or series"""
        times = list(sorted(list(coordinates.values())[-1].keys()))
        """mutatory list"""
        mutationlist = []
        """loop through the mutations"""
        for i in range(mutation_number):
            """randomly chosen size"""
            size = random.choice(range(mutation_size))
            """parameter selected for mutation"""
            parameter = random.choice(self.model.control)
            """time of the pulse"""
            pulse = random.choice(range(len(times)))
            """pulse go left or right"""
            chosensequence = []
            if random.random() <= 0.5:
                """go left along the number line [(),(),()...<- (chosen),(),(),()]"""
                start_index = pulse-size
                if start_index < 0:
                    start_index = 0
                for t in range(start_index,pulse,1):
                    try:
                        chosensequence.append(times[t])
                    except IndexError:
                        break
            else:
                """go right along the number line [(),(),(),(chosen) ->...,(),(),()]""" 
                end_index = pulse+size
                if end_index > len(times)-1:
                    end_index = len(times)
                for t in range(pulse,end_index,1):
                    try:
                        chosensequence.append(times[t])
                    except IndexError:
                        break
            """size of the index space"""
            size = random.choice(range(pindex))
            """append mutation"""    
            if len(chosensequence) != 0:
                mutationlist.append((parameter,chosensequence,size))

        """tracking mutations"""
        for parameter,sequence,index in mutationlist:
            for i in sequence:
                """update the coordinate vector"""
                coordinates[parameter][i] = index 
        """enter coordinates and simulates"""
        for k,v in forced_mutation.items():
            pass
        
        measurement = simulate_measurement(self.model,name = name + "RandomPulse",dt = 1,coordinates = (coordinates,space),conditions= conditions,store = store)
        
        if store == True:
            sns.set()
            """plot data"""
            fig = plt.figure(figsize=(10,5))
            ax = plt.subplot(1,2,1)
            for k,v in measurement.profile.items():
                ax.plot(v,label = k)
                ax.scatter(measurement.measurement_time,measurement.rawdata[k])
    
            ax.set_title("Measurement Data",size = 14)
            ax.set_xlabel("Time ({})".format(desired_time),size = 14)
            ax.set_ylabel("Concentration (M)",size = 14)
            ax.legend(fancybox = True)
    
            ax = plt.subplot(1,2,2)
            """subplot with the conditions of the experiment"""
            if measurement.temperature and measurement.flowrate:
                    textstr = '\n'.join(('Temperature: {}'.format(measurement.temperature,),
                                         'Flow rate: {}'.format(measurement.flowrate,),))
                    props = dict(boxstyle='round', facecolor='grey')
                    # place a text box in upper left in axes coords
                    ax.text(0.08, 0.97, textstr, transform=ax.transAxes, 
                            fontsize=14,verticalalignment='top', bbox=props) 
                    
            """plot bar plot to show condition with which experiment was performed"""
            labels,values = [],[]
            
            """the measurement conditions"""
            for k,v in measurement.conditions.items():
                labels.append(k)
                values.append(v) 
            xrn = range(len(values))
            
            """logarithmic values of the conditions that were used"""
            logvalues = []
            for i in values:
                if type(i) == list or type(i) == tuple:
                    print ("""A list or tuple has been past as a number
                    to the array There is probably a delimiter issue at the parsing section of the data
                    adjust it to ; or , and check if it solves the issue""")
                if i != 0:
                    logvalues.append(math.log10(float(i)))
                else:
                    logvalues.append(0)
                        
                """optional scale to plot the different conditions"""
                if scale == "log":
                    ax.bar(xrn, logvalues, align='center',edgecolor = 'k',alpha = 0.7)
                else:
                    ax.bar(xrn, values, align='center',edgecolor = 'k',alpha = 0.7)
                """set the xticks of the experiment"""
                ax.set_xticks(xrn)
                ax.set_xticklabels(labels,rotation=90)
                ax.set_title("Control Parameters")
                plt.tight_layout()
                plt.savefig(self.directory + "SM " + name + ".png",dpi = 600)  
                plt.close()
                
        return measurement
            
    def optimized_pulse_measurement(self,#which parameters are pulsed
                         time_dependent_parameters  = [],
                         initial_control_conditions = {},
                         conditions = {},
                         #what is measuremd
                         observables = [],
                         #number of generations and agents per generation
                         generations = 100,
                         agents = 10,
                         #length of pulse time of experiment, start of sequence number of indices
                         time = (0,200),
                         dt = 1,
                         pulsestart = 2,
                         plength = 1,
                         #size of control region
                         pindex = 10,
                         #number of start samples
                         multistart = 10,
                         #scoring functions
                         sf = "D_Fisher",
                         name = '',
                         store = False,
                         gif = False):
        
        from OptimizeExperiment import OptimizeInformation
        """optimize the measurement pulse sequence to extract the most information possible"""
        optimization = OptimizeInformation(self.model,
                                          time_dependent_parameters = time_dependent_parameters,
                                          initial_control_conditions = initial_control_conditions,
                                          conditions = conditions,
                                          observables = observables,
                                          generations = generations,
                                          agents = agents,
                                          time = time,
                                          dt = dt,
                                          pulsestart = pulsestart,
                                          plength = plength,
                                          pindex = pindex,
                                          multistart = multistart,
                                          sf = sf)
        
        """create gif of the optimization"""
        if gif:
            optimization.gif()           
        """fittest in the set, i.e. the manually controlled input set"""
        coordinates = optimization.fittest.coordinate_track[optimization.fittest.fittest]
        """create a measurement object with this vector"""
        measurement = simulate_measurement(self.model,name = name + "OptimizedPulse",dt = 1,coordinates = (coordinates,optimization.fittest.pulsespace),store = store)
        return measurement
    
    
#