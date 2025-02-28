 # -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:49:05 2020

@author: Bob
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 06 10:43:21 2019

@author: huckg
"""
from OptimizationOperations import *
from Scoring import *
from Operations import *
from Distributions import *
from DataTransform import *
from PlotData import *

import pandas
import pandas as pd
import seaborn as sns
import time as timesystem


from operator import itemgetter
from itertools import chain

import matplotlib.pylab as plt
import statsmodels.api as sm
import os
import scipy
import random
import numpy
import math
import copy

class Convergence:
    def __init__(self,):
        """the highest score"""
        self.lowest_SD  = 1*10**10
        """the best simdata"""
        self.simdata     = []
        """iteration where the best simdata happened"""
        self.iteration   = []
        """best position"""
        self.position    = []
        """best coordinates"""
        self.coordinates = []
        """scores"""
        self.scores      = []
        """converge"""  
        self.converge    = []
        
        """update count, tracks how often the data has been updated"""
        self.update_count = 1
        
    def update_convergence(self,agents,it,threshold = -3):
        """update the counter"""
        self.update_count += 1
        """agents and the agent number ID"""
        nmbr,agents = zip(*agents.items())
        """minimal values"""
        SD_values = [i.SD for i in agents]
        """smallest value"""
        fittest_agent = min(SD_values)
        """get index of smallest"""
        index = SD_values.index(fittest_agent)
        if min(SD_values) < self.lowest_SD:
            self.lowest_SD = copy.copy(fittest_agent)
            """print the score and iteration it was updated"""
            print(str(round(self.lowest_SD,2)) + '%','in Iteration ' + str(it))
            
            """update simdata, coordinates, position, iteration etc."""
            self.iteration.append(copy.copy(it))
            self.scores.append(copy.copy(fittest_agent))
            self.position.append(copy.copy(agents[index].fitset))
            self.coordinates.append(copy.copy(agents[index].coordinates))
            self.simdata.append(agents[index].simdata)
            self.fittest_agent = copy.copy(agents[index])
     
        """there are 3 miminum termination criteria, no improvement in 20 generations
           atleast 30 generations
           the fit needs to be atleast 0.1 devation on average if you want to stop earlier"""
        terminate = False
        indices   = nmbr
        """get list agents"""
        booleans,difference = [],[]
        """set two different termination criteria"""
        if self.update_count > 50:
            try:
                print(self.iteration[-1] - self.update_count)
            except:
                pass
            if len(self.iteration) > 1:
                """no improvement over 100 generations = no movement"""
                if abs(self.iteration[-1] - self.update_count) > 20:
                    terminate = True

            
        for agent in agents:
            """booleans of the convergence"""
            boolean = agent.accepted_movements[threshold:]
            if any(boolean) == True:
                booleans.append(False)
            else:
                booleans.append(True) 
            difference.append(float(agent.SD)/float(self.lowest_SD))
            
        """overwrite agents if they havent moved or if the difference between worst and best is >5"""
        for i in range(len(difference)):
            if difference[i] > 1.05 and it > 30 and it%15 == 0:
                print(difference,'The unfit Agents have been re-initialized')
                """superimpose sensitivities and mutation probililities:"""
                """fit = max(list(self.fittest_agent.forward_sensitivity.keys()))
                agents[indices[i]].update_mutation_probability(self.fittest_agent.forward_sensitivity[fit],it)"""
                agents[indices[i]].SD              = copy.copy(self.lowest_SD)
                agents[indices[i]].height[0].score = copy.copy(self.lowest_SD)
                """superimpose position and coordinates"""
                agents[indices[i]].coordinates  = copy.copy(self.coordinates[-1])
                agents[indices[i]].position     = copy.copy(self.position[-1])
                """superimpose position and coordinates in the tracker"""
                agents[indices[i]].coordinate_track[agents[indices[i]].fittest] = copy.copy(self.coordinates[-1])                
                agents[indices[i]].track[agents[indices[i]].fittest]            = copy.copy(self.position[-1])
                    
                    
        """modify converge"""
        self.converge.append(self.lowest_SD)
        """redict the agents"""
        agents = {nmbr[i]:agents[i] for i in range(len(agents))}
        return terminate,agents,self.converge

class Agent:
    def __init__(self,measurements,p_space,include,pos = []):
        """add measurement and model to agent"""
        self.measurements = measurements #dictionary with the state and the data (ineterpolated and the time)   
        self.modelnames   = list(set([self.measurements[i].model for i in range(len(self.measurements))]))
        self.models       = [self.measurements[i].model for i in range(len(self.measurements))]
        
        """the boundaries in the system"""
        self.boundaries,self.observables,self.control,self.fixed = {},[],[],[]
        for model in self.models:
            """unpack the information in the models"""
            for k,v in model.boundaries.items():
                self.boundaries[k] = v
            for i in model.observables:
                self.observables.append(i)
            for i in model.control:
                self.control.append(i)
        """remove duplicates"""
        self.observables = list(set(self.observables))
        self.control = list(set(self.control))
        
        self.include = include
        if len(include) == 0:
            self.include = [i for i in self.boundaries.keys() if i not in self.control] 
        for i in self.control:
            self.fixed.append(i)
        for parameter,value in self.boundaries.items():
            if i not in include and include != []:
                self.fixed.append(parameter)
        
        """the needed information for the parameters and and parameter coordinates for mutation"""
        self.p_space = p_space
        self.spacing = len(list(self.p_space.values())[-1])
        self.forward = True
        
        """ inverted dictionary for parametres"""
        self.c_space = {}
        for p,dct in self.p_space.items():
            self.c_space[p]  = {v:k for k, v in self.p_space[p].items()}
        
        """these are the continuously updated parameter sets, its corresponding coordinates
        and the score sequence. They are only updated if the score is better"""
        self.position = pos
        self.coordinates = {k:self.c_space[k][v] for k,v in self.position.items()}
        self.height = None
        
        self.pareto_parameters = []
        self.last_forward_update = 0
        """information of the EA et al. is still information on the kinetics
        of the model and range of potential behavior thus we store the data
        and the transition moreover we can track their status over time"""
        self.track    = {}
        self.coordinate_track = {}
        
        """the list of accepted mutations/recombinations"""
        self.progress = {}
        """"forward motion acceptance"""
        self.accepted_movements = []
        """track agents which keep track of the evolution itself"""
        self.forward_sensitivity  = {}
        
    def initial_height(self,scores,simdata):
        """"this is the initial height of the system as obtained by the
        monte carlo before we start with the optimizaton routine"""
        self.height  = scores  
        self.simdata = simdata
        """"these will get overwritten (probably) as the evolution starts"""
        self.fittest = 0
        """track first element in optimization"""
        self.progress[self.fittest] = (scores,copy.copy(self.position))
        """set initial position of the agent"""
        self.position = self.track[self.fittest]
        """initial average deviation"""
        self.fitset = self.position
#        self.average_SD,self.median_SD,self.largest = calculate_averageSD(self.simdata,self.measurement)
        self.SD = calculate_deviation(self.simdata,self.measurements)
        
    def update_position(self,coordinates,it):
        """new coordinate from the recombinations"""
        self.coordinates = coordinates
        for k,v in self.coordinates.items():
            self.position[k] = self.p_space[k][v]
        """update the tracker"""
        self.track[it] = copy.copy(self.position)
        self.coordinate_track[it] = copy.copy(self.coordinates)
        
    def random_mutation(self,it,attr = []):
        """Mutate the parameters based on the current state of the"""
        number,size,rnd = mutationrange(self.spacing)
        """store the mutations that are made"""
        mutations = []
        
        for i in range(number):
            step,p = random.choice(range(size)),random.choice(self.include)
            left,right = (range(self.coordinates[p]),range(self.coordinates[p],len(self.p_space[p]),1))
 
            """move left or right"""
            sign = random.choice([-1,1]) 
            
            """move to the left or right"""
            if sign < 0:
                if len(left) < step:
                    if left:
                        step = random.choice(left)
                    else:
                        step = 0
                mutations.append((p,self.coordinates[p] + (step*sign)))
            else:
                if len(right) < step:
                    if right:
                        step = random.choice(range(len(right)))
                    else:
                        step = 0
                mutations.append((p,self.coordinates[p] + (step*sign)))
            if rnd:
                mutations.append((p,random.choice(range(self.spacing)))) 
                
        """overwrite the coordinates and position by moving"""    
        for p,c in mutations:
            if c == self.spacing:
                c -= 1
            self.coordinates[p] = c
            self.position[p] = self.p_space[p][c]
            
        """update the tracker"""
        self.track[it] = copy.copy(self.position)
        self.coordinate_track[it] = self.coordinates
        
    def micromutation(self,it,mutation_number = 2):
        """5% change"""
        size = int(self.spacing*0.03)
        """Mutate the parameters based on the current state of the"""
        number = random.choice(range(1,mutation_number,1))
        """mutate the parameters in line with the model"""
        if random.random() > 0.33:
            number = len(self.models)
        """store the mutations that are made"""
        mutations = []
        for i in range(number):
            """get the to be mutated parameter"""
            step = random.choice(range(1,size+1,1))
            try:
                p = random.choice(list(self.mutation.keys()))
            except:
                p =  random.choice(self.include)
            """get the coordinates"""
            left,right = (range(self.coordinates[p]),range(self.coordinates[p],len(self.p_space[p]),1))
            """sign of the mutation"""
            try:
                sign = self.mutation[p] 
            except:
                sign = random.choice([-1,1])
            if p in self.pareto_parameters:
                sign = random.choice([-1,1])
            """move to the left or right"""
            if sign < 0:
                if len(left) < step:
                    if left:
                        step = random.choice(left)
                    else:
                        step = 0
                mutations.append((p,self.coordinates[p] + (step*sign)))
                
            else:
                if len(right) < step:
                    if right:
                        step = random.choice(range(len(right)))
                    else:
                        step = 0
                mutations.append((p,self.coordinates[p] + (step*sign)))     
                
        """overwrite the coordinates and position by moving"""
        for p,c in mutations:
            if c == self.spacing:
                c -= 1
            self.coordinates[p] = c
            self.position[p] = self.p_space[p][c]

        """update the tracker"""
        self.track[it] = copy.copy(self.position)
        self.coordinate_track[it] = copy.copy(self.coordinates)
                        
    def update_mutation_probability(self,sensitivities,it):
        """set the new sensitivities"""
        self.forward_sensitivity[it] = sensitivities
        """this is the time of the last update"""
        self.last_forward_update = it
        """senstivity of the parameter"""
        self.factor = {i:[0,0] for i in self.position.keys() if i not in self.fixed}
        """sensitivity to states"""
        for number,m in self.measurements.items():
            """observables and data"""
            for observable,data in m.profile.items():
                lsq = data-self.simdata[number][observable]
                """loop through all the sensitivities and add number to factor"""
                for parameter in self.position.keys():
                    try:
                        fct = sum(lsq*self.forward_sensitivity[it][number][parameter][observable])
                        if fct <= 0:
                            self.factor[parameter][0] += fct
                        else:
                            self.factor[parameter][-1] += fct
                    except KeyError:
                        """the parameter is not in the model thus is cannot count towards the factor"""
                        pass
                        
        restricted = {}
        """the boundaries and the parameters"""
        for parameter,boundaries in self.boundaries.items():
            """"the lower and upper boundaries"""
            lower,upper = boundaries
            """"check and sort restricted"""
            restricted[parameter] = float(upper/lower)
        """sensitivity parameters"""
        self.mutable_parameters, sensitivity = zip(*self.factor.items())
        """"the negative and positive sensitivities"""         
        neg,pos = zip(*sensitivity)
        """we need to sort the parameters"""
        #######################################################################
        rneg = ranklist(numpy.array(neg))
        rpos = ranklist(numpy.array(pos))
        #######################################################################        
        """update the list of mutatory candidates"""
        self.mutation_candidates = []
        
        """mneg and mpos"""
        self.mneg,self.mpos = [],[]
        for i in range(int(len(self.mutable_parameters)/2)): 
            factors = self.factor[self.mutable_parameters[rneg[i]]]
            if abs(factors[-1]) < abs(factors[0]):
                self.mneg.append((-1,self.mutable_parameters[rneg[i]]))
        for i in list(reversed(range(  int(len(self.mutable_parameters)/2),len(self.mutable_parameters),1))):
            factors = self.factor[self.mutable_parameters[rpos[i]]]
            if abs(factors[-1]) > abs(factors[0]):
                self.mpos.append((1,self.mutable_parameters[rpos[i]]))
        
        """set mutation candidates and prune them"""
        self.mutation_candidates = self.mpos + self.mneg
        
        """"these can be deleted in the beginning if there is an impact"""
        self.pareto_parameters = []
        """count the number of parameters with a paretolike impact"""
        for sign,i in self.mutation_candidates:
            neg,pos = self.factor[i][0],self.factor[i][1]
            if neg == 0. and pos > 0:
                pass
            elif neg < 0. and pos == 0.:
                pass
            else:
                self.pareto_parameters.append(i)
                   
        """the boundaries that have no impact on actually changing rate"""
        self.mutation_candidates = [(s,i) for s,i in self.mutation_candidates if restricted[i] >= 1.5]

        self.mutation = {}
        """increase the range and update sign"""
        for i in range(len(self.mutation_candidates)):
            sign,parameter = self.mutation_candidates[i]
            self.mutation.update({copy.copy(parameter):copy.copy(sign)})

    def update_height(self,scores,simdata,it,sf = "standardeviation"):
        ''' This function updates the score of the agent
        it tracks if the score and or position will be updated'''
        y = simdata
        if not it:
            self.height,self.fittest,self.simdata,self.fitset = scores,it,y,self.position
        """loop through the list of the heights and assess if scores are better than previous scores"""
        cpr = [i for i in range(len(scores)) if scores[i].score < self.height[i].score]
        """"assess if it is an accepted movement if so, update current position agent"""
        self.accepted = False
        if len(cpr) >= 1:        
            self.accepted = True
            self.height   = scores
            self.fittest  = it
            self.simdata  = y
            self.SD       = self.height[0].score
            self.fitset   = self.track[self.fittest]

        """re-initialize the position and coordinate system"""
        self.position    = self.track[self.fittest]
        self.coordinates = self.coordinate_track[self.fittest]
        """accepted movement which needs to maintain"""
        self.accepted_movements.append(copy.copy(self.accepted))
          
    def convergence_track(self,threshold = -3):
        """default is false"""
        self.forward = True
        if len(self.accepted_movements) > 3:
            self.boolean = self.accepted_movements[threshold:]
            """assess the boolean of the list"""
            self.forward = True
            if any(self.boolean) == True:
                self.forward = False
            if self.fittest <= self.last_forward_update:
                self.forward = False
                
class Optimization:
    def __init__(self,measurements,
                 #measurements and parameters that can be evolved over the course of the optimization
                 include     = [],
                 #number of generations and agents per generation
                 generations = 250,
                 agents = 10,
                 #time of the simulation
                 time = (0,100),
                 dt   = 1,
                 #startsample number for multistart
                 startsamples = 100,
                 sobol = False,
                 forward = False,
                 #integerization of the parameter space
                 probability_random_mutation = 0.25,
                 startset = [],
                 sf = "leastsquaresfit"):   
        
        '''to estimate the parameters you need the model the measurement you wish to fit to and 
        the parameters you wist to include i.e. those that can be mutated as an input'''
        self.measurements = {i:measurements[i] for i in range(len(measurements))}
        self.models       = {i:self.measurements[i].model for i in range(len(self.measurements))}
        
        """combined parameter space"""
        p_space = {}
        """create a cohesive parameter space"""
        for model in self.models.values():
            for k,v in model.p_space.items():
                p_space[copy.copy(k)] = copy.copy(v)
            
        '''store data of the fittest mutations'''
        self.data = {}
        '''The multistart function randomly samples n number of models using sobol sequence
        and solves them to give the appropriate starting conditions for modelfitting'''
        
        if startset == []:
            multistart = build_global_search(p_space,include = include,samples = startsamples, order = int(math.log10(model.spacing)),sobol = False)   
            """multistart for the system"""
            explore = multimodalstart(multistart,self.measurements,include,sf = sf)
            """initial position"""
            ipos  = [multistart[i] for i in select(explore)]
        else:
            ipos = [{k:v for k,v in random.choice(startset).items() if k in p_space.keys()}]
            ipos.append(copy.copy(ipos[0]))
            for n in ipos:
                for k,v in n.items():
                    arr = (numpy.array([p_space[k][i] for i in range(len(p_space[k]))]) - v)**2
                    minimum = numpy.argmin(arr)
                    n[k] = p_space[k][minimum]    

        """This gives us an output of n number of parents these will be recombined creating a total of n children + parents""" 
        self.agents = {i:Agent(self.measurements,p_space,include,pos = recombine(ipos)) for i in range(agents)}
        for i in range(len(ipos)):
            self.agents[i] = Agent(self.measurements,p_space,include,pos = ipos[i])
        
        """convergence of the system"""
        self.convergence = Convergence()
        """calculate the height of the agents at position the 0th iteration"""
        scores = {}
        import time as trackseconds
        for i in range(len(self.agents)):
            """scores and simutions"""
            start = trackseconds.time()
            scores[i],simdata,fwd = simulate_measurements(self.measurements,self.agents[i].position,sf = sf)
            end = trackseconds.time()
            print(start-end)
            """update the parameter vector and score for vector"""
            self.agents[i].update_position(self.agents[i].coordinates,0)
            """"the initial height of the agents in the landscape"""
            self.agents[i].initial_height(scores[i],simdata)
        """check initial convergence"""
        terminate,self.agents,converge = self.convergence.update_convergence(self.agents,0)
          

        optimization_time = timesystem.time()
        """loop through the generation and update the agent structure which contains the parameter sets which are needed """
        for it in range(generations):
            """update the current coordinate position of all the agents cpos with the fittest position fpos"""
            scores = {}
            """if the iteration is zero we need to instantiate agents from multifit"""
            if not it:
                cpos = [self.agents[i].coordinates for i in range(len(ipos))] 
                
            for move in range(len(self.agents)):
                """turn of forward sensitivities"""
                if forward == False:
                    self.agents[move].forward = False

                """update the position of the agent through recombination and mutation"""
                if random.random() < 0.25:
                    self.agents[move].update_position(recombine(cpos),it)
                if random.random() < probability_random_mutation:
                    self.agents[move].random_mutation(it)
                else:
                    if it != 0:
                        self.agents[move].micromutation(it)
                """simulate and score the measurements for this agent with this simcscore function"""
                scores[move],simdata,forward_vector = simulate_measurements(self.measurements,self.agents[move].position,forward = self.agents[move].forward, sf = sf)
 
                """update position of the agents in the fitness landscape"""
                self.agents[move].update_height(scores[move],simdata,it,sf = sf)
                """assess the overall convergence and determine more sensitivity info is needed"""
                if self.agents[move].forward:
                    self.agents[move].update_mutation_probability(forward_vector,it)
                    
                """the convergence is needed """
                self.agents[move].convergence_track()
            
            """check the convergence of the agents as a group"""
            terminate,self.agents,converge = self.convergence.update_convergence(self.agents,it)
            if terminate and it > 150:
                break

            '''select n parents to be selected based on rank based selection''' 
            cpos = [self.agents[i].coordinates for i in select(scores,select = len(cpos))]
            self.optimalset = copy.copy(self.convergence.position[-1])

        """ find and name the fittest agent in themic""" 
        self.fittest = self.agents[select({i:self.agents[i].height for i in range(len(self.agents))},select = 1)[0]]
 
        """calculate the time of the optimization"""
        ft = (timesystem.time()-optimization_time)/60               
        print("Estimation Accomplished optimization statistics are \n \t \t score {0} \n \t \t opimization time {1} \n \t \t  genaration  {2} ".format(self.convergence.lowest_SD,ft,it))   
        
    def gif(self,):
        """the path of the desktop folder"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')    
        """directory where the files are stored"""
        directory = desktop + "\\__FitGif__\\"
        if not os.path.exists(directory):            
            os.makedirs(directory)    
        filenames = []

        """create gif folder"""
        gif_folder = directory +"Gif_{}\\".format(str(filecount(directory))) 
        if not os.path.exists(gif_folder):            
            os.makedirs(gif_folder)          
    
        """"Loop through the generations and fix this"""
        for generation in range(len(self.convergence.simdata)):
            """plot the data and the resulting fits"""
            fig = plt.figure(figsize=(14,7))
            
            """subplot 1"""
            ax = plt.subplot(1,2,1)
            """measurements and convergence to measurement"""
            for n in range(len(self.measurement)):
                """"the observables of the measurement"""
                for state in self.measurement[n].observables:
                     """the lsq measurements"""
                     if n == 0:
                         ax.plot(self.convergence.simdata[generation][n][state],label = state)
                         ax.scatter(self.measurement[n].time,self.measurement[n].profile[state],alpha = 0.08)
                     else:
                         ax.plot(self.convergence.simdata[generation][n][state])
                         ax.scatter(self.measurement[n].time,self.measurement[n].profile[state],alpha = 0.08)
     
                ax.set_xlabel("Time")
                ax.set_ylabel("Concentration")
                ax.legend(fancybox = True)
                """set it to be a multiple of the largest measurement scale"""
                largest = max([max(i) for i in self.measurement[n].profile.values()])
            """set the yscale proportionate to the measurement"""
            ax.set_ylim([0,largest*1.25])
            
            """subplot 2"""                
            ax = plt.subplot(1,2,2)
            """scatter the plot of the estimates"""
            ax.plot(self.convergence.iteration[0:generation],self.convergence.scores[0:generation])
            ax.scatter(self.convergence.iteration[0:generation],self.convergence.scores[0:generation],s = 100,color = 'DarkBlue')
            ax.legend(fancybox = True)
            ax.set_xlabel("Algorithm Iteration",size = 14)
            ax.set_ylabel("Average Standard Deviation",size = 14)
            ax.set_xlim([0,max(self.convergence.iteration)+1])
            ax.set_ylim([0,max(self.convergence.scores)+1])
            """ensure the layout is not altered"""
            plt.tight_layout()
            """the path of the stored images, needs to be in gif specific directory"""
            path = gif_folder + str(filecount(gif_folder))
            filenames.append(copy.copy(path))
            """save the figure to the folder"""
            plt.savefig(path,dpi = 460)
            
        """create the gif"""
        import imageio
        images = []
        for filename in filenames:
            images.append(imageio.imread(filename+".png"))
        imageio.mimsave(gif_folder + "movie.gif", images,duration = 0.8)
        
    def show(self,):
        sns.set()
        """plot the data and the resulting fits"""
        fig = plt.figure(figsize=(14,7))
        """subplot 1"""
        ax = plt.subplot(1,2,1)
        """measurements and convergence to measurement"""
        for nmbr,measure in self.measurements.items():
            for state in measure.observables:
                clr = [random.random() for i in range(3)]
                ax.scatter(range(len(measure.profile[state])),measure.profile[state],label = state, color = clr)
                ax.plot(self.convergence.simdata[-1][nmbr][state],color = clr)
                ax.set_xlabel("Time",size = 14)
                ax.set_ylabel("Concentration",size=14)
                ax.legend(fancybox = True)   

        """subplot 2"""                
        ax = plt.subplot(1,2,2)
        """plot the change the scores over the number of iterations"""
        ax.plot(self.convergence.iteration,self.convergence.scores)
        ax.legend(fancybox = True)
        ax.set_xlabel("Algorithm Iteration",size = 14)
        ax.set_ylabel("Average Standard Deviation",size = 14)
        plt.tight_layout()
        plt.show()
            

class BoxplotOptimization:
    def __init__(self,
                 measurements,
                 #measurements and parameters that can be evolved over the course of the optimization
                 include     = [],
                 name        = '',
                 fixed       = {},
                 #number of generations and agents per generation
                 generations = 250,
                 agents = 10,
                 #time of the simulation
                 time = (0,100),
                 dt   = 1,
                 #startsample number for multistart
                 startsamples = 1000,
                 sobol = False,
                 forward = False,
                 #integerization of the parameter space
#                 """number of optimizations in the box"""
                 probability_random_mutation = 0.25,
                 optimization_number = 20,
                 constrained = False,
                 sf = "standardeviation",
                 startset = [],
                 storedata = ''):   
        
        '''to estimate the parameters you need the model the measurement you wish to fit to and 
        the parameters you wist to include i.e. those that can be mutated as an input'''
        self.measurements = measurements
        self.models       = [self.measurements[i].model for i in range(len(self.measurements))]
        """store the number of optimizations"""
        self.optimizations = []
        
        """the path of the desktop folder"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')    
        """directory where the files are stored"""
        self.directory = desktop + "\\__SimulatedExperiments__\\"
        """include this in the optimization other parameters are not mutated"""
        """the path of the desktop folder"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')    
        """directory where the files are stored"""
        self.directory = desktop + "\\__SimulatedExperiments__\\"
        """include this in the optimization other parameters are not mutated"""
        self.include,self.control = [],[]
        for i in self.models:
            for j in i.control:
                self.control.append(copy.copy(j))
            for j in i.boundaries.keys():
                if j not in self.control:
                    self.include.append(copy.copy(j))

        """include the values from the database"""
        if include != []:
            self.include = include

        """overwrite the fixed values"""
        if fixed != {} and type(fixed) != list:
            self.include = [i for i in self.include if i not in list(fixed.keys())]
            for pID, value in fixed.items():
                for m in range(len(self.measurements)):
                    try:
                        self.measurements[m].model.fixed[pID] = value
                    except:
                        print("parameter not present in model")
                       
        if type(fixed) == list:
            self.include = [i for i in self.include if i not in fixed]

        print(fixed)
        """reduce listsize"""     
        self.include = list(set(self.include))
                     
        if storedata != '':
            """store the optimization in the right place!"""
            fit    = storedata + '\\Parameters\\'
            if not os.path.exists(fit):
                os.makedirs(fit)
            newfit = storedata + '\\New Parameters\\'
            if not os.path.exists(newfit):
                os.makedirs(newfit)
            simdata = storedata + '\\Simdata\\'
            if not os.path.exists(simdata):
                os.makedirs(simdata)
            msnt = storedata + '\\Measurement\\'
            if not os.path.exists(msnt):
                os.makedirs(msnt)
                
            """get model information"""
            modelinformation = {'Models':[]}
            """which models are used"""
            boundaries       = {}
            for i in range(len(measurements)):
                modelinformation['Models'].append(measurements[i].model.name)
                for ID,bnd in measurements[i].model.boundaries.items():
                    boundaries[ID] = bnd
            """the pandas dataframe"""     
            df = pandas.DataFrame(boundaries)
            df.to_csv(fit + 'boundaries' + ".txt", sep='\t', na_rep="none")        
            df = pandas.DataFrame(modelinformation)
            df.to_csv(fit + 'models' + ".txt", sep='\t', na_rep="none")        
                    
                    
        tstart = timesystem.time()
        for opt in range(optimization_number):
            """optimization of model to simulated measuremement or real measurement"""
            optimization = Optimization(self.measurements,
                 #measurements and parameters that can be evolved over the course of the optimization
                 include     = self.include,
                 #number of generations and agents per generation
                 generations = generations,
                 agents = agents,
                 #time of the simulation
                 time = time,
                 dt   = dt,
                 #startsample number for multistart
                 startsamples = startsamples,
                 sobol = sobol,
                 forward = forward,
                 startset = startset,
                 #integerization of the parameter space
                 probability_random_mutation= probability_random_mutation,
                 sf = sf)

            # optimization.show()
            """append to list"""
            # self.optimizations.append(optimization)
            """open the relevant files and check if the data is already stored somewhere"""
            values = {}

            if storedata != '':
                score      = optimization.convergence.lowest_SD
                parameters = optimization.optimalset
                parameters['score'] = score
                """optimization data does not yet exist, please check again"""
                with open(newfit + '{}.pickle'.format(str(len(os.listdir(newfit)))), 'wb') as handle:
                    pickle.dump(parameters, handle, protocol=pickle.HIGHEST_PROTOCOL)   
            
            if storedata != '':
                try:
                    path = fit + 'parameters.csv'
                    dataframe = pd.read_csv(path)
                    t_value = dataframe.to_dict()
                    
                    """the score and the convergence"""
                    score      = optimization.convergence.lowest_SD
                    parameters = optimization.optimalset
                    print(parameters,optimization.convergence)
                    
                    for pID,value in fixed.items():
                        parameters[pID] = value
                        
                    """the t value are the fitted """
                    for ID,value in parameters.items():
                        t_value[ID][len(t_value[ID])] = value
                    t_value['score'][len(t_value['score'])] = score
                    data = pandas.DataFrame(t_value)
                    data.to_csv(fit + 'parameters.csv', index = False)
                    print(parameters,'append')

                except:
                    score      = optimization.convergence.lowest_SD
                    parameters = optimization.optimalset
                    print(parameters,'pretest')
                    parameters['score'] = score
                    
                    """optimization data does not yet exist, please check again"""
                    data = pandas.DataFrame(LDtoDL([parameters]))
                    data.to_csv(fit + 'parameters.csv', index = False)
                    print(parameters,'except')
                
            if opt == 0:
                print("""{} minutes per optimization""".format(str((timesystem.time()-tstart)/60.)))
                print("""{} (hours) time remaining for result""".format(str((timesystem.time()-tstart)*optimization_number/3600.)))
                
            del optimization
            
        """"List the overall number of generator objects""" 
        self.data  = LDtoDL([i.convergence.position[-1] for i in self.optimizations])
        self.score = [self.optimizations[i].convergence.lowest_SD for i in range(len(self.optimizations))]
        self.likelihood = {i:{} for i in self.data.keys()}
        for parameter,values in self.data.items():
            self.likelihood[parameter] = {i:sorted(values)[i] for i in range(len(values))}

        """create full dataframes, all the data, fittest data and """
        self.dataframe = pd.DataFrame.from_dict(self.data,orient='index').transpose()  
        self.dataframe['score'] = self.score
     
        """the data of the optimizations"""
        self.simdata = [pd.DataFrame.from_dict(self.optimizations[i].convergence.simdata[-1],orient='index').transpose() for i in range(len(self.optimizations))]
        
        """the means and deveation for the info"""
        self.mean  = {k:numpy.mean(v) for k,v in self.data.items()}
        self.std   = {k:numpy.std(v) for k,v in self.data.items()}
        self.var   = {k:numpy.var(v) for k,v in self.data.items()}
        """pickle objects"""
        self.pickle= (self.data,self.score)
        
    def analyse(self,name = '',directory = '',store = True,show = False):
        index = self.score.index(min(self.score))
        m = [i.convergence.position[-1] for i in self.optimizations]
        maximum = m[index]
        """sort the keys from big to small"""
        sortkeys = []
        """the boundaries and the parameters"""
        for i in range(len(self.measurements)):
            for parameter,boundaries in self.measurements[i].model.boundaries.items():
                if parameter in self.include:
                    """"the lower and upper boundaries"""
                    lower,upper = boundaries
                    """"check and sort restricted"""
                    sortkeys.append((str(float(upper/lower)),parameter))
        sortkeys = list(set(sortkeys))
        val,ID  = zip(*sortkeys)
        sortval = numpy.argsort([eval(i) for i in val])
        
        """create column identifiers"""
        columns = []
        for i in reversed(sortval):
            columns.append(ID[list(sortval).index(i)])
        columns = [i for i in columns if i in self.data.keys()]
        
        
        """"create the dataframe"""
        dataframe = pd.DataFrame.from_dict(self.data,orient='index').transpose()
        if directory != '':
            self.directory = directory

        """"create the figure and subplot"""
        f, axes = plt.subplots(1, 2,figsize=(14,7))
        """plot the results"""
        sns.set(style="whitegrid")
        """First plot"""
        sns.boxenplot(data=dataframe,ax = axes[0],order=columns)
        sns.stripplot(data=dataframe,size=4, jitter=True, color="gray",ax = axes[0],order=columns)
        axes[0].set_ylabel("Parameter Value",size = 14)
        axes[0].set_yscale('log')
        """rotate second labels"""
        for tick in axes[0].get_xticklabels():
            tick.set_rotation(90)
        
        """Second plot"""        
        sns.distplot(self.score,ax = axes[1],hist  = True)
        axes[1].set_ylabel("Convergence Quality",size = 14)
        axes[1].set_xlim([0,1.2*max(self.score)])       
        axes[1].set_xlabel("Average Standard Deviation Data-Model Fit",size = 14)
        plt.tight_layout()
        """store the figure in a directory"""
        if store:
            plt.savefig(self.directory + "Uncertainty" + name + ".png",dpi = 600)  
        if show:
            plt.show()
        plt.close()
         
        """create the dataframe"""
        median = numpy.median(self.score)
        """get indices"""
        indices = [i for i in range(len(self.score)) if self.score[i] < median]
        """parameters"""
        parameters = [i.convergence.position[-1] for i in self.optimizations]
        self.median_data = LDtoDL([parameters[i] for i in indices])
        dataframe = pd.DataFrame.from_dict(self.median_data,orient='index').transpose()
        """create the figure and subplot"""
        f, axes = plt.subplots(1, 2,figsize=(14,7))
        """plot the results"""
        sns.set(style="whitegrid")
        
        """First plot"""
        sns.boxenplot(data=dataframe,ax = axes[0],order=columns)
        sns.stripplot(data=dataframe,size=4, jitter=True, color="gray",ax = axes[0],order=columns)
        dataframe = pd.DataFrame.from_dict(maximum,orient='index').transpose()
        sns.stripplot(data=dataframe,size=4, jitter=True, color="red",ax = axes[0],order=columns)
        axes[0].set_ylabel("Parameter Value",size = 14)
        axes[0].set_yscale('log')
        """rotate second labels"""
        for tick in axes[0].get_xticklabels():
            tick.set_rotation(90)
        
        """Second plot"""        
        sns.distplot([self.score[i] for i in indices],ax = axes[1],hist  = True)
        axes[1].set_ylabel("Convergence Quality",size = 14)
        axes[1].set_xlim([0,1.2*max([self.score[i] for i in indices])])
        axes[1].set_xlabel("Average Standard Deviation Data-Model Fit",size = 14)
        plt.tight_layout()
        """store the figure in a directory"""
        if store:
            plt.savefig(self.directory + "UncertaintyOutlier" + name + ".png",dpi = 600)  
        if show:
            plt.show()
        plt.close()        
        """plot the fit to the data"""
        fig = plt.figure(figsize=(14,7))
        colors = [numpy.random.rand(3) for i in range(100)]
        sns.set(style="whitegrid")
        for i in range(len(self.optimizations)):
            data = self.optimizations[i].convergence.simdata[-1]
            """set a count to ensure each state gets the right color"""
            count = 0
            for n in range(len(data)):
                for state,simulation in data[n].items():
                    if state in self.measurements[n].profile.keys():
                        c = colors[count]
                        if i == 0:
                            plt.plot(simulation,label = state,color = c)
                        else:
                            plt.plot(simulation,color = c)                            
                        plt.scatter(range(len(self.measurements[n].profile[state])),self.measurements[n].profile[state],color = c,alpha = 1)
                        count += 1
        """plot the legend and set a tight layout"""
        plt.legend(fancybox = True)
        plt.tight_layout()
        """store the figure in a directory"""
        if store:
            plt.savefig(self.directory + "Fits" + name + ".png",dpi = 600)  
        if show:
            plt.show()
        plt.close()
            

            
            
    