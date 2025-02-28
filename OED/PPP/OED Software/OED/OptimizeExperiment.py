# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:17:35 2019

@author: huckg
"""

""""turn yur script in an autorun module that way you can combine python 2 and 3 in 1 program
Note that your boundaries cannot be zero or fixed, you can set your fixed to zero if needed and not add them
to the "include" variable"""
import time
import numpy
import random
import math
import matplotlib.pylab as plt
import os
import seaborn as sns
import copy

from PlotData import *
from Parse import *
from Operations import *
from Measurements import *
from OptimizationOperations import *
from ExperimentOptimizationOperations import *
from Scoring import __score__

         
class Convergence:
    def __init__(self,):
        """the highest score"""
        self.collinearity = 1*10**25
        """the best simdata"""
        self.simdata      = []
        """iteration where the best simdata happened"""
        self.iteration    = []
        """best coordinates"""
        self.coordinates  = []
        """scores"""
        self.scores       = []
        """converge"""  
        self.converge     = []

    def update_convergence(self,agents,it,threshold = -3):
        """agents and the agent number ID"""
        nmbr,agents = zip(*agents.items())
        """minimal values"""
        collinearities = [i.average for i in agents]
        """smallest value"""
        fittest_agent = min(collinearities)
        """get index of smallest"""
        index = nmbr[collinearities.index(fittest_agent)]
        if fittest_agent < self.collinearity:
            self.collinearity = fittest_agent
            """update simdata, coordinates, position, iteration etc."""
            self.iteration.append(it)
            self.scores.append(fittest_agent)
            self.coordinates.append(agents[index].coordinates)
            self.simdata.append(agents[index].simdata)
            self.fittest_agent = copy.copy(agents[index])
                    
        """there are 3 miminum termination criteria, no improvement in 20 generations
           atleast 30 generations
           the fit needs to be atleast 0.1 devation on average if you want to stop earlier"""
        terminate = False
        if it > 10:
            """set two different termination criteria"""
            if len(self.iteration) > 2:
                """no improvement over 100 generations = no movement"""
                if self.iteration[-1] - self.iteration[-2] > 20:
                    terminate = True
        """modify converge"""
        self.converge.append(self.collinearity)
        """redict the agents"""
        agents = {nmbr[i]:agents[i] for i in range(len(agents))}
        return terminate,agents,self.converge
            
class ExperimentAgent:
    def __init__(self,
                 model,
                 #Coordinates and pulsespace
                 coordinates,
                 pulsespace,
                 #the initial conditions where stuff is not pulsed and or changed
                 initial     = {},
                 #observables and include
                 observables = [],
                 include     = [],
                 #parameters in pulsepattern that cannot be mutated
                 fixed       = [],
                 #time and coordinate spacing of these things
                 time = (0,200),
                 pulsestart = 50,
                 spacing = 4,
                 pindex  = 10,
                 #show the pattern
                 show = True):
        
        """the model"""
        self.model = model
        """assign the coordinates"""
        self.coordinates = coordinates
        """definition of the pulse space"""
        self.pulsespace = pulsespace
        """if not include than all parameters can be mutated"""
        if not include:
            include = list(self.coordinates.keys())
        self.control_interval = list(list(self.coordinates.values())[-1].keys())
        
        """ the state variables we are scoring against"""
        self.observables = observables
        self.include = include
        self.fixed = fixed
    
        """time vector"""
        self.start,self.end = time
        """control index"""
        self.pindex = pindex
        """"first control parmameters change"""
        self.pulsestart = pulsestart
        '''these are the continuously updated parameter sets, its corresponding coordinates
        and the score sequence. They are only updated if the score is better'''
#        self.initial = initial
#        if not initial:
#            for parameter,times in self.coordinates.items():
#                self.initial[parameter] = {(self.start,self.pulsestart):random.choice(range(self.pindex))}
#                self.coordinates[parameter].update(self.initial[parameter])
        """get pulsesequence"""
        self.pulsesequence = list(sorted(list(self.coordinates.values())[-1].keys()))
        """set a selection treshold"""
        self.selection_treshold = 1
        """set the overall number of mutations by pulse size we 
        state that a mutation can be upto 25% of the space"""
        self.msize = int(len(self.pulsesequence)*random.choice([0.1+(i*0.01) for i in range(2)]))
        """initial height"""
        self.height = None
        """coordinate spacing"""
        self.spacing = spacing
        '''information of the EA et al. is still information on the kinetics
        of the model and range of potential behavior thus we store the data
        and the transition moreover we can track their status over time'''
        self.track,self.coordinate_track = {},{}
        '''the list of accepted mutations/recombinations'''        
        
    def initial_height(self,scores,simdata):
        """"this is the initial height of the system as obtained by the
        monte carlo before we start with the optimizaton routine"""
        self.height  = scores  
        self.simdata = simdata
        """"these will get overwritten (probably) as the evolution starts"""
        self.fittest = 0
        """set initial position of the agent"""
        self.coordinate_track[self.fittest] = self.coordinates
        self.average = self.height[-1].score

    def update_position(self,coordinates,it):
        self.coordinates = coordinates
        """track the crd"""
        self.coordinate_track[it] = coordinates
        
    def mutate(self,it,attr = [],mutation_number = 2):
        self.mutationlist = []
        """number of mutations"""
        number_of_mutations = random.choice(range(1,mutation_number + 1,1))
        for i in range(number_of_mutations):
            """randomly chosen size"""
            size = random.choice(range(self.msize))
            """parameter selected for mutation"""
            parameter = random.choice(self.include)
            """time of the pulse"""
            pulse = random.choice(range(len(self.pulsesequence)))
            """pulse go left or right"""
            chosensequence = []
            if random.random() <= 0.5:
                """go left along the number line [(),(),()...<- (chosen),(),(),()]"""
                start_index = pulse-size
                if start_index < 0:
                    start_index = 0
                for t in range(start_index,pulse,1):
                    try:
                        chosensequence.append(self.pulsesequence[t])
                    except IndexError:
                        print("mutation is failing")
                        break
            else:
                """go right along the number line [(),(),(),(chosen) ->...,(),(),()]""" 
                end_index = pulse+size
                if end_index > len(self.pulsesequence)-1:
                    end_index = len(self.pulsesequence)
                for t in range(pulse,end_index,1):
                    try:
                        chosensequence.append(self.pulsesequence[t])
                    except IndexError:
                        print("mutation is failing")
                        break
            """size of the index space"""
            size = random.choice(range(self.pindex))
            """append mutation"""    
            if len(chosensequence) != 0:
                self.mutationlist.append((parameter,chosensequence,size))
            
        self.index = True
        """tracking mutations"""
        if random.random() < 0.25:
            self.index = 0
        for parameter,sequence,index in self.mutationlist:
            if self.index == True:
                self.index = index
            for i in sequence:
                """update the coordinate vector"""
                self.coordinates[parameter][i] = copy.copy(self.index)
        
    def update_height(self,scores,simdata,it):
        """need to update the tracker everytime you change the coordinates"""
        self.coordinate_track[it] = copy.deepcopy(self.coordinates)
        ''' This function updates the score of the agent
        it tracks if the score and or position will be updated'''
        y = simdata
        if it == 1:
            self.height  = scores
            self.fittest = it
            self.simdata = y
            
        """if the mutation is accepted than the vector need to be updated"""
        cpr = [i for i in range(len(scores)) if scores[i].score < self.height[i].score] #CHEKC THIS LINE, the list is ordered but created from dict.items() call
        """assess if mutation is accepted"""
        self.accepted = False
        if len(cpr) >= self.selection_treshold:
            self.accepted = True
            self.height  = scores
            self.fittest = it
            self.simdata = y
            self.average = self.height[-1].score
        """position needs to be updated"""
        self.coordinates = self.coordinate_track[self.fittest]
        

class OptimizeInformation:
    def __init__(self,model,
                 #wchich parameters are pulsed
                 time_dependent_parameters = [],
                 initial_control_conditions= {},
                 conditions                = {},
                 #what is measuremd
                 observables = [],
                 #number of generations and agents per generation
                 generations = 30,
                 agents = 5,
                 #length of pulse time of experiment, start of sequence number of indices
                 time = (0,300),
                 dt =1,
                 pulsestart = 10,
                 plength = 20,
                 #size of control region
                 pindex = 10,
                 #number of start samples
                 multistart = 10,
                 #scoring functions
                 sf = "D_Fisher"): 
        

        """"define model"""
        self.model = model
        """find control points"""
        self.control = time_dependent_parameters
        if len(self.control) == 0:
            self.control = model.control
        self.observables = observables
        if len(self.observables) == 0:
            self.observables = model.observables
        start,end = time
        self.pindex = pindex
        """set pulse space and initial_pulse_vector"""
        self.pulsespace = {}
        """"here we will generate a number of spaces where whe can sample and create
        a input space in a MIMO system, we subsequently define some functions such as mutate
        funciton to allow this input patternt to evolve, the class needs to maintain an
        input structure i.e. is it on or off and what size the input should have followed by the 
        lenght of the input, there is a super switch which monitors if a pulse is active or not
        recombination exchanges vectors"""
        """
#        self.initiation = {}
#        if len(initial_control_conditions) > 0:
#            self.initiation =  {(start,pulsestart):initial_control_conditions}"""
        pulsevector = {}
        for i in self.control:
            controlbounds = {i:model.boundaries[i]}
            """define boundaries of the pulse space"""
            size,exclude  = define_parameter_boundaries(controlbounds,spacing = pindex) 
            """update size of pulse space"""
            self.pulsespace.update(size)
            if not self.pulsespace:
                raise "The model needs to have defined control parameters for this to work"
            """"size of the pulse being given"""  
            pulsevector[i] = {(j,j+plength):int(pindex/2) for j in range(pulsestart,end,plength)}
            pulsevector[i].update({(start,pulsestart):int(pindex/2)})
            
        """pulse space transferred to the model"""
        self.model.pulsespace = self.pulsespace
        """generate random pulse sequence"""
        c_set = {}
        for i in range(multistart):
            coordinates     = copy.deepcopy(pulsevector)
            new_coordinates = gerenate_mutated_pattern(coordinates,pindex,mutation_number = 1)
            c_set[i]        = copy.copy(new_coordinates)
    
        """self.agents contains all the pulsepattern objects"""
        self.agents = {}
        for i,c in c_set.items():
            self.agents[len(self.agents)] = ExperimentAgent(self.model,c,self.pulsespace)

        """convergence tracker"""
        self.convergence = Convergence()

        scores ={}
        import time as trackseconds
        """the number of agents which can evolve patterns"""
        for index,agent in self.agents.items():
            print(index)
            """simulate the data"""
            score,simdata,forward = simulate_experiment(self.model,dt = dt,conditions = conditions,time = time,time_dependent_parameters = agent.coordinates,forward = True, sf = sf)
            """update the scores"""
            scores[index] = score
            """set initial position in landscape"""
            self.agents[index].initial_height(score,simdata)
        """update convergence"""
        self.convergence.update_convergence(self.agents,0)
            
        """agents selected to run the simulation"""
        self.agents = {i:self.agents[select(scores,select = agents)[i]] for i in range(agents)}
        
        """"Loop through the generation and evolve the pattern"""
        for it in range(1,generations,1):
            print(it, 'Generation')
            """define scores"""
            scoreset = {}
            """update the fittest members of the generation"""
            if it == 1:
                cpos = [self.agents[i].coordinates for i in [0,1]]
            for move in range(len(self.agents)):
                print(move)
                scores = []
                """move agents by updating and mutating the pattern"""
                if random.random() < 0.75:
                    self.agents[move].update_position(recombine_pulsepattern(cpos),it) 
                self.agents[move].mutate(it)
                """define conditions"""
                score,simdata,forward = simulate_experiment(self.model,dt = dt,conditions = conditions,time = time,time_dependent_parameters = agent.coordinates, forward = True,sf = sf)
                scoreset[move] = score
                """update input vector"""  
                self.agents[move].update_height(scoreset[move],simdata,it)
                
            """check the convergence of the agents as a group"""
            terminate,self.agents,converge = self.convergence.update_convergence(self.agents,it)
            if terminate:
                break
            if it > 150:
                if converge[it] == converge[it-150]:
                    break
                
            """select the fittest"""
            cpos = [self.agents[i].coordinates for i in select(scoreset,select = len(cpos))]
        """show the fittest members"""
        self.fittest = self.agents[select({i:self.agents[i].height for i in range(len(self.agents))},select = 1)[0]]
        
        
    def gif(self,):
        try:
            """the path of the desktop folder"""
            desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop')    
            """directory where the files are stored"""
            directory = desktop + "\\__ExperimentalGif__\\"
            if not os.path.exists(directory):            
                os.makedirs(directory)    
            filenames = []
            """create gif folder"""
            gif_folder = directory +"Gif_{}\\".format(str(filecount(directory))) 
            if not os.path.exists(gif_folder):            
                os.makedirs(gif_folder)      
    
    
            new_scores = []
            for i in range(len(self.convergence.scores)):
                new_scores.append(max(self.convergence.scores)-self.convergence.scores[i])
            new_scores = [i/max(new_scores) for i in new_scores]  
    
            sns.set() 
            """"Loop through the generations and fix this"""
            maximum = 0
            for mutation in range(len(self.convergence.simdata)):
                """"simdata and search for highest value"""
                for k,v in self.convergence.simdata[mutation].items():
                    if maximum < max(v):
                        maximum = max(v)         
            
            for mutation in range(len(self.convergence.simdata)):
                """plot the fittest pulsepattern with expected concentration"""
                fig = plt.figure(figsize=(10,5))
    
                """first plot"""            
                ax = plt.subplot(1,2,1)
                
                """"simdata and search for highest value"""
                for k,v in self.convergence.simdata[mutation].items():
                    ax.plot(v,label = k)  
                ax.set_ylim([0,maximum*1.1]) 
                ax.set_xlabel("Time")
                ax.set_ylabel("Concentration")
                plt.legend(fancybox = True)
                
                """second plot"""
                ax = plt.subplot(1,2,2)  
                ax.set_ylim([0,self.pindex+1])
                
                """"time and control parameter input"""
                time,control = {},{}
                for parameter,t in self.convergence.coordinates[mutation].items():
                    seq,coordinate = zip(*t.items())
                    """sort the seq"""
                    seq = list(sorted(seq))
                    if parameter not in control.keys():
                        control[parameter] = []
                        time[parameter]    = []
                    for i in range(len(seq)):
                        """append time"""
                        control[parameter].append(t[seq[i]])
                        """append time vector"""
                        start,end = seq[i]
                        time[parameter].append(start)
                        
                """plot the control parameters"""
                for parameter,value in control.items():
                    ax.plot(time[parameter],value,label = parameter)
                ax.legend(fancybox = True)
                ax.set_xlabel("Time")
                ax.set_ylabel("Control Parameter Value")
                ax.set_title("Information: {}".format(round(new_scores[mutation],4)))
                """the path of the stored images, needs to be in gif specific directory"""
                path = gif_folder + str(filecount(gif_folder))
                filenames.append(copy.copy(path))
                """save the figure to the folder"""
                plt.savefig(path,dpi = 460)
                """update generation"""
                plt.close(fig)
                
            """import the imageio toolbox and make gif"""
            import imageio
            images = []
            for filename in filenames:
                images.append(imageio.imread(filename+".png"))
            imageio.mimsave(gif_folder + "movie.gif", images,duration = 0.8)      
        except:
            print("GIF maker did not work")