# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:26:43 2019

@author: huckg
"""

import numpy
#import tellurium as te
import os
import pickle
import matplotlib.pylab as plt
from Operations import *

class RateLaw:
    def __init__(self,rate,states,parameters):
        """rateterm representing a reaction"""
        self.rate = rate
        """the states present in this ratelaw"""
        self.states = []
        for i in states:
            if i in self.rate:
                self.states.append(i)
        """the parameters present in this rate"""
        self.parameters = [] 
        for i in parameters:
            if i in rate:
                self.parameters.append(i)
        """state index"""
        self.state_indices = [states.index(i) for i in self.states]

class experimental_conditions:
    def __init__(self,model,conditions):
        self.simulated_conditions = copy.copy(conditions)
        """"set intitial of the conditions of the """
        self.initial = {}
        self.control = {}
        for i,j in self.simulated_conditions.items():
            if i not in model.states and i not in model.boundaries.keys():
                del self.simulated_conditions[i]
            if i in model.states:
                self.initial[i] = j
            elif i in model.boundaries.keys():
                self.control[i] = j
                
        """assess whether control parameters that act as initial conditions actually exist"""
        self.initial_control = []
        if model.initial_control:
            parameter,states      = zip(*model.initial_control) #unzip parameters and states
            boolean               = [i for i in states if i in model.states] #check if they exist at all in the current model
            self.initial_control  = [model.initial_control[states.index(i)] for i in boolean] #include them in initial control if they exist
            
    def initial_conditions(self,model,parameters):
        """note this does wqork if your control and initial conditions are the same but 
        seperated such that control supercedes initial conditions, now it is possible to
        set initial conditoins but seperate control variables"""
        ic = numpy.zeros((len(model.states)))
        
        """add the initial conditions reflecting the point if initial control"""
        for parameter,state in self.initial_control:
            ic[model.map[state]] = parameters[model.map[parameter]]
            
        """add the right values to the states"""
        for state,value in self.initial.items():
            ic[model.map[state]] = value 
        return ic 



def solve_mxl(model,condition,observables = [],dt = 0.1,time = 10,show = False):
        conditions = experimental_conditions(model,condition)    
        """"calculate the maximum likelihood of the network starting with the parameters vector"""
        p     = numpy.zeros(len(model.fixed))
        """set the fixed values of the parameters"""
        for fix,value in model.fixed.items():
            p[model.map[fix]] = value
        """the measurement conditions"""
        for pvar,value in conditions.control.items():
            p[model.map[pvar]] = value   
        """initial control contains all the control parameters denoted as states and vice versa"""
        for pvar,state in conditions.initial_control:
            if state in conditions.initial.keys():
                p[model.map[pvar]] = conditions.initial[state] 
        """simulate the network"""
        ic = conditions.initial_conditions(model,p)  
        maximum_likelihood,t  = solve_system(model,p,ic = ic,t = time,dt = dt)     
        if show:
            for i in range(len(model.states)):
                plt.plot(t,maximum_likelihood[:,i],label = model.states[i])
            plt.legend(fancybox = True)
            plt.xlabel("Time")
            plt.ylabel("Concentration")
            plt.title("Simulation of a single user defined parameter set")
            plt.show()
        return maximum_likelihood
    
def find_elements(s, ch):
    return [(i,ch) for i, ltr in enumerate(s) if ltr == ch and s[i+1] == ' ']

def define_ratelaw(s):
    ind = list(sorted(find_elements(s,"+")+find_elements(s,"-")))
    rate_equations = []
    for i in range(len(ind)-1):
        fi,fs = ind[i]    #the index of the first sign + or -
        ni,ns = ind[i+1]  #the index of the second sign + or -
        """stringcut which defines the rate"""
        rate_equations.append(s[fi:ni])
    """the final cut"""
    li,ls = ind[-1]
    """"the rate equations to strip"""
    rate_equations.append(s[li:len(s)])
    return rate_equations
                
def model_matrix_format(equations,states,fixed):
    staterate = {}
    for equation in range(len(equations)):
        eq = equations[equation]
        y  = states[equation]
        """define ratelaws"""
        ratelaws = define_ratelaw(eq)  
        """class of formalized rate laws"""
        formalized_rates = []      
        for rate in ratelaws:
            r,statelist = [],[]
            """"find the sign"""
            sign = rate[0]
            """"loop through parameters and append if in rate"""
            for i in fixed.keys():
                if i in rate:
                    r.append(i)
            """"loop through states and append if in rate"""                
            for i in states:
                if i in rate:
                    statelist.append(i)
            if y not in staterate.keys():
                staterate[y] = []
            staterate[y].append((sign,rate[1:len(rate)]))

    """fluxvector"""
    fluxes = list(set([item.strip() for sublist in staterate.values() for sign,item in sublist]))
    """create stoichiometric and flux vector this can subsequently be used as a
    tool for creating the nessecary tellurium models"""
    S = numpy.zeros((len(states),len(fluxes)))
    """fill up the vector with fluxes"""
    for i in range(len(states)):
        """unzip the fluxes"""
        signs,ratelaws = zip(*staterate[states[i]])
        for j in range(len(ratelaws)):
            """set the flux"""
            flux = ratelaws[j].strip()
            """ratelaws added to stoichiometry"""
            for n in range(len(fluxes)):
                if flux == fluxes[n].strip():
                    if signs[j].strip() == "-":
                        S[i][n] = -1
                    elif signs[j].strip() == "+":
                        S[i][n] = 1    
                        
    """update ratelaws """
    ratelaws = {i:RateLaw(i,states,fixed) for i in fluxes}
    return fluxes,S,ratelaws,staterate

def define_ratelaw(s):
    ind = list(sorted(find_elements(s,"+")+find_elements(s,"-")))
    rate_equations = []
    for i in range(len(ind)-1):
        fi,fs = ind[i]    #the index of the first sign + or -
        ni,ns = ind[i+1]  #the index of the second sign + or -
        """stringcut which defines the rate"""
        rate_equations.append(s[fi:ni])
    """the final cut"""
    li,ls = ind[-1]
    """"the rate equations to strip"""
    rate_equations.append(s[li:len(s)])
    return rate_equations

#_________________Parse a model to antimony and SBML________________*
def antimony_parser(states,parameters,fluxes,S,ratelaws,staterate,fixed,mapping):
    """find all indices"""
    row,column = numpy.nonzero(S)
    """the paired indices of the S matrix"""
    stoch = [(row[i],column[i]) for i in range(len(row))]
    """find all paired fluxes between states etc."""
    connections = {i:[] for i in range(len(fluxes))}
    for i,j in stoch:
        sign = S[i,j]
        connections[j].append((sign,i)) #state to flux
    fluxlist = []
    """single incoming and outgoing reaction"""
    for flux,inv in connections.items():
        for sign,state in inv:
            reaction = ""
            if sign == -1:
                reaction += " " + states[state] + " -> ; " + fluxes[flux]
            elif sign == 1:
                reaction += " -> {} ; ".format(states[state]) + fluxes[flux]
            fluxlist.append(reaction)
    fluxlist = [i.replace("  "," ") for i in fluxlist]
            
    """create a model list"""
    model = ''
    for i in range(len(fluxlist)):
        """add flux equations"""
        model += "J{}: ".format(i) + fluxlist[i] + ' ; \n'
    model += ''
    """add initial conditions"""
    for i in states:
        model += " "+i+" = 0 ; \n"
    model += ''
    
    """add parameters to the system"""
    for k,v in fixed.items():
        model += " " + k + ' = {} ; \n'.format(str(v))
        
    """the y and p vector replacing the states"""
    y = {states[i]:'y{}'.format(mapping[states[i]]) for i in range(len(states))}
    p = {parameters[i]:'k{}'.format(mapping[parameters[i]]) for i in range(len(parameters))}
    
    """overwrite the states and parameters"""
    for k,v in y.items():
        k = " " + k + " "
        model = model.replace(k,v)
    for k,v in p.items():
        k = " " + k + " "
        model = model.replace(k,v)
        
    """replace the delimiters"""
    model = model.replace("**","^")
    model = model.replace("|+|","+")
    model = model.replace("|-|","-")
    return model
  
"""store the SBML file"""
def storeSBML(sbml,name = 'test'):
        """"set the path of the reactions"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') 
        folder = desktop + '\\__SBMLmodels__\\'
        """build folder on desktop if none exists"""
        if not os.path.exists(folder):
            os.makedirs(folder)
        path = folder + name + ".xml"
        """sbml file with model"""
        sbml_file = open(path,'w')
        sbml_file.write(sbml)
        sbml_file.close()

"""pickle file name"""
def pickle_store(path,data):
    with open(path, 'wb') as handle:
        pickle.dump(data, handle)
    return None

""" direct conversion to SBML"""
def directSBMLConversion(antimony_model,name):
        """if convergence then convert to SBML"""
        SBML = te.tellurium.antimonyToSBML(antimony_model)
        """store the sbml file"""
        path = storeSBML(SBML,name = name)    
        
def SBMLConversion(models,control):
    """convert the models in sbml format"""
    for number,model in models.items():
        """test conditions"""
        condition = control[number].conditions[-1]
        """transform the matrix"""
        fluxes,S,ratelaws,staterate = model_matrix_format(model.splitmodel,model.states,model.fixed)
        parsed_model = antimony_parser(model.states,model.parameters,fluxes,S,ratelaws,staterate,model.fixed,model.map) 
        """Load the lot"""
        antimony = te.loada(parsed_model)
        for ID,value in condition.items():
            if ID in model.states:
                ID = "y{}".format(str(model.map[ID]))
            if ID in model.parameters:
                ID = "k{}".format(str(model.map[ID]))            
            exec("antimony.{0} = {1}".format(ID,str(value)))
        """if convergence then convert to SBML"""
        SBML = te.tellurium.antimonyToSBML(parsed_model)
        """store the sbml file"""
        storeSBML(SBML,name = model.name)           
        
def convertModeltoSBML(model):
        """model matrix format which can be parsed into a antimony parser"""
        fluxes,S,ratelaws,staterate = model_matrix_format(model.splitmodel,model.states,model.fixed)
        """parsed model in atimony format which can be parsed into an SBML file"""
        antiparse = antimony_parser(model.states,model.parameters,fluxes,S,ratelaws,staterate,model.fixed,model.map)   
        """if convergence then convert to SBML"""
        SBML = te.tellurium.antimonyToSBML(antiparse)
        """store the sbml file"""
        storeSBML(SBML,name = model.name)      
    
    

        