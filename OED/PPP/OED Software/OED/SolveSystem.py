# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:39:08 2019

@author: huckg
"""

import re
import copy 
import numpy
import scipy.integrate  
import numpy
import matplotlib.pylab as plt

from pyDOE import *
from Distributions import *
from DataEvaluation import *


import libsbml
import importlib
import amici
import os
import sys
import numpy as numpy

class AmiciDataObject:
    def __init__(self,model,data,show = True):
        """rawdate ffrom the amici simulation"""
        self.rawdata = data
        """the time and rdata"""
        self.time    = data["ts"]
        self.rdata   = data["x"]
        self.ydata   = data["y"]        

        self.amiciStateIDs = model.Cmodel.getStateIds()
        """collect simdata by observable"""
        self.simdata = {}
        for i in range(len(self.amiciStateIDs)):
            self.simdata[model.invertedSBML[self.amiciStateIDs[i]]] = self.rdata[:,i]
        """"sbml parameter IDs"""
        parameterIDs = model.Cmodel.getParameterIds()
        """"final state of the jacobian"""
        self.jacobian = data["J"]
        """initial_conditions for the solver as the pattern changes"""
        self.initial_state = self.rdata[-1] 
        """extract the sensitivties"""
        self.forward_sensitivities = {model.invertedSBML[i]:{} for i in parameterIDs}
        
        """assess existence of FWD matrix"""
        fwd = data["sx"]
        
        try:
            if fwd == None:
                sensitivity = False
        except ValueError:
            sensitivity = True
            
        if sensitivity:
            self.initial_forward_sensitivity = fwd[-1]
            for i in range(len(parameterIDs)):
                for j in range(len(self.amiciStateIDs)):
                    parameter,state = model.invertedSBML[parameterIDs[i]],model.invertedSBML[self.amiciStateIDs[j]]
                    """add forward sensitivities to the state """
                    self.forward_sensitivities[parameter][state] = fwd[:,i,j]

               
    def TDI_update(self,model,data):
        self.time = numpy.append(self.time,data['ts']+self.time[-1])
        """add simulation to the vector of observable data"""
        for i in range(len(self.amiciStateIDs)):
            state = model.invertedSBML[self.amiciStateIDs[i]]
            """add simdata"""
            self.simdata[state] = numpy.append(self.simdata[state],data["x"][:,i])

        """"sbml parameter IDs"""
        amiciStateIDs = model.Cmodel.getStateIds()
        parameterIDs  = model.Cmodel.getParameterIds()
        """"final state of the jacobian"""
        self.jacobian = data["J"]
        """initial_conditions for the solver as the pattern changes"""
        self.initial_state = data['x'][-1] 
        
        """assess existence of FWD matrix"""
        fwd = data["sx"]
        try:
            if fwd == None:
                sensitivity = False
        except ValueError:
            sensitivity = True

        """update the sensitivities"""
        if sensitivity:
            self.initial_forward_sensitivity = fwd[-1]
            for i in range(len(parameterIDs)):
                for j in range(len(self.amiciStateIDs)):
                    parameter,state = model.invertedSBML[parameterIDs[i]],model.invertedSBML[amiciStateIDs[j]]
                    """add forward sensitivities to the state """
                    self.forward_sensitivities[parameter][state] = numpy.append(self.forward_sensitivities[parameter][state],fwd[:,i,j])           
            
    def sensitivity_update(self,measurement):
        """senstivity of the parameter"""
        sensitivity_factor = {}
        """distance from goal"""
        lsq_factor = {}
        """sensitivity to states"""
        for number,m in measurement.items():
            """observables and data"""
            for obs,data in m.profile.items():
                lsq_factor[obs] = numpy.array([(self.simdata[obs][i]-data[obs][i]) for i in range(len(data))])
            
        """sensitivity_factor"""
        observables = lsq_factor.keys()
        for ID,vector in self.forward_sensitivities.items():
            state,parameter = ID
            if state in observables:
                sensitivity_factor[(state,parameter)] = numpy.sum(numpy.transpose(lsq_factor[obs])*self.forward_sensitivities[(state,parameter)])

class PythonDataObject:
    def __init__(self,model,data,time):
        self.name  = model.name
        self.time  = time
        self.rdata = data
        """collect simdata by observable"""
        self.simdata = {}
        for i in range(len(self.rdata[-1])):
            self.simdata[model.states[i]] = data[:,i]

        """initial_conditions for the solver as the pattern changes"""
        self.initial_state = self.rdata[-1] 

    def testplot(self,):
        for k,v in self.simdata.items():
            plt.plot(self.time,v,label = k)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.legend(fancybox = True)
        folder  = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') + "\\TestPlot\\"
        """build folder on desktop if none exists"""
        if not os.path.exists(folder):
            os.makedirs(folder) 
        plt.savefig(folder + "testplot_{0}.png".format(self.name))
        plt.close()

class TelluriumDataObject:
    def __init__(self,model,data,time):
        self.name  = model.name
        self.time  = time
        self.rdata = data
        """collect simdata by observable"""
        self.simdata = {}
        
        """collect simdata by observable"""
        self.simdata = {}
        for i in model.states:
            ID ='['+ model.StringIDToSBML[i] +']'
            self.simdata[i] = self.rdata[ID]
        """initial state for a new simulation"""
        self.initial_state = {k:v[-1] for k,v in self.simdata.items()}
        
    def testplot(self,):
        for k,v in self.simdata.items():
            plt.plot(self.time,v,label = k)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.legend(fancybox = True)
        folder  = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') + "\\TestPlot\\"
        """build folder on desktop if none exists"""
        if not os.path.exists(folder):
            os.makedirs(folder) 
        plt.savefig(folder + "testplot_{0}.png".format(self.name))
        plt.close()
        
class ModelSolver:
    def __init__(self,model,
                 #parameters and pulse patterns
                 parameters   = [],
                 #set the initial conditions as a numpy array
                 variables = [],
                 #the forward problem true or false
                 forward_sensitivity = False,
                 adjoint_sensitivity = False,
                 measurement_time = [],
                 #the simulation times
                 simtime = (0,40),
                 dt = 1,
                 #tolerances
                 atol = 1e-12,
                 rtol = 1e-8,
                 #number of steps the solver makes
                 manual_TDI = {},
                 nsteps  =10**6,
                 pythonintegrator = "lsoda"):
        
        """assess if a pulsepattern is present and see if a manual overwrite is needed"""
        time_dependent_parameters = variables.time_dependent_parameters
        if len(manual_TDI) != 0:
            variables.amici_TDI = manual_TDI
            """update the boolean"""
            time_dependent_parameters = True

        """start and end time"""
        start,end = simtime
        """steps"""
        steps = int(end/dt)
        """set time"""
        time = numpy.linspace(start, end, steps)
        if len(measurement_time) != 0:
            time = measurement_time
        """set the timepoints"""
        """see if there is a compiled version of the model
        if not then use the python version to solve the equation"""
        try:
            Cmodel = model.Cmodel
            compilation = True
        except:
            compilation = False
    
        """the compilation is not present"""
        if not compilation:
            import tellurium as te
            try:
                equations = model.antiparse
                SBML   = True
                python = False
            except:
                SBML   = False
                python = True
                
            if SBML:
                equations = model.antiparse
                """M load models"""
                r = te.loada(equations)
                r.setIntegrator('cvode')
#               """set the parameters of the model from model variables"""
                for ID,value in variables.pID.items():
                    r.setValue(model.StringIDToSBML[ID], value)
                
                for state,value in variables.initial.items():
                    r.setValue(model.StringIDToSBML[ID],value)
                """simulate the object with tellurium which partially compiles the equations on the fly"""
                try:
                    results = r.simulate(start,end,steps)
                    time    = results['time']
                    self.data = TelluriumDataObject(model,results,time)
                    r.resetToOrigin()
                except RuntimeError:
                    python = True
                    print("Cvode was not able to solve the system, python's LSODA will be tried")
 
            python = False
            if python:
                equations = model.equations
                if forward_sensitivity:
                    equations = model.forward_sensitivities

                """solve the equation for what it is, i.e. no time dependent inputs"""
                def ODEfunc(t,y,p):
                    eq = []
                    for i in equations:
                        eq += [eval(i)]
                    return eq
                
                """ for solve the equations for altered states"""
                p = variables.p
            
                """set the timepoint"""
                t0 = start
                t  = end
                
                """"points and time points stored frm the integration"""
                points     = []
                timepoints = []
                
                """"set IC"""
                y = [i + 1*10**-1 for i in variables.ic]
    
                system = scipy.integrate.ode(ODEfunc).set_integrator(pythonintegrator,atol=atol,rtol=rtol,nsteps=nsteps,with_jacobian = False)
                system.set_f_params(p)
                """ Loop through the integration windows"""
                system.set_initial_value(y, t0)
                while system.successful() and system.t < t:
                    values = system.integrate(system.t + dt)
                    timepoints.append(system.t)
                    """update the pulse pattern"""
                    if time_dependent_parameters:
                        system.set_f_params(variables.python_TDI[int(t)])
                    """points value"""
                    points.append(values)   
                points = numpy.vstack(points) 
                t = numpy.array(timepoints) 
                """store in data object"""
                self.data = PythonDataObject(model,points,t)
            
        if compilation:
            Cmodel.setAllStatesNonNegative()
            Cmodel.setTimepoints(time)
            """get the solver"""
            solver = Cmodel.getSolver()   

            """"check if forward sensitivities are a thing"""
            if forward_sensitivity:
                solver.setSensitivityMethod(amici.SensitivityMethod_forward)
                solver.setSensitivityOrder(1)

            """set the solver tolerance"""
            solver.setAbsoluteTolerance(atol)
            solver.setRelativeTolerance(rtol)
            # solver.setMaxStepsBackwardProblem(10**6)
            solver.setMaxSteps(nsteps)
    
            """overwrite potentiall differing compiled values"""
            for ID,value in model.fixed.items():
                Cmodel.setParameterById(model.StringIDToSBML[ID],value)
    
            """set initial conditions"""
            indices = list(Cmodel.getStateIds())
            
            ic = [i + 0.0001 for i in variables.ic] 
            for k,v in variables.initial.items():
                    ic[indices.index(model.StringIDToSBML[k])] = v

            """if there is no pulse then we can """
            if not time_dependent_parameters:
                """dump minivalues in here to get the solver to work"""
                for ID,value in variables.pID.items():
                    Cmodel.setParameterById(model.StringIDToSBML[ID],value)
                """the Cmodel with the right initial states"""
                Cmodel.setInitialStates(ic)
                """run the simulation"""
                data = amici.runAmiciSimulation(Cmodel, solver)
                """"store in amici data object to be used in dat analysis"""
                self.data = AmiciDataObject(model,data)
                
            """if we give a pulse we need to stop start the solver"""
            if time_dependent_parameters:
                """get pulse length and find the interval for this pulse pattern
                get pulse length and define start to end"""
                times  = sorted(variables.amici_TDI.keys())
                """initial"""
                init = times[0]
                """number of steps"""
                steps = int(end/dt)
                """set the timepoints"""
                Cmodel.setTimepoints(numpy.linspace(start, end, steps))
                """set iniital conditions dump minivalues in here to get the solver to work"""
                ic = [i + 0.0001 for i in variables.ic] 
                for k,v in variables.initial.items():
                    ic[indices.index(model.StringIDToSBML[k])] = v

                """set first sequence of initial states"""
                Cmodel.setInitialStates(ic)
                """set parameters"""
                for ID,value in variables.pID.items():
                    Cmodel.setParameterById(model.StringIDToSBML[ID],value)
                for t in times:
                    vector = variables.amici_TDI[t]
                    """start and end time"""
                    start,end = t
                    """steps of data"""
                    steps = int((end-start)/dt)
                    """set the timepoints"""
                    Cmodel.setTimepoints(numpy.linspace(1,end-start, steps))
                    """set parameter vector with the new pulse scheme"""
                    for ID,value in vector.items():
                        Cmodel.setParameterById(model.StringIDToSBML[ID],float(value))
                    """ redefine the set of equations"""
                    data = amici.runAmiciSimulation(Cmodel, solver)
                    """append dat object"""
                    if t == init: 
                        self.data = AmiciDataObject(model,data)                        
                        Cmodel.setInitialStates(self.data.initial_state) 
                    else:
                        self.data.TDI_update(model,data)
                        Cmodel.setInitialStates(self.data.initial_state)  
                    """set the initial state of the system"""
                    
                
    def __getData__(self,):
        return self.data
    
def create_flux_vectors(n_sub):
    flux_vector = defaultdict(list)
    inv_flux_vector = defaultdict(list)    
    for i,re in enumerate(n_sub.re):
        reactants,poducts = re
        a,b = reactants
        if a != None:
            flux_vector[i].append(a)
        if b != None:
            flux_vector[i].append(b)
    for key,values in flux_vector.items():
        for value in values:
            inv_flux_vector[value].append(values)
    return flux_vector,inv_flux_vector

def gillespie(state,n_sub, max_t):
    flux_vector,rev_flux_vector = create_flux_vectors(n_sub)
    flux_constants = n_sub.k
    S = n_sub.S
    
    t = 0
    state_history = [np.array(state)]
    time_history = [t]

    p = np.array(flux_constants)
    for i, states in flux_vector.items():
        for s in states:
            p[i] *= state[s]

    i = -1
    while t < max_t:
        p_cumsum = np.cumsum(p)
        p_sum = p_cumsum[-1]
        p_cumsum /= p_sum

        # Pick reaction
        rand = random.random()
        r = np.searchsorted(p_cumsum, rand)

        # Update
        dstate = S[:, r]
        state += dstate

        # Update dirty propensities
        for si in np.where(dstate)[0]:
            for pi in rev_flux_vector[si]:
                p[pi] = flux_constants[pi]
                for s in flux_vector[pi]:
                    p[pi] *= state[s]
        t += 1.0 / p_sum * np.log(1.0 / rand)

        state_history.append(np.array(state))
        time_history.append(t)
    return np.array(state_history), np.array(time_history)

def stochastic_solver(n_sub):     
    state = np.ones(len(n_sub.S))
    s, t = gillespie(state,reactions,100)
    s_i = 20
    tp = np.arange(0.0, 100.0, 1.0 / 60.0)
    a = np.interp(tp, t, s[:, s_i])
    return tp,a
    


