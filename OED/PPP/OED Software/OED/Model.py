# -*- coding: utf-8 -*-
"""
Created on Thu Feb 07 12:48:40 2019

@author: huckg
"""
#

import libsbml
import importlib
import amici
import os
import sys
import numpy as numpy
from Parse import *
#import tellurium as te

from importlib import reload
import datetime
import sympy
import pandas as pd
import copy
from time import sleep
from DataTransform import * 
from Operations import *

def runModel(model = '',models = {},
             #which subfunctions should be on
             derivatives = True,
             matrix      = True,    
             compilation = True,
             #for the sensitvities and derivatives
             forward_observables    = [],
             sensitivity_parameters = [],
             #show the reaction generated by antimony and C++
             show_reactions = False):
    
    if model:        
        """äctivate the following subfunctions for the model"""
        if derivatives:
            model.Derivatives(forward_observables = forward_observables,sensitivity_parameters=sensitivity_parameters)
        if compilation:
            model.SBMLconversion()
            model.PytoCompile(show_reactions = show_reactions)
        """return the model to benchmark or your module module"""
        return model

    if models:
        activated_models = {}
        """loop through the library of models"""
        for number,model in models.items():
            """äctivate the following subfunctions for the model"""
            if derivatives:
                model.Derivatives(forward_observables = forward_observables,sensitivity_parameters=sensitivity_parameters)
            if compilation:
                model.SBMLconversion()
                model.PytoCompile(show_reactions = show_reactions)
            """return the model to dict"""
            activated_models[number] = model
        return activated_models

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

class ModelVariables:
    def __init__(self,model,conditions = {},modification = {}):
        self.model = model
        """set time dependent parameters"""
        self.time_dependent_parameters = False
        """assess how the conditions relate to the parameters"""
        """"set intitial of the conditions of the """
        self.initial,control = {},{}
        for i,j in conditions.items():
            if i not in model.states and i not in model.boundaries.keys():
                print("There is a namespace error, some {} string in the measurement and or conditions does not match any defined in the model".format(i))
            """see if it is a state or a parameter"""
            if i in model.states:
                self.initial[i] = j
            elif i in model.boundaries.keys():
                control[i] = j
                
        """assess whether control parameters that act as initial conditions actually exist"""
        self.StateAndControl = []
        if model.initial_control:
            parameter,states = zip(*model.initial_control) #unzip parameters and states
            boolean = [i for i in states if i in model.states] #check if they exist at all in the current model
            self.StateAndControl  = [model.initial_control[states.index(i)] for i in boolean] #include them in initial control if they exist
        
        """ID value vector"""
        pID = {}
        """parameter vector"""
        p_vector = numpy.zeros(len(model.fixed))   
        
        """set the fixed values of the parameters"""
        for ID,value in model.fixed.items():
            """set vector"""
            p_vector[model.map[ID]] = value 
            """set ID"""
            pID[ID] = value
            
        """set the sampled parameter values of the parameters"""
        for ID,value in modification.items():
            """set vector"""
            p_vector[model.map[ID]] = value 
            """set ID"""
            pID[ID] = value
                    
        """set the control values of the parameters"""
        for ID,value in control.items():
            """set vector"""
            p_vector[model.map[ID]] = value 
            """set ID"""
            pID[ID] = value

        """"the parameter vector and dictionary"""        
        self.p = p_vector
        self.pID = pID       
        
        """note this does wqork if your control and initial conditions are the same but 
        seperated such that control supercedes initial conditions, now it is possible to
        set initial conditoins but seperate control variables"""
        self.ic = numpy.zeros((len(model.states)))
        
        """add the initial conditions reflecting the point if initial control"""
        for parameter,state in self.StateAndControl:
            self.ic[model.map[state]] = self.pID[parameter]
            
        """add the right values to the states"""
        for state,value in self.initial.items():
            self.ic[model.map[state]] = value 

        """python TDI (time dependendent input)"""
        self.python_TDI = {}
        """amici TDI"""    
        self.amici_TDI  = {}
            
    def set_time_dependent_variables(self,time_dependent_variables,space):
        self.time_dependent_parameters = True
        for time in list(time_dependent_variables.values())[-1]:
            """"create the amici and python tdi input vectors"""
            dct = {}
            """loop"""
            for parameter in time_dependent_variables.keys():
                dct[parameter] = space[parameter][time_dependent_variables[parameter][time]]
            self.amici_TDI[time] = dct
        
        """assess the model and check if a compiled monniker exists"""
        try:
            Cmodel =self. model.Cmodel
            compilation = True
        except:
            compilation = False            
        """the compilation is not true create pythonvectors"""
        if compilation == False:
            print("For this to run effectively you need to compile the code, otherwise you'll be simulating for ages")
        
    def update(self,varlist):
        sc = dict(zip(*self.StateAndControl))
        """update the parameter dict dynamically i.e. in case of a local sensitivity change"""
        for i in varlist:
            for ID,value in i.items():
                if ID in self.model.states:
                    self.ic[self.model.map[ID]] = value
                else:
                    self.p[self.model.map[ID]] = value
                    self.pID[ID] = value
                if ID in sc.keys():
                    self.ic[self.model.map[sc[ID]]] = value
                
class ModelObject:
    def __init__(self,model,states,boundaries,fixed,
                 experimental_conditions   = [],
                 initial_conditions        = [],
                 control_parameters        = [],
                 diffusion_parameters      = [],
                 initial_control           = [],
                 likelihood                = [],
                 ratetypes                 = [],
                 sensitivity_parameters    = [],
                 
        
                 include                   = [],
                 plot                      = [],
                 observables               = [],
                 lumped                    = {},
                 spacing = 1000,
                 name = ''):
 
        """name of the model"""
        self.name        = name          #the name of the model used for storage
        """"bool for the compilation of a model to C, if the compilation subfunction has been activated it will be True"""
        self.compilation = False 
        """splot the model into sections"""
        self.splitmodel = re.split('; |,| \n ',model)
        """replace the delimeters for the fluxes"""
        model = model.replace("|+|","+")
        model = model.replace("|-|","-")
        
        """These are the 4 user defined features of a model
        including a model in stringformat, its correspnding states indexed
        according to the model equations, the measured observables
        and the boundaries""" 
        self.stringmodel = model         #model in stringformat
        self.states      = states        #states in the model in order
        self.observables = observables   #observables we measure
        self.lumped      = lumped        #grouped observables/states we measure
        self.boundaries  = boundaries    #boundaries for rates in the models
        self.control     = control_parameters
        self.spacing     = spacing       #spacing of the parameter vector
        self.diffusion   = diffusion_parameters
        
        """the sensitivity parameters"""
        self.sensitivity_parameters = sensitivity_parameters
        if len(self.sensitivity_parameters) == 0:
            self.sensitivity_parameters = [i for i in fixed.keys() if i not in self.control]
               
        """Build the parameter space needed to obtain a space that can be sampled for simulation"""
        self.u_space,exclude = define_parameter_boundaries(self.boundaries,spacing = self.spacing) 
        self.fixed  = fixed
        
        """this includes information to be included in possible model analyses"""
        self.include = include       #include these parameters in the analyses
        if not self.include:
            self.include = self.boundaries.keys()
        """ïf the boundaries of a parameter set are (0,0) they are not taken into account in the analysis"""
        self.include = [i for i in self.include if i not in exclude]
        self.plot = plot
        if not self.plot:
            self.plot = self.include     #include these parameters in the plotting
        """set space for priors"""
        self.priors = []
        """set the pulsespace"""
        self.pulsespace = {i:self.u_space[i] for i in self.control}
        """these are the basic units of a model needed to run it
        this includes a set of equations taht are compiled to within the class
        to the model object, we include a map for the observables and 
        the parameter vector to ensure we have a viable system"""
        self.strq = build_equations(self.stringmodel,self.states)


        """"parameter vector...............................................
        SORTING THE PARAMETER VALUES IS IMPERATIVE WHEN INDEXING< THE SBML FILES
        ARE GENERATED SEPERATELY IN PYTHON 2.7 WHERE INDICES WILL BE GIVEN TO THE
        SPECIES AND PARAMETERS, IF THE DICT CALL SCRAMBLES THEM DIFFERENTLY SIMULATIONS
        WILL NOT MATCH AND PARAMETER ASSIGNMENTS WILL BE WRONG, ITS "DANGEROUS" BUT
        MARKED UNDER WILL MAYBE FIX LATER"""
        keys = sorted(list(self.boundaries.keys()))
        vector = {keys[i]:"p[{}]".format(i) for i in range(len(keys))}                                    
        self.equations,self.uncompiled_equations = build_base_model(self.strq,vector)   
        print(self.uncompiled_equations)
        
        '''complete map of all the indices '''
        self.map = {}
        for i in vector.keys():
            self.map[i] = strip_string(vector[i])
        for i in self.states:
            self.map[i] = self.states.index(i)
        self.obs = [(i,self.map[i]) for i in self.observables]
        
        """append indices of lumped observables in list"""
        self.lumped = lumped #grouped observables/states we measure
        for k,v in self.lumped.items():
            n = []
            for i in v:
                n.append(self.map[i])
            self.lumped[k] = n

        """vector with initial conditions"""
        self.initial_conditions = initial_conditions
        """The initial conditions need to be set according to an input parameter"""
        self.initial_control    = initial_control
        """complete map of rate types"""
        self.ratetypes  = ratetypes  
        """define the maximum likelihood"""
        self.maximum_likelihood = {}
        """"reset the entire space by including likelihoods, uniform if nothing is given
        the sampled likelihood if there is an actual likelihood"""
        self.p_space            = self.u_space
        self.maximum_likelihood = self.fixed
        
        """add prior likelihoods to the model parameters"""
        self.likelihood = likelihood
        self.priors.append((self.likelihood,self.p_space))
        """sample the likelihood"""
        if self.likelihood:
            self.p_space = sample_likelihood(self.likelihood,spacing = self.spacing)
        if self.likelihood:
            """"define maximum likelihood"""
            for parameter,probabilities in self.likelihood.items():
                if probabilities:
                    p,xvals = probabilities
                self.maximum_likelihood[parameter] = xvals[list(p).index(max(list(p)))]
        """add differentials 1st, second, parameter, and sensitvities"""
        """if the parameter is assigned in a dict but left out of the analysis
        the map will assign indices wrongly to the differentiation vector
        ergo there will be 27 parameters/indices but 26 diff objects meaning 
        p27 will be overwritten by p2 should there be 27-26 parameters"""
        self.parameters = sorted(list(set(list(self.boundaries.keys()) + list(self.fixed.keys()))))

        """In order to COMPILE the model its needs to be converted to SBML, to go from
         a textbased model to SBML will require a Model in matrix format"""
        self.staterate = {}
        """loop through the equations and find rate laws denoted by + or -, |+| means internal to rate"""
        for equation in range(len(self.splitmodel)):
            eq = self.splitmodel[equation]
            y  = self.states[equation]
            """define ratelaws"""
            ratelaws = define_ratelaw(eq)  
            """class of formalized rate laws"""    
            for rate in ratelaws:
                r,statelist = [],[]
                """"find the sign"""
                sign = rate[0]
                """"loop through parameters and append if in rate"""
                for i in self.fixed.keys():
                    if i in rate:
                        r.append(i)
                """"loop through states and append if in rate"""                
                for i in self.states:
                    if i in rate:
                        statelist.append(i)
                if y not in self.staterate.keys():
                    self.staterate[y] = []
                self.staterate[y].append((sign,rate[1:len(rate)]))
        
        
        """fluxvector and the stoichiometry"""
        self.fluxes = list(set([item.strip() for sublist in self.staterate.values() for sign,item in sublist]))
        """create stoichiometric and flux vector this can subsequently be used as a
        tool for creating the nessecary tellurium models"""
        self.S = numpy.zeros((len(self.states),len(self.fluxes)))
        """fill up the vector with fluxes"""
        for i in range(len(self.states)):
            """unzip the fluxes"""
            signs,ratelaws = zip(*self.staterate[self.states[i]])
            for j in range(len(ratelaws)):
                """set the flux"""
                flux = ratelaws[j].strip()
                """ratelaws added to stoichiometry"""
                for n in range(len(self.fluxes)):
                    if flux == self.fluxes[n].strip():
                        if signs[j].strip() == "-":
                            self.S[i][n] = -1
                        elif signs[j].strip() == "+":
                            self.S[i][n] = 1    
        """update ratelaws """
        self.ratelaws = {i:RateLaw(i,self.states,self.fixed) for i in self.fluxes}     
#        *__________________Directories where model results can be stored_____________
        """this function builds a map and stores the results of your model analysis automatically
        this includes a folder with all the figures, per model,  and the classses"""
        desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') 
        figures = '\\Figures\\'
        classes = '\\Classes\\'
        directory = desktop + '\\' + "__Models__"  + '\\' 
        modelfolder = directory + name
        if not os.path.exists(directory):
            os.makedirs(directory)
        if not os.path.exists(modelfolder):
            os.makedirs(modelfolder) 
        for i in [figures,classes]:
            if not os.path.exists(modelfolder+i):
                os.makedirs(modelfolder+i)        
        self.figure_folder = modelfolder + figures
        self.class_folder  = modelfolder + classes
        self.figure_path   = self.figure_folder + self.name
        self.class_path    = self.class_folder + self.name

    def Derivatives(self,forward_observables = [],sensitivity_parameters = []):
        ratemap   = {}
        """equations"""
        equations = copy.copy(self.stringmodel)
        """"ratemap"""
        terms = {i.strip():j for i,j in self.map.items()}
        for i in self.parameters:
            ratemap[i] = " p{} ".format(self.map[i.strip()])
        iratemap = {v:k for k,v in ratemap.items()}
        
        """create a state map with a p1 to p_n term a sympy exec function can compile
        this way cannonicalization is a less "bad" coding practice""" 
        for k,v in ratemap.items():
            k = k.strip() 
            k = " "+ k +" "
            equations = equations.replace(k,v) 
            
        """define the symbols for the parameters"""  
        psym = [i for i in ratemap.values()]
        for i in psym:
            i = sympy.symbols(i)
    
        """create a state map with a y1 to yn a sympy exec function can compile
        this way cannonicalization is a less "bad" coding practice"""   
        statemap   = {}
        """"statemap"""
        for i in self.states:
            statemap[i] = " y{} ".format(terms[i.strip()])
    
        """loop through the equation"""
        for k,v in statemap.items():
            k = k.strip()
            k = " "+ k +" "
            equations = equations.replace(k,v)
    
        """its either defined in different lines or seperated by comma's"""
        if '\n' in equations:
            equations = equations.split("\n")
        if ',' in equations:
            equations = equations.split(',')
     
        sym = [statemap[i] for i in self.states]
        """execute the definition of the sympy sympols ALWAYS y1 to yn"""
        for i in sym:
           i = sympy.symbols(i.strip())
        """empty hessian and jacobian"""
        self.jacobian = []; self.hessian = []
        """loop through equations and differentiate, loop through states and
        differentiate, for the hessian loop through the states twice"""
    
        for equation in equations:
            """jacobi row"""
            jacobian_row = []
            """hessian function row"""
            hessian_function = []
            for si in sym:
                """differentiate dfdx i.e. the function with respect to state x"""
                dfdx = sympy.diff(equation,si)
                jacobian_row.append(str(dfdx))
                
                """"define hessian row"""
                hessian_row = []
                for sj in sym:
                    """differentiate dfdx again with respect to state y"""
                    dfdxdy = sympy.diff(dfdx,sj)
                    hessian_row.append(str(dfdxdy))
                for i in range(len(hessian_row)):
                    for k,v in statemap.items():
                        hessian_row[i].replace(v,k)     
                """update hessian"""
                hessian_function.append(hessian_row)
               
            """revert back to original format"""
            for i in range(len(jacobian_row)):
                for k,v in statemap.items():
                    jacobian_row[i].replace(v,k)
                    
            """update jacobian"""
            self.jacobian.append(jacobian_row)
            self.hessian.append(hessian_function)
    
        """the states order is fixed the parameters is not, therefore 
        we need to return the parameter derivates wit a map signifiing what
        the jth column was differentiated towards"""
        pd = {i:[] for i in ratemap.values()}
        for equation in equations:
            for pi in psym:
                dfdp = sympy.diff(equation,pi)
                pd[str(pi)].append(str(dfdp))
                
        """reformat the equations to suit an index, we reverse the list so p[1] does not replace p[11] e.g. p1 subbed in p1...1
        this needs to be done because sympy does not take/differiatate indexed vector calls"""
        vector = {"y{}".format(i):"y[{}]".format(i) for i in range(len(self.states))}
        vector.update({"p{}".format(i):"p[{}]".format(i) for i in range(len(self.parameters))})
        """reversed elements in the matrices"""
        elements = reversed(sorted(vector.keys()))
        """"loop through the matrices and confirm the differentials"""
        for i in elements:
            """alter the jacobian"""
            for function in range(len(self.jacobian)):
                for state in range(len(self.jacobian[function])):
                    self.jacobian[function][state] = self.jacobian[function][state].replace(i,vector[i])
            """alter parameter derivatives"""
            for p,fxdp in pd.items():
                for term in range(len(fxdp)):
                    pd[p][term] = pd[p][term].replace(i,vector[i])
                    
            """alter hessian to include indices"""
            for function in range(len(self.hessian)):
                for row in range(len(self.hessian[function])):
                    for column in range(len(self.hessian[function][row])):
                        self.hessian[function][row][column] = self.hessian[function][row][column].replace(i,vector[i])        
    
        self.c_jacobian = []
        for i in range(len(self.jacobian)):
            line = []
            for j in range(len(self.jacobian[-1])):
                line.append(compile(self.jacobian[i][j],'jacobian','eval'))
            self.c_jacobian.append(line)
    
        """derivatives of the system"""
        self.derivatives = {}
        for k,v in pd.items():
            self.derivatives[iratemap[k]] = v

        """define a set of forward sensitivity equations and obtain"""
        self.sensitivity_equations,self.sensitivity_map,self.sensitivity_states = {},{},[]

        observables = self.states
        """sensitivity equations"""
        if not sensitivity_parameters:
            self.sensitivity_parameters = [i for i in list(set(list(self.boundaries.keys())+list(self.fixed.keys())))]
        """a set of sensitivity states added to the state vector"""
        index = {self.states[i]:i for i in range(len(self.states))}
        for i in observables:
            for j in self.sensitivity_parameters:
                self.sensitivity_states.append((i,j))
                
        """add parameter derivatives to the equations"""        
        for i in range(len(observables)):         
            for j in range(len(self.sensitivity_parameters)):
                """add term to equation"""
                if self.derivatives[self.sensitivity_parameters[j]][index[observables[i]]] != str(0):
                    self.sensitivity_equations[(observables[i],self.sensitivity_parameters[j])] = self.derivatives[self.sensitivity_parameters[j]][index[observables[i]]]
                    
        """all combinations of states and parameters"""
        combinations = []
        """add state derivatives to the equations"""
        for i in range(len(self.sensitivity_parameters)):
            """parameter term in the sensitivity vector""" 
            parameter = self.sensitivity_parameters[i]
            for function in range(len(observables)):
                fncobs = observables[function]
                """build equations combining vectors"""
                equation = ''
                for j in range(len(observables)):
                    """states of the network"""
                    state = observables[j]
                    """derivative of the jacobian"""
                    j_derivative = self.jacobian[index[observables[function]]][j]
                    """assemble equation with derivative"""
                    if j_derivative != str(0):
                        equation +=  "+ ( " + j_derivative + " ) " + " * {} ".format(state+"%&"+parameter)
                        combinations.append(state+"%&"+parameter)
                    if (state,parameter) not in self.sensitivity_equations.keys():
                        self.sensitivity_equations[(fncobs,parameter)] = ''
                self.sensitivity_equations[(fncobs,parameter)] += equation
                
        """filter the sensitivity equations all null terms"""        
        self.sensitivity_equations = {k:v for k,v in self.sensitivity_equations.items() if v != "0"}
        
        count = 0
        """update the map if the state is bigger than zero"""
        for i in self.sensitivity_states:
            self.map[i] = len(self.states) + count    
            count += 1
        
        """loop throught the equations and add a y-vector"""
        for state,parameter in self.sensitivity_equations.keys():
            for k,v in self.sensitivity_equations.items():
                self.sensitivity_equations[k] = v.replace(" "+ state+"%&"+parameter +" " ,"y[{}]".format(self.map[(state,parameter)]))

        """if there is no derivative add zero term"""
        for i in self.sensitivity_states:
            if not self.sensitivity_equations[i]:
                self.sensitivity_equations[i] = '0'
        
        from math import log
        self.theoretical_FSE = {}
        """append list of equations of fully observed compiled sensitivity equations"""
        for k,v in self.sensitivity_equations.items():
            self.theoretical_FSE[k] = compile(v.strip(),'sensitivity_equations','eval')
        self.observed_FSE = self.equations + [self.theoretical_FSE[i] for i in self.sensitivity_states]
        
    def AntimonyConversion(self,show_reactions = False):
        """Create antimony version of system normally reserved for SBML conversion but simulatable with tellurium"""
        self.antiparse = antimony_parser(self.states,self.parameters,self.fluxes,self.S,self.ratelaws,self.staterate,self.fixed,self.map)
        """get proper ID's"""
        self.StringIDToSBML = {}
        for k,v in self.map.items():
            if type(k) != tuple:
                if k in self.parameters:
                    self.StringIDToSBML[k] = "k{}".format(self.map[k])
                elif k in self.states:
                    self.StringIDToSBML[k] = "y{}".format(self.map[k])
        self.invertedSBML = {}
        for k,v in list(self.StringIDToSBML.items()):
            self.invertedSBML[v] = k
         
    def SBMLconversion(self,show_reactions = False):
        self.compilation = True
        self.SBMLfolder  = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') + "\\__SBMLmodels__\\"
        """build folder on desktop if none exists"""
        if not os.path.exists(self.SBMLfolder):
            os.makedirs(self.SBMLfolder)        
        """get a list of the sbml models currently stored in the SBML folder
        if no models are present will create a SBML formatted model"""
        SBMLmodels =  [f for f in os.listdir(self.SBMLfolder) if os.path.isfile(os.path.join(self.SBMLfolder, f))]
        """"convert this model to matrix format"""
        self.filecheck_sbml = False
        print(self.name)
        if self.name + ".xml" not in SBMLmodels:
            import tellurium as te
            print(""""The SBML version of this model does not exist because this library
                  does not take text thus create this first""")
            """parsed model in atimony format which can be parsed into an SBML file"""
            self.antiparse = antimony_parser(self.states,self.parameters,self.fluxes,self.S,self.ratelaws,self.staterate,self.fixed,self.map)
            """if convergence then convert to SBML"""
            self.SBML = te.tellurium.antimonyToSBML(self.antiparse)
            """store the sbml file by starting to set set the path of the reaction list"""
            desktop = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') 
            folder = desktop + '\\__SBMLmodels__\\'
            path = folder + self.name + ".xml"
            """sbml file with model"""
            sbml_file = open(path,'w')
            sbml_file.write(self.SBML)
            sbml_file.close()
            """boolean which signals file exists"""
            print("""the script will be rerun and reloaded so the libraries are in the proper 
                  order this will keep happenening untill all model files have been generated""")
            os.execv(sys.executable, [sys.executable] + sys.argv)
            
        """model location"""
        self.locSBMLmodel = self.SBMLfolder+ "{}.xml".format(self.name)
        """the parameters should all be in line with the index map generated, in this manner 
        we can which parameters are which and which states are which"""
        self.StringIDToSBML = {}
        for k,v in self.map.items():
            if type(k) != tuple:
                if k in self.parameters:
                    self.StringIDToSBML[k] = "k{}".format(self.map[k])
                elif k in self.states:
                    self.StringIDToSBML[k] = "y{}".format(self.map[k])
        self.invertedSBML = {}
        for k,v in list(self.StringIDToSBML.items()):
            self.invertedSBML[v] = k
        
    """in this function we convert the models using tellurium and SBML 
    to create standardized files to be used in AMICI and Pesto"""
    def PytoCompile(self,show_reactions = False):
        self.compilation = True
        """check if there is a module
        import the model by extending the
        path in the python environment"""
        self.AMICImodelfolder  = os.path.join(os.path.join(os.path.expanduser('~')), 'Desktop') + "\\__AMICImodels__\\"
        if not os.path.exists(self.AMICImodelfolder):
            os.makedirs(self.AMICImodelfolder)
        set_dir = self.AMICImodelfolder  + "{}\\".format(self.name) 
        sys.path.insert(0, set_dir)
        """assess if we can import the PytoC++ module"""
        try:
            modelModule = importlib.import_module(self.name)
            CplusPy = True
        except ModuleNotFoundError:
            CplusPy = False

        """check if the compiled version of the model already exists in the working dir
        if so then do not compile this system it will take to long""" 
        if not CplusPy:
            """set sbml model name"""
            sbmlname = self.locSBMLmodel
            print(self.locSBMLmodel)
            """get the sbml document from folder"""
            document = libsbml.readSBML(self.locSBMLmodel)
            """the sbml file"""
            model = document.getModel()
            """show model reactions"""
            if show_reactions:
                print('\nReactions:')
                for reaction in model.getListOfReactions():
                    reactants = ' + '.join(['%s %s'%(int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '', r.getSpecies()) for r in reaction.getListOfReactants()])
                    products  = ' + '.join(['%s %s'%(int(r.getStoichiometry()) if r.getStoichiometry() > 1 else '', r.getSpecies()) for r in reaction.getListOfProducts()])
                    reversible = '<' if reaction.getReversible() else ''
                    print('%3s: %10s %1s->%10s\t\t[%s]' % (reaction.getId(),
                                        reactants,
                                        reversible,
                                        products,
                                        libsbml.formulaToL3String(reaction.getKineticLaw().getMath())))
            
            """the sbml importer """
            sbml_importer = amici.SbmlImporter(sbmlname)
            """compile the model with the sbml to amici step"""
            sbml_importer.sbml2amici(self.name,
                                     output_dir = set_dir,
                                     verbose=True)
            """import the model by extending the path in the python environment"""                
            set_dir = self.AMICImodelfolder + "{}\\".format(self.name) 
            sys.path.insert(0, set_dir)
            try:
                print(self.AMICImodelfolder + "{}".format(self.name))
                modelModule = importlib.import_module(self.AMICImodelfolder + "{}".format(self.name))
               
            except ModuleNotFoundError:
                print("This module has been created but has not been loaded, rerun the script")
            """rerun the script so modules can be loaded properly"""
            print("the script will be rerun and reloaded so the libraries are in the proper order")
            os.execv(sys.executable, [sys.executable] + sys.argv)
                
        try:
            """the C++ python module of the model"""
            self.Cmodel = modelModule.getModel()
            """set the states to be nonnegative we cannot get negative concentrations"""
            self.Cmodel.setAllStatesNonNegative() 
            """create a map for the states and parameters in the sbml file and the human readable parameters"""
            self.AmiciParameterIDs  = self.Cmodel.getParameterIds()
            self.AmiciObservableIDs = self.Cmodel.getObservableIds()
            self.AmiciStateIDs      = self.Cmodel.getStateIds()
        except UnboundLocalError:
            print("There is a local unbound error downloading them model")

        