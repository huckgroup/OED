"""
Created on Thu Dec 13 13:37:30 2018

@author: huckg
"""


import re
import copy 
import numpy
import scipy.integrate  
import itertools
import os
import pickle
import math
import sobol_seq
import random
import matplotlib.pylab as plt
import sympy
from itertools import tee
from pyDOE import *
from Distributions import *
from DataEvaluation import *


def derivatives(states,parameters,strq,terms):
    ratemap   = {}
    """equations"""
    equations = copy.copy(strq)
    """"ratemap"""
    terms = {i.strip():j for i,j in terms.items()}
    for i in parameters:
        ratemap[i] = " p{} ".format(terms[i.strip()])
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
    for i in states:
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
        
 
    sym = [statemap[i] for i in states]
    """execute the definition of the sympy sympols ALWAYS y1 to yn"""
    for i in sym:
       i = sympy.symbols(i.strip())
    """empty hessian and jacobian"""
    jacobian = []; hessian = []
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
        jacobian.append(jacobian_row)
        hessian.append(hessian_function)

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
    vector = {"y{}".format(i):"y[{}]".format(i) for i in range(len(states))}
    vector.update({"p{}".format(i):"p[{}]".format(i) for i in range(len(parameters))})
    """reversed elements in the matrices"""
    elements = reversed(sorted(vector.keys()))
    """"loop through the matrices and confirm the differentials"""
    for i in elements:
        """alter the jacobian"""
        for function in range(len(jacobian)):
            for state in range(len(jacobian[function])):
                jacobian[function][state] = jacobian[function][state].replace(i,vector[i])
        """alter parameter derivatives"""
        for p,fxdp in pd.items():
            for term in range(len(fxdp)):
                pd[p][term] = pd[p][term].replace(i,vector[i])
                
        """alter hessian to include indices"""
        for function in range(len(hessian)):
            for row in range(len(hessian[function])):
                for column in range(len(hessian[function][row])):
                    hessian[function][row][column] = hessian[function][row][column].replace(i,vector[i])        

    compiled_jacobian = []
    for i in range(len(jacobian)):
        line = []
        for j in range(len(jacobian[-1])):
            line.append(compile(jacobian[i][j],'jacobian','eval'))
        compiled_jacobian.append(line)

    """derivatives of the system"""
    derivatives = {}
    for k,v in pd.items():
        derivatives[iratemap[k]] = v
    return compiled_jacobian,jacobian,hessian,derivatives

def clustering(distance):
    from scipy.spatial import distance as sqr
    from scipy.cluster.hierarchy import linkage 
    y = sqr.squareform(distance)
    clstr = linkage(y, method='single', metric='euclidean')
    import scipy
    scipy.cluster.hierarchy.dendrogram(clstr)
    return clstr

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

def get_fixed_conditions(model):
    return {i:model.fixed[i] for i in model.control}

class collinearity_object:
    def __init__(self,parameters,collinearity,observables):
        self.parameters   = parameters
        self.collinearity = collinearity
        self.observables  = observables 

def build_constructs(promoters,ORF,tags = []):
    constructs = []
    for i in promoters:
        for ii in ORF:
            if tags:
                for iii in tags:
                    constructs.append(i+"."+ii+":"+iii)
            constructs.append(i +'.'+ ii)
    return constructs 

def inputvector(pset,inputs,boundaries,time):
    if not inputs:
        return []
    prf = {i:[] for i in range(time+(time*3))}
    parameters = [(p,tstart,tend) for p,value,tstart,tend in inputs]
    for i in range(len(prf)):
        for p,start,end in parameters:
            if start <= i >= end:
                prf[i].append((p,pset[p]))
    for p,value,tstart,tend in inputs:
        for i in range(tstart,tend,1):
            """list of values that can be used to build prf"""
            if value == "Max":  
                prf[i].append((p,boundaries[p][-1]))
            if value == "Min":  
                prf[i].append((p,boundaries[p][0]))
            if type(value) == int or type(value) == float:
                prf[i].append((p,value)) 
            try:
                if "Factor" in value:
                    value = float(value.split(" ")[-1])
                    prf[i].append((p,current * value))
            except:
                pass
    return prf

def time_units():
    nmc = ["min","hour","hrs","sec","hours","mins"]
    return nmc

def measured_states():
    nmc = ["observable"]

def unitfinder(IDs):  
    lower_case_IDs = [i.lower() for i in IDs]
    standard_case_IDs = IDs
    for i in time_units():
        tracker = 0
        for j in lower_case_IDs:
            if i == j:
                return standard_case_IDs[tracker]
            tracker += 1
    return "There is not unit check your datafile"


def stringfinder(IDs,base):
    lower_case_IDs = [i.lower() for i in IDs]
    standard_case_IDs = IDs;identifier = base; observable_ID = []
    tracker = 0
    for i in lower_case_IDs:
        if identifier in i:
             observable_ID.append(standard_case_IDs[tracker])
        tracker += 1
    return observable_ID   
        
def getcolumn(df,header):
    return str(df[header].dropna()[0])

def dfix(triset,biset):
    p = [i for i in triset if i not in biset]
    return p[0]

def build_array(indexinfo,dimensions,data):
    arrays  = {i:[] for i in dimensions}
    for i in indexinfo:
        for j in arrays.keys():
            arrays[j].append(data[j][i])
    return tuple([numpy.array(arrays[i]) for i in dimensions])


def set_ic(states,ic,conditions,p_map,p):#states,initial conditions mapped onto parameters,parameter values
    if len(ic) != len(states):
        y0 = numpy.zeros((len(states))) 

    if conditions:
        si = {states[i]:i for i in range(len(states))}
        y0 = numpy.zeros((len(states))) + numpy.array(ic)
        for parameter,state in conditions:
            y0[si[state]] = p[p_map[parameter]]
        return y0
    else:
        return numpy.array(ic)


def count_letters(molecule, searched = 'C'):
    counter = 0
    for letter in molecule:
        if letter == searched:
            counter += 1 
    return counter 

#create novel reaction class
def create_novel_reaction_class(rt,boundaries = (0.001,1000),spacing = 1000):
    from NetworkSpace import ParameterClass
    """add the probability for the non-existed space"""
    u_space,exclude = define_parameter_boundaries({rt:boundaries},spacing = spacing)
    """"the probability generated the """
    probability = (numpy.array([1/float(spacing) for i in range(len(u_space[rt].values()))]),list(sorted(u_space[rt].values())))
    """"the smarts or "parameter class" in the network"""
    return ParameterClass(rt,rt,boundaries,probability)


#regression statistics for datasets
def reg_m(y, x):
    import statsmodels.api as sm
    import numpy as np
    ones = np.ones(len(x[0]))
    X = sm.add_constant(np.column_stack((x[0], ones)))
    for ele in x[1:]:
        X = sm.add_constant(np.column_stack((ele, X)))
    results = sm.OLS(y, X).fit()
    return results

def ParameterInclusion(model):
    include = []
    """check the boundaries to auto include/exlude
    parameters from the optimization to improve convergence times"""
    boundaries = model.boundaries
    for k,v in boundaries:
        minimum,maximum = v
        if minimum > maximum:
            factor = minimum/maximum
        if maximum > minimum:
            factor = maximum/minimum
        if factor > 1.1:
            include.append(copy.copy(k))
    return include

#build your parameter space given a set of p_types
def number_distance(a, b):
    if (a <= 0) and (b <= 0) or (a >= 0) and (b >= 0):
        sign = abs( abs(a) - abs(b) )
        return sign
    if (a <= 0) and (b >= 0) or (a >= 0) and (b <= 0):
        sign = abs( abs(a) + abs(b) )
        return sign
    
def check_distribution(boundaries):
    """lower and upper boundaries"""
    lower,upper = boundaries   
    try:
        lo = int(math.log10(lower))
    except ValueError:
        "There is a domain error in the lower bound! 0 does not have an order, change the boundaries of {}"
    try:
        lu = int(math.log10(upper))
    except ValueError:
        "There is a domain error in the upper bound! 0 does not have an order, change the boundaries of {}"    
        
    order = number_distance(lo,lu)
    if order > 2:
        return "log"
    else:
        return "uniform"
    
def define_parameter_boundaries(boundaries,spacing = 1000, distribution = "log_uniform_dist"):
    """"parameters that are supposed to be zero"""
    exclude = []
    """build the space"""
    space = {}
    for p,boundary in boundaries.items():
        if type(boundary) == str:
            boundary = eval(boundary)
        lower,upper = boundary
        
        """filter out the zero values set by accident"""
        if not lower and upper:
            lower = upper/100.
        if not upper and lower:
            upper = lower*100
        if not upper and not lower:
            lower,upper = 1,100
            exclude.append(p)
            
        try:
            lo = int(math.log10(lower))
        except ValueError:
            "There is a domain error in the lower bound! 0 does not have an order, change the boundaries of {}".format(p)
        try:
            lu = int(math.log10(upper))
        except ValueError:
            "There is a domain error in the upper bound! 0 does not have an order, change the boundaries of {}".format(p)    
            
        order = number_distance(lo,lu)
        if order == None:
            order = 1
        if upper == lower:
            space[p] = upper
        if order > 1.5:
            space[p] = d_func[distribution](lower,upper,spacing)
        else:
            space[p] = d_func["uniform_dist"](lower,upper,spacing)
    
    return space,exclude

#count the number of files in a directory
def countfiles(path):
    return len([name for name in os.listdir('.') if os.path.isfile(name)])
def filecount(folder):
    total = 0
    for root, dirs, files in os.walk(folder):
        total += len(files)
    return total

# build the equations, i.e. replace the states with indexed vectors e.g. y[0], y[1] etc
def build_equations(sm,states):
    for j in range(len(states)):
        states[j] = states[j].strip()
        sm = sm.replace(" "+states[j]+" ",' y[{}] '.format(j))
    if ',' in sm:
        parsed = sm.split(',')    
    elif '\n' in sm:
        parsed = sm.split('\n')
    parsed = [i for i in parsed if i]
    return parsed

#order a modelset by size
def OrderbySize(models):  
    modelset = {}
    number_of_states = [len(i.states) for i in models]
    """"the model and the length of the states"""
    for i in range(len(model)):
        modelset[i] = models[numpy.argsort(number_of_states).index(i)]
    """"the model and the length of the states"""
    return modelset

#build the complete model by updating the parameters using string replacement not p[1],p[2] etc.
def build_models(equations, p_set):
    model = []
    for i in equations:
        for k,v in p_set.items():
            i = i.replace(k,str(v))
        model.append(i)
    return model

#list the ranks of a list
def ranklist(array):
    """temp of the overall array"""
    temp = array.argsort()
    """"get the ranks of the array"""
    ranks = numpy.empty_like(temp)
    """fill the array"""
    ranks[temp] = numpy.arange(len(array))
    return ranks

#build a model with the parameters sets defined as indexed items in the network
#make sure to compile the model to prevent massive slowdowns (i.e. eval would have to eval every integration step instead of only once)
def build_base_model(equations,p_set):
    model = []
    for i in equations:
        for k,v in p_set.items():
            k = k.strip()
            i = i.replace(" "+k+" ",v)
        model.append(i.strip())
    stringmodel = copy.copy(model)
    for i in range(len(model)):
        model[i] = compile(model[i],'state_equations','eval')  
    return model,stringmodel

def plot_datavector(data):
    figure = plt.figure(figsize=(10,10))
    for ID,vector in data.items():
        plt.plot(vector,label = ID)
        plt.ylabel("[C]")
        plt.xlabel("Time")
    plt.legend(fancybox = True)
    plt.show()
    

def start_integration(func,p,ic,t0,t,integrator,atol,rtol,nsteps):
    r = scipy.integrate.ode(func).set_integrator(integrator,atol=atol,rtol=rtol,nsteps=nsteps,with_jacobian = False)
    r.set_f_params(p)
    points     = []
    timepoints = []
    r.set_initial_value(ic, t0)
    while r.successful() and r.t < t:
        values = r.integrate(r.t + 1)
        timepoints.append(r.t)
        points.append(values)   
    points = numpy.vstack(points) 
    t = numpy.array(timepoints) 
    return points,t

def update_parameters(model,parameters,event):
    pID = model.p_space.keys()
    for k,v in event.items():
        if parameters in pID:
            index = model.p_map[k]
            parameters[index] = v 
    return parameters
            
def update_ic(model,ic,event):
    for k,v in event.items():
        if k in model.states:
            ic[model.ic_map[k]] = v
    return ic
            
def filecount(path):
    return str(len(os.listdir(path)))

def strip_string(mystring):
    return int(re.search(r'\d+', mystring).group()) 

def all_subsets(ss):  
    from itertools import chain, combinations
    return chain(*map(lambda x: combinations(ss, x), range(0, len(ss)+1)))

def build_local_search(space,fixed,include = [],samples = 25):
    samplerange = [i*int(len(list(space.values())[-1])/samples) for i in range(samples)] 
    if not include:
        include = space.keys()
    local_space = {k:list(sorted(list(v.values()))) for k,v in space.items() if k in include}
    for k,v in local_space.items():
        local_space[k] = [v[i] for i in samplerange]
    return local_space

#set up a latin hypercube sample
def build_global_search(space,include = [],samples = 10000,order = 2,sobol = False): 
    if not include:
        include = space.keys()
    print(include)
    parameters = []     
    y = numpy.array(lhs(len(include),samples)*(10**order),dtype = int) # list van lists
    if sobol:
        y = numpy.array(sobol_seq.i4_sobol_generate(len(include), samples)*(10**order),dtype = int)
    for i in y:
        pset = {k:None for k in space.keys() if k in include}
        maps = list(pset.keys())
        for j in range(len(i)):
            pset[maps[j]] = space[maps[j]][i[j]]
        parameters.append(pset)
    return parameters

def pareto_frontier(Xs, Ys, maxX = True, maxY = True):
# Sort the list in either ascending or descending order of X
    myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
# Start the Pareto frontier with the first value in the sorted list
    p_front = [myList[0]]    
# Loop through the sorted list
    for pair in myList[1:]:
        if maxY: 
            if pair[1] >= p_front[-1][1]: 
                p_front.append(pair) 
        else:
            if pair[1] <= p_front[-1][1]: 
                p_front.append(pair) 
# Turn resulting pairs back into a list of Xs and Ys
    p_frontX = [pair[0] for pair in p_front]
    p_frontY = [pair[1] for pair in p_front]
    return p_frontX, p_frontY

def build_single_phase(space,fixed,include = [],samples = 4):
    if not include:
        include = space.keys()
    samplerange = [i*(len(space.values()[-1])/samples) for i in range(samples)] 
    phase = [i for i in space.keys() if i in include]
    combinations = list(itertools.product(samplerange,repeat = len(phase)))
    parameters = []
    for i in range(len(combinations)):
        template = copy.deepcopy(fixed)
        for j in range(len(combinations[i])):
            template[phase[j]] = space[phase[j]][combinations[i][j]]
        parameters.append(template)
    return parameters

def build_complete_phase(boundaries,fixed, include = [],samples = 4):
    parameters = [i for i in boundaries.keys() if i in include]
    phaseparameters = itertools.combinations(parameters,3)
    phase = {}    
    for i in phaseparameters:
        localbound = {j:boundaries[j] for j in i}
        space,exluded = define_parameter_boundaries(localbound)
        phase[i] = build_single_phase(space,fixed,include = [],samples = samples)
    return phase
        
def translate_periodicities(periodicities):
        return {i: 2*(1./i) for i in periodicities}

def sobol_sequence_sampling(samples = 1000,dimensions = 2):
    data = sobol_seq.i4_sobol_generate(dimensions, samples)
    random.shuffle(data)
    return data

def latin_hypercube_sampling(samples = 1000,dimensions = 2):
    data = lhs(dimensions,samples,criterion = 'corr')
    return data

def plot_input_database(data):
    for i in data.values():
        plt.plot(numpy.arange(0,len(i)),i)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
    plt.show()
    return

def sort_axes(x,y):
    """"get ranks"""
    ranks = numpy.argsort(x)
    """ sort x values"""
    new_x = list(sorted(x))
    """ sort y values """
    new_y = [y[i] for i in ranks]
    return new_x,new_y


def iproportional(periodicities = (50,100),samples = 10, boundaries = (0,1),time = 400,pattern = 'binary'):
    pi = {}
    minimum,maximum = periodicities
    periods = translate_periodicities(numpy.linspace(minimum,maximum,samples))
    t = numpy.arange(0,time,1)
    for period in numpy.linspace(minimum,maximum,samples):
        if pattern == 'binary':
             y = 0.5*signal.square(periods[period]*numpy.pi*t)+0.5
             pi[period] = y
        elif pattern == 'oscillating':    
             y = 0.5*(numpy.sin(periods[period]*numpy.pi*t))+0.5
             pi[period] = y            
        elif pattern == 'sawtooth':
             y = 0.5*signal.sawtooth(periods[period]*numpy.pi*t)+0.5
             pi[period] = y             
    return pi

def inputgeneration(ispace,time,fixed = [],period = (10,10),bounds = (0,1),pattern = 'binary', imin= 1,imax = 1, samples = 1, order = 100,random = False):
    '''inputs are always time dependent i..e an operational time where we modify the conditions of reactor or cell
    as it displays its dynamic behaviour'''
    ip = {i:[] for i in range(time)}
    if fixed:
        for item,ivec in ip.items():
            for mod in ivec:
                window,value = mod
                start,end = window
                for i in range(start,end,1):
                    ip[i].append((item,value))

    '''there are three different input conditions which can be modified to achieve a certain output'''
    if period:
        profile = iproportional(periodicity = period, samples = samples, boundaries = bounds,time = time,pattern = pattern)
        for k,v in ispace.item():
            order = int(math.floor(math.log10(len(v))))
            profile *= 10**order
            for i in range(len(profile)):
                ip[i].append((k,v[int(profile[i])]))
    return ip

def pickle_store(path,file_name,data):
    with open(path+file_name, 'wb') as handle:
        pickle.dump(data, handle,pickle.HIGHEST_PROTOCOL)
    return None

def store_class(data_obj,path = "~\\OneDrive\\Desktop\\",filename = "analysis"):
    path = os.path.expanduser(path)
    pickle_store(path,"result_{}",data_obj)
    return


def invert_dict_nonunique(d):
    newdict = {}
    for k, v in d.iteritems():
        newdict.setdefault(v, []).append(k)
    return newdict


def build_constructs(promoters,ORF,tags = []):
    constructs = []
    for i in promoters:
        for ii in ORF:
            if tags:
                for iii in tags:
                    constructs.append(i+"."+ii+":"+iii)
            constructs.append(i +'.'+ ii)
    return constructs

def invert_dict(item):
    new = defaultdict(list)
    for k,value in index.items():
        for v in value:
            new[v].append(k)  
    return new             

def flatten(l):
      return [item for sublist in l for item in sublist]     

def extractstring(l):
    l = list(set(l))
    stm = ''
    for i in l:
        stm += i
    return stm   

def strspace(s):
    return " " + s + " "   

def storedata(directory = '', name = '', data = ''):
    import pickle
    with open(directory + name, 'w') as handle:
        pickle.dump(data, handle,pickle.HIGHEST_PROTOCOL)
    return None

def create_commandlist(measurement):
    pulsescheme = list(measurement.values())[-1].time_dependent_parameters
    commandlist = {}
    for times,commands in pulsescheme.items():
        for control,value in commands.items():
            try:
                commandlist[control].update({times:value})
            except:
                commandlist[control] = {}
                commandlist[control].update({times:value}) 
    return commandlist 

def single_value_decomposition(sensitivity_vector,observables,pair = False):
    """store the variables"""
    store,cumsense = [],[]
    """combinations of all parameter classes"""
    if not pair:
        combinations = all_combinations(list(sensitivity_vector.keys()))
    """the combinations in the parameter set"""
    if pair:
        combinations = itertools.combinations(list(sensitivity_vector.keys()),2)
    """calculate the sensititivity vector"""
    for combination in combinations:
        S_k = numpy.transpose(tuple([sensitivity_vector[i] for i in combination]))
        
        arr,sv,uni = numpy.linalg.svd(S_k)
        try:
             collinearity = 1./float(abs(min(sv)))
        except ZeroDivisionError:
             collinearity = False
        store.append(collinearity_object(combination,collinearity,observables))
    for parameter,data in sensitivity_vector.items():
        cumsense.append(math.sqrt(numpy.sum(data*data)))
    return store,cumsense

def Fisher_D(sensitivity_vector,observables,pair = True,):
    pass
    
def Fisher_E(sensitivity_vector,observables,pair = True,):
    pass
    
def call_python_version(Version, Module, Function, ArgumentList):
    """call  python 2 function"""
    import execnet
    gw = execnet.makegateway("popen//python=python%s" % Version)
    channel = gw.remote_exec("""
        from %s import %s as the_function
        channel.send(the_function(*channel.receive()))
    """ % (Module, Function))
    channel.send(ArgumentList)
    return channel.receive()

def find_elements(s, ch):
    return [(i,ch) for i, ltr in enumerate(s) if ltr == ch and s[i+1] == ' ']

def maximum_likelihood_simulation(model,condition,observables = [],dt = 1,time = (0,50),show = False):
        from SolveSystem import ModelSolver
        from Model import ModelVariables
      
        """create the parameter vectors to solve the system"""
        variables = ModelVariables(model,modification = model.fixed,conditions = condition)
        """simulate the network"""

        solution = ModelSolver(model,variables = variables,simtime = time,dt = dt)
        data = solution.__getData__()

        if show:
            for i in range(len(model.states)):
                plt.plot(data.time,data.rdata[:,i],label = model.states[i])
            plt.legend(fancybox = True)
            plt.xlabel("Time")
            plt.ylabel("Concentration")
            plt.title("Simulation of a single user defined parameter set")
            plt.show()
        return data
    
def all_combinations(any_list):
    return itertools.chain.from_iterable(itertools.combinations(any_list, i + 1)for i in range(len(any_list)))
    
def window(iterable, size):
    iters = tee(iterable, size)
    for i in range(1, size):
        for each in iters[i:]:
            next(each, None)
    return zip(*iters)



