# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 16:05:07 2023

@author: huckg
"""

#build the reaction for the network
def categorize_reactions(network):
    import copy
    
    """template reaction reference"""
    template = {'R'  :[],
                'P'  :[],
                'E'  :False,
                'I'  :[],
                'A'  :[],
                'K'  :False}
    
    """reaction network"""
    reaction_network = []
    
    for i in network:
        """RPEIA, the individual components of the dataset"""
        E,rev,K,R,P,A,I = i
        
        print("A",A,'I',I)
        """define the reaction"""
        r = copy.deepcopy(template)
        
        """the overall reaction RPEAI"""
        for n in R:
            r['R'].append(n)
        for n in P:
            r['P'].append(n)
        r['E'] = E
        if K:
            r['K'] = K
        for n in I:
            r['I'].append(n)
        for n in A:
            r['A'].append(n)
        
        """appending the reaction to the reaction network"""
        reaction_network.append(copy.deepcopy(r))
        
        if rev == True:
            """define the reaction"""
            r = copy.deepcopy(template)
            """the overall reaction RPEAI"""
            for n in R:
                r['P'].append(n)
            for n in P:
                r['R'].append(n)
            r['E'] = E
            if K:
                if len(P) > 1 and len(R) > 1:
                    r['K'] = K
                elif len(P) == 1:
                    r['K'] = 'MM'
                elif len(P) > 1 and len(R) == 1:
                    r['K'] = 'GH'
                else:
                    r['K'] = K              
            for n in I:
                pass
                # r['I'].append(n)
            for n in A:
                pass
                # r['A'].append(n)
            
            """appending the reaction to the reaction network"""
            print(r)
            reaction_network.append(copy.deepcopy(r))
    return reaction_network

def construct_model_fluxterms(network,reactormix = {}):
    """import the modules"""
    from BuildModelModule import ReactionKinetics
    from BuildModelModule import ReactorKinetics
    
    """Get all the kinetics that could apply"""
    members = [func for func in dir(ReactionKinetics) if callable(getattr(ReactionKinetics, func)) if 'kinetic' in func]

    """This funciton takes in a list of reaction equations to build model"""
    reaction_kinetics = []

    """Build the reaction network"""
    for r in network:
        flux = ReactionKinetics(r['E'],r['R'],r['P'],r['I'],r['A'])
        if r["K"] == False:
            for kinetic in members:
               exec('flux.{}()'.format(kinetic))
        else:
            exec('flux.kinetic_{}()'.format(r["K"]))
        reaction_kinetics.append(flux)
        
    """Extract the states from the reaction"""
    states = []
    for i in reaction_kinetics:
        states += i.states
    states = list(set(sorted([i for i in states if i])))
    
    """flux equations for controlled species"""
    reactor_kinetics = ReactorKinetics(states,reactormix=reactormix)

    """Get all the kinetics that could apply"""
    members = [func for func in dir(ReactorKinetics) if callable(getattr(ReactorKinetics, func)) if 'reactor' in func]
    for reactor in members:
        exec('reactor_kinetics.{}()'.format(reactor))
    return (reactor_kinetics,reaction_kinetics)

def construct_equation(reactants,products,fluxterm):
    """this function takes the products and reactants
    it subsequently ties them to the flux to build an equation
    the reactants and products are defined a [S1,S2],[P1,P2],
    the category S or P defines whether its a negative or postive flux
    the flux term is a string (a + b*c / b + etc etc.)""" 
    equations = {i:[] for i in reactants+products}
    """define the terms"""
    for i in reactants:
        equations[i].append( ' - {0} '.format(fluxterm))
    for i in products:
        equations[i].append( ' + {0} '.format(fluxterm))        
    return equations

def construct_kinetic_combinations(reactor_kinetics,reaction_kinetics,name = 'Test_1', reactor_type = 'reactor_predefined_control_single_reactor'):
    """import itertools to make model combinations"""
    import itertools as it
        
    """Get the combination of reaction kinetics"""
    combinations = {}
    for i in range(len(reaction_kinetics)):
        combinations[i] = []
        for n in reaction_kinetics[i].equations.keys():
            combinations[i].append(n)
            
 
    """The metamodel set gives the list of the dataset"""
    metamodel_set = []
    
    """Posssible kinetic combinations for the reactor"""
    c = {}
    for i in range(len(combinations)):
        c[i] = []
        for n in combinations[i]:
            c[i].append((i,n))

    """Get the combinations of the IT.product"""
    chosen_models = list(it.product(*list(c.values()))) 
    print('There are {} possible models: please proceed with caution, is this correct if so press 1 \n'.format(str(len(chosen_models))))

    """The manually defined input"""
    if len(chosen_models) == 1:
        """The reactor and kinetics"""
        fluxterms = len(reaction_kinetics)
        for i in range(len(chosen_models)):
            """This is the model,The control parameters, Kinetic parameters"""
            model_placeholder   = []
            """3 classes of parameters, control, known parameters and kinetic parameters, to be sorteed by parameter class family in modelbuild module"""
            kinetic_parameters  = []
    
            """Append the individual fluxes within a single dict"""
            for r,rtype in chosen_models[i]:
                model_placeholder.append(reaction_kinetics[r].equations[rtype])
                kinetic_parameters.extend(reaction_kinetics[r].kinetic_parameters[rtype])
            """the model in the dataset"""
            model_placeholder.append(reactor_kinetics.equations[reactor_type])
            """Add the control parameters form the reactors"""
            kinetic_parameters.extend(reactor_kinetics.kinetic_parameters[reactor_type])    

            """Initialize an empty result dictionary"""
            model_dict = {}
            
            """Loop through the dictionaries and append lists to the result_dict"""
            for d in model_placeholder:
                for key, value in d.items():
                    if key in model_dict:
                        model_dict[key].extend(value)
                    else:
                        model_dict[key] = value
                       
            """Define the base model in the class"""
            from BuildModelModule import base_Model 
            base = base_Model(chosen_models[i],
                              model_dict,
                              kinetic_parameters,
                              name     = name)
            
            """Modelset that contains all the models that have been generated"""
            return base  
    else:
        print('please check if you generate singular model or all possible variants, if you want to build all variants press 1')
        build_all_variants = input()
        
        if build_all_variants:
            """The reactor and kinetics"""
            fluxterms = len(reaction_kinetics)
    
        """The reactor and kinetics"""
        fluxterms = len(reaction_kinetics)

        modelset = [] 
        for i in range(len(chosen_models)):
            """This is the model,The control parameters, Kinetic parameters"""
            model_placeholder   = []
            """3 classes of parameters, control, known parameters and kinetic parameters, to be sorteed by parameter class family in modelbuild module"""
            kinetic_parameters  = []
    
            """Append the individual fluxes within a single dict"""
            for r,rtype in chosen_models[i]:
                model_placeholder.append(reaction_kinetics[r].equations[rtype])
                kinetic_parameters.extend(reaction_kinetics[r].kinetic_parameters[rtype])
                
            """the model in the dataset"""
            model_placeholder.append(reactor_kinetics.equations[reactor_type])
            """Add the control parameters form the reactors"""
            kinetic_parameters.extend(reactor_kinetics.kinetic_parameters[reactor_type])    

            """Initialize an empty result dictionary"""
            model_dict = {}
            
            """Loop through the dictionaries and append lists to the result_dict"""
            for d in model_placeholder:
                for key, value in d.items():
                    if key in model_dict:
                        model_dict[key].extend(value)
                    else:
                        model_dict[key] = value
                       
            """Define the base model in the class"""
            from BuildModelModule import base_Model 
            base = base_Model(chosen_models[i],
                              model_dict,
                              kinetic_parameters,
                              name     = name + str(i))
            modelset.append(base)
        return modelset
