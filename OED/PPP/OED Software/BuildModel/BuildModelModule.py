# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 15:48:14 2023

@author: huckg
"""

"""Equation builder function"""
from ModelBuilderOperations import construct_equation

"""Import the python modules"""
import math

class ReactorKinetics:    
    def __init__(self,states,reactormix = {}):
        """every state gets an inflow and outflow term,
        the input file takes these inflow terms and matches
        them to the inputs that are present!"""
        self.controlled_kinetics = {}
        
        """get the states"""
        self.states =  states
        
        """store the parameters that are generated:
            we make a distinction between the parameters and the control parameters
            in this instance the parameters are actually known prior to the start of the
            experiment like volume and stock concentration"""
        self.kinetic_parameters =  {}

        """plug the equations into the reactor"""
        self.equations = {}
        
        """if there is a reactormix please define it"""
        self.reactormix = reactormix
                    
    def reactor_uniform_control(self,):
        """add category for parameters"""
        self.kinetic_parameters['reactor_uniform_control'] = []
        
        """inflow and outflow"""
        inflow = {}
        outflow= {}
        
        """self.states in the data"""
        for n in self.states:
            """The n species with in the channel"""
            kf   = 'kflow'
            kin  = '{0}(in)'.format(n)
            term = ' ( {0} * {1} ) '.format(kf,kin)

            """The reaction inflow"""
            inflow[n]  = term
            
            """Premix outflow"""
            k_out      = ' ( kflow * {0} ) '.format(n)
            
            """the outflow of te reaction mixture"""
            outflow[n] = k_out
            
            """the channel plus or minus 1"""
            self.kinetic_parameters['reactor_uniform_control'].extend([kf,kin])
            
        """For this function we will write the equations already"""
        equations = {}
        
        """Inflow terms"""
        for state,flux in inflow.items():
            if state not in equations:
                equations[state] = []
            equations[state].append( ' + ' + flux)
                
        """Outflow terms"""
        for state,flux in outflow.items():
            if state not in equations:        
                equations[state] = []
            equations[state].append( ' - ' + flux) 
        
        """update the equations"""
        self.equations['reactor_uniform_control'] = equations
            
    def reactor_individual_control(self,):
        """add category for parameters"""
        self.kinetic_parameters['reactor_individual_control'] = []

        """Every reactor has two states pre-mixing and post_mixing:
        premixed states"""
        premixed = {i:'pre_{}'.format(i) for i in self.states}
        
        """These premixed states travel trough a channel and combine
        to react within the reactor: premixed terms going in, premix flowing out in thereactor and outflow"""
        premix_inflow  = {}
        premix_outflow = {}
        
        """inflow and outflow"""
        inflow  = {}
        outflow = {}
        
        
        total_flow = '( '
        
        """start of the flowterm"""
        start = True
        
        """channel that flows into the reactor"""
        channel = 1
        for n in self.states:
            c = 'channel_{}'.format(channel)
            """Map the flow into the reactor"""
            if not start:
                total_flow += '|+| kflow({0}) '.format(c)
            else:
                total_flow += 'kflow({0}) '.format(c)
            
            """update the channel"""
            channel += 1
            
            """one iteration has been had, it is no longer the start"""
            start = False
            
        """the flow terms for the individual species"""
        channel = 1
        for n in self.states:
            """the channel of the reaction i.e. syrynge"""
            c = 'channel_{}'.format(channel)
            
            """parameters of the controlled parameters"""
            flow   = 'kflow({0})'.format(c)
            stock  = '{0}(in)'.format(n)
            volume = 'volume_{}'.format(c)
            
            """The n species with in the channel"""
            term   = '( ( {0} * {1} ) / {2} )'.format(flow,stock,volume)
            
            """The reaction inflow"""
            premix_inflow[premixed[n]]  = term
            
            """Premix outflow"""
            term = '( ( {0} * {1} ) / {2} )'.format(flow,premixed[n],volume)
            premix_outflow[premixed[n]] = term
            
            """the channel plus or minus 1"""
            self.kinetic_parameters['reactor_individual_control'].extend([flow,volume,stock])

            """the channel plus or minus 1"""
            channel += 1
            
        """Close the total flow"""
        total_flow += ' )'
        
        """the secondary reactor also has a reactor volume"""
        reactor_volume = 'reactor_volume'
        """Map the flow out of the reactor, its the sum of the individual flow rates"""
        for state in self.states:
            outflow[state] = ' ( ( {0} * {1} ) / {2} ) '.format(total_flow ,state, reactor_volume)
            
        """add this parameter"""
        self.kinetic_parameters['reactor_individual_control'].append(reactor_volume)         
        
        """Get the inflow of the mixed species"""
        for state in self.states:
            inflow[state] = premix_outflow[premixed[state]]
            
        """for this function we will write the equations already"""
        equations = {}
        
        """inflow terms"""
        for i in [inflow,premix_inflow]:
            for state,flux in i.items():
                if state not in equations:
                    equations[state] = []
                equations[state].append( ' + ' + flux)
                
        """outflow terms"""
        for i in [outflow,premix_outflow]:
            for state,flux in i.items():
                if state not in equations:        
                    equations[state] = []
                equations[state].append( ' - ' + flux) 
              
        """update the equations"""
        self.equations['reactor_individual_control'] = equations
       
    def reactor_predefined_control_single_reactor(self,):
        """if the reactor mix is not given it defaults to a 
        system where each species has its own syringe"""
        if not self.reactormix:
            self.reactormix = {i:[self.states[i]] for i in range(len(self.states))}

        """add category for parameters"""
        self.kinetic_parameters['reactor_predefined_control_single_reactor'] = []

        """These premixed states travel trough a channel and combine
        to react within the reactor, premixed terms going in, premix flowing out in thereactor and outflow"""
        inflow  = {}
        outflow = {}
        
        """the flowterm for the total flow within the reactor"""
        total_flow = '( '
        
        """start of the flow reaction"""
        start = True
        for channel,species in self.reactormix.items():
            """the channel of the reaction i.e. syrynge"""
            c = 'channel_{}'.format(channel)

            """two parametes the volume in channel and flow"""
            flow   = 'kflow({0})'.format(c)
            
            """the channel plus or minus 1"""
            self.kinetic_parameters['reactor_predefined_control_single_reactor'].extend([flow])
            
            """Map the flow into the reactor"""
            if not start:
                total_flow += '|+| {} '.format(flow)
            else:
                total_flow += '{} '.format(flow)
            start = False
            
            """the flow terms for the individual species"""
            for n in species:
                """parameters of the controlled parameters"""
                stock  = '{0}(in)'.format(n)
                volume = 'reactor_volume'
                
                """Premixed inflow"""
                term = '( ( {0} * {1} ) / {2} )'.format(flow,stock,volume)
                inflow[n] = term
                
                """update the stock concentration parameters"""
                self.kinetic_parameters['reactor_predefined_control_single_reactor'].append(stock)
                               
        """Close the total flow"""
        total_flow += ' )'

        """Map the flow out of the reactor, its the sum of the individual flow rates"""
        reactor_volume = 'reactor_volume'
        for state in self.states:
            outflow[state] = ' ( ( {0} * {1} ) / {2} ) '.format(total_flow ,state,reactor_volume)
          
        "update the volume of the reactor (parameter)"
        self.kinetic_parameters['reactor_predefined_control_single_reactor'].append(reactor_volume)
        
            
        """for this function we will write the equations already"""
        equations = {}
        
        """inflow terms"""
        for state,flux in inflow.items():
            if state not in equations:
                equations[state] = []
            equations[state].append( ' + ' + flux)
                
        """outflow terms"""
        for state,flux in outflow.items():
            if state not in equations:        
                equations[state] = []
            equations[state].append( ' - ' + flux) 
                
        """update the equations"""
        self.equations['reactor_predefined_control_single_reactor'] = equations        
       
        
    def reactor_predefined_control_staged_reactor(self,):
        """if the reactor mix is not given it defaults to a 
        system where each species has its own syringe"""
        if not self.reactormix:
            self.reactormix = {i:[self.states[i]] for i in range(len(self.states))}
            """Every reactor has two states pre-mixing and post_mixing:
            premixed states:"""
            premixed = {i:'pre_{}'.format(i) for i in self.states}
        else:
            """premixed labels"""
            premixed = {}
            for channel,species in self.reactormix.items():
                for state in species:
                    premixed[state] = 'pre_{}'.format(state) 

        """add category for parameters"""
        self.kinetic_parameters['reactor_predefined_control_staged_reactor'] = []

        """These premixed states travel trough a channel and combine
        to react within the reactor, premixed terms going in, premix flowing out in thereactor and outflow"""
        premix_inflow  = {}
        premix_outflow = {}
        inflow  = {}
        outflow = {}
        
        """the flowterm for the total flow within the reactor"""
        total_flow = '( '
        
        """start of the flow reaction"""
        start = True
        for channel,species in self.reactormix.items():
            """the channel of the reaction i.e. syrynge"""
            c = 'channel_{}'.format(channel)

            """two parametes the volume in channel and flow"""
            flow   = 'kflow({0})'.format(c)
            volume = 'volume_{}'.format(c)
            
            """the channel plus or minus 1"""
            self.kinetic_parameters['reactor_predefined_control_staged_reactor'].extend([flow,volume])
            
            """Map the flow into the reactor"""
            if not start:
                total_flow += '|+| {} '.format(flow)
            else:
                total_flow += '{} '.format(flow)
            start = False
            
            """the flow terms for the individual species"""
            for n in species:
                """parameters of the controlled parameters"""
                stock  = '{0}(in)'.format(n)

                """The n species with in the channel"""
                term = '( ( {0} * {1} ) / {2} )'.format(flow,stock,volume)
                
                """The reaction inflow"""
                premix_inflow[premixed[n]] = term
                
                """Premix outflow"""
                term = '( ( {0} * {1} ) / {2} )'.format(flow,premixed[n],volume)
                premix_outflow[premixed[n]] = term
                
                """update the stock concentration parameters"""
                self.kinetic_parameters['reactor_predefined_control_staged_reactor'].append(stock)
                               
        """Close the total flow"""
        total_flow += ' )'

        """Map the flow out of the reactor, its the sum of the individual flow rates"""
        reactor_volume = 'reactor_volume'
        for state in self.states:
            outflow[state] = ' ( ( {0} * {1} ) / {2} ) '.format(total_flow ,state,reactor_volume)
          
        "update the volume of the reactor (parameter)"
        self.kinetic_parameters['reactor_predefined_control_staged_reactor'].append(reactor_volume)
        
        """Get the inflow of the mixed species"""
        for channel,species in self.reactormix.items():
            for state in species:
                inflow[state] = premix_outflow[premixed[state]]
            
        """for this function we will write the equations already"""
        equations = {}
        
        """inflow terms"""
        for i in [inflow,premix_inflow]:
            for state,flux in i.items():
                if state not in equations:
                    equations[state] = []
                equations[state].append( ' + ' + flux)
                
        """outflow terms"""
        for i in [outflow,premix_outflow]:
            for state,flux in i.items():
                if state not in equations:        
                    equations[state] = []
                equations[state].append( ' - ' + flux) 
                
        """update the equations"""
        self.equations['reactor_predefined_control'] = equations
        
class ReactionKinetics:
    def __init__(self,E,R,P,I,A):
        """this class simply stores the information from 
        the modelbuilder, since the same reaction can have different
        variants we include this model here"""
        self.reactants  = R
        self.products   = P
        self.inhibitors = I
        self.activators = A
        self.enzyme     = E
        
        """Get the kinetic parameters"""
        self.kinetic_parameters = {}
        
        """Get the kinetic flux and set the equation"""
        self.kinetic   = {}
        self.equations = {}
     
        """States in this equations"""
        self.states = [self.enzyme] + self.reactants + self.products + self.inhibitors + self.activators
        
        """Extract some general reaction rules from the vector"""
        enzymatic_reaction = True
        if self.enzyme == False:
            k = 'k_fwd_'
            for S in self.reactants:
                """if there is no enzymatic reactoin we will assume there is a 
                general mass action mechanism that is taking place to drive the reaction"""
                k += S + '_' 
            """delete last underscore"""
            k = k[:-1]
            
            term = k + ' * '
            for S in self.reactants:
                """the flux term equation for mass action"""
                term += S + ' * '
            """delete last multiplication marker"""
            term = term[:-3]
            
            """update the kinetic parameters"""
            if 'MA' not in self.kinetic_parameters:
                self.kinetic_parameters['MA'] = []
            self.kinetic_parameters['MA'].append(k)
            
            """update the kinetic"""
            self.kinetic['MA'] = term
    
            """update the equations"""
            self.equations['MA']  = construct_equation(self.reactants,self.products,term)

        self.allosteric_terms = []
        """the allosteric interaction
           inhibition term:"""
        self.a_inh = '' 
        if self.inhibitors:
            """build allostery equation"""
            for i in self.inhibitors:
                k_inh   = 'k_inh_({0}:{1})'.format(self.enzyme,i)
                k_N_inh = 'k_N_inh({0}:{1})'.format(self.enzyme,i)
                term    = " ( ( 1 |+| ( {0} / {1} ) ) ** {2} ) ".format(i,k_inh,k_N_inh)
                """Track allosteric parameters"""
                self.allosteric_terms.append(k_inh)
                self.allosteric_terms.append(k_N_inh)
                
                """update the term"""
                self.a_inh += ' * ' + term  
         
        """activation term"""
        self.a_act = '' 
        if self.activators:
            for i in self.activators:
                act     = 'k_act_({0}:{1})'.format(self.enzyme,i)
                k_N_act = 'k_N_act({0}:{1})'.format(self.enzyme,i)
                term    = " ( 1 |+| ( {0} / {1} ) ) ".format(i,act,k_N_act)
                
                """Track allosteric parameters"""
                self.allosteric_terms.append(act)
                self.allosteric_terms.append(k_N_act)
                
                """update the term"""
                self.a_act += ' * ' + term  
                
    def kinetic_MM(self,):
        if len(self.reactants) == 1 and self.enzyme != False:  
            """unpack the substrate"""
            S1 = self.reactants[0]
            
            """define two parameters in MM equation"""
            km1  = 'k_m_({0}:{1})'.format(self.enzyme,S1)
            kcat = 'k_cat_({0}:{1})'.format(self.enzyme,S1)
            
            """update parameter set"""
            if 'MM' not in self.kinetic_parameters:
                self.kinetic_parameters['MM'] = []
            self.kinetic_parameters['MM'].extend([km1,kcat])
            if self.allosteric_terms:
                self.kinetic_parameters['MM'].extend(self.allosteric_terms)
            
            """define the flux term"""
            MM_numerator   = ' ( ( {1} * {2} ) * {3} ) {0}  '.format(self.a_act,self.enzyme,km1,S1)
            MM_denominator = ' ( {0} {1} |+| {2} )  '.format(km1,self.a_inh,S1)         
            
            """flux term"""
            self.MM_flux = '( (' + MM_numerator + ')' + ' / ' + '(' + MM_denominator + ') )'
            
            """update kinetoic"""
            self.kinetic['MM'] = self.MM_flux
            
            """update the equations"""
            self.equations['MM']  = construct_equation(self.reactants,self.products,self.MM_flux)
  
    def kinetic_GH(self,):
        """bi-reactant flux equations generalized hill equations"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
        
            """the km1 and km2 and catalysis rate"""
            km1  = 'k_m_({0}:{1})'.format(self.enzyme,S1) 
            km2  = 'k_m_({0}:{1})'.format(self.enzyme,S2)  
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2) 
            
            """update parameter set"""
            if 'GH' not in self.kinetic_parameters:
                self.kinetic_parameters['GH'] = []
            self.kinetic_parameters['GH'].extend([km1,km2,kcat])
            if self.allosteric_terms:
                self.kinetic_parameters['GH'].extend(self.allosteric_terms)
            
            """the nominator and denominator terms for the flux"""
            GH_numerator   = '( ( {1} * {2} ) * {3} * {4} ) {0} '.format(self.a_act,self.enzyme,kcat,S1,S2)
            GH_denominator = ' ( ( {0} * {1} ) |+| ( {3} * {4} ) |+| ( {5} * {6} ) |+| ( {7} * {8} ) ){2} '.format(km1,km2,self.a_inh,km2,S1,km1,S2,S1,S2)
              
            """flux term"""
            self.GH_flux = ' (' + GH_numerator + ')' + ' / ' + '(' + GH_denominator + ') '
            
            """update kinetoic"""
            self.kinetic['GH'] = self.GH_flux
            
            """update the equations"""
            self.equations['GH'] = construct_equation(self.reactants,self.products,self.GH_flux)
            
    def kinetic_GHB(self,):
        """bi-reactant flux equations generalized hill equations"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
        
            """the km1 and km2 and catalysis rate"""
            km1  = 'k_m_({0}:{1})'.format(self.enzyme,S1) 
            km2  = 'k_m_({0}:{1})'.format(self.enzyme,S2)  
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2) 
            kia  = 'k_inh_({0}:{1})'.format(self.enzyme,S1)
            kib  = 'k_inh_({0}:{1})'.format(self.enzyme,S2)
            
            """update parameter set"""
            if 'GHB' not in self.kinetic_parameters:
                self.kinetic_parameters['GHB'] = []
            self.kinetic_parameters['GHB'].extend([km1,km2,kcat,kib,kia])
            if self.allosteric_terms:
                self.kinetic_parameters['GHB'].extend(self.allosteric_terms)
            
            """the nominator and denominator terms for the flux"""
            GHB_numerator   = ' ( ( ( {1} * {2} ) * {3} * {4} ) {0} ) '.format(self.a_act,self.enzyme,kcat,S1,S2)
            GHB_denominator = ' ( ( {0} * {1} ) * ( 1 |+| ( {2} / {0} ) ) * ( ( 1 |+| ( {3} / {1} ) ) ) * ( ( 1 |+| ( {2} / {4} ) ** 2) ) * ( ( 1 |+| ( {3} / {6} ) ** 2 ) ) {5} )  '.format(km1,km2,S1,S2,kia,self.a_inh,kib)
              
            """flux term"""
            self.GHB_flux = ' (' + GHB_numerator + ')' + ' / ' + '(' + GHB_denominator + ') '
            
            """update kinetic"""
            self.kinetic['GHB'] = self.GHB_flux
            
            """update the equations"""
            self.equations['GHB']  = construct_equation(self.reactants,self.products,self.GHB_flux)
            
    def kinetic_MWC(self,):
        """bi-reactant flux equations generalized hill equations"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
        
            """the km1 and km2 and catalysis rate"""
            km1  = 'k_m_({0}:{1})'.format(self.enzyme,S1) 
            km2  = 'k_m_({0}:{1})'.format(self.enzyme,S2)  
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2) 
            Na   = 'k_N_({0}:{1})'.format(self.enzyme,S1)
            Nb   = 'k_N_({0}:{1})'.format(self.enzyme,S2)
            
            """update parameter set"""
            if 'MWC' not in self.kinetic_parameters:
                self.kinetic_parameters['MWC'] = []
            self.kinetic_parameters['MWC'].extend([km1,km2,kcat,Na,Nb])
            if self.allosteric_terms:
                self.kinetic_parameters['MWC'].extend(self.allosteric_terms)
            
            """the nominator and denominator terms for the flux"""
            MWC_numerator   = ' ( ( ( {1} * {2} ) * {3} * {4} ) {0} ) '.format(self.a_act,self.enzyme,kcat,S1,S2)
            MWC_denominator = ' ( ( {0} * {1} ) *  ( 1 |+| ( {2} / {0} ) ** {4} ) * ( ( 1 |+| ( {3} / {1} ) ** {6} ) ) {5} )  '.format(km1,km2,S1,S2,Na,self.a_inh,Nb)
              
            """flux term"""
            self.MWC_flux = ' (' + MWC_numerator + ')' + ' / ' + '(' + MWC_denominator + ') '
            
            """update kinetoic"""
            self.kinetic['MWC'] = self.MWC_flux
            
            """update the equations"""
            self.equations['MWC']  = construct_equation(self.reactants,self.products,self.MWC_flux)
            

           

    def kinetic_EO(self,order = ()):
        """bi-reactant flux equations, equilibrium ordered"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
            if order != ():
                S1,S2 = order
            
            """the kai and kb and catalysis rate"""
            kai  = 'k_ai_({0}:{1}_{2})'.format(self.enzyme,S1,S2)   
            kb   = 'k_b_({0}:{1})'.format(self.enzyme,S1)   
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2)  
            
            """update parameter set"""
            if 'EO' not in self.kinetic_parameters:
                self.kinetic_parameters['EO'] = []
            self.kinetic_parameters['EO'].extend([kai,kb,kcat])
            if self.allosteric_terms:
                self.kinetic_parameters['EO'].extend(self.allosteric_terms)
            
            """the nominator and denominator terms for the flux"""
            EO_numerator   = ' ( ( {1} * {2} ) * {3} * {4} ) {0} '.format(self.a_act,self.enzyme,kcat,S1,S2)
            EO_denominator = ' ( ( {0} * {1} ) |+| ( {3} * {4} ) |+| ( {5} * {6} ) ){2} '.format(kai,kb,self.a_inh,kb,S1,S1,S2)
              
            """flux term"""
            self.EO_flux = ' (' + EO_numerator + ')' + ' / ' + '(' + EO_denominator + ') '
            
            """update kinetoic"""
            self.kinetic['EO'] = self.EO_flux
            
            """update the equations"""
            self.equations['EO']  = construct_equation(self.reactants,self.products,self.EO_flux)
            
    def kinetic_PP(self,):
        """bi-reactant flux equations"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
            
            """the ka and kb and catalysis rate"""
            ka   = 'k_a_({0}:{1})'.format(self.enzyme,S1)
            kb   = 'k_b_({0}:{1})'.format(self.enzyme,S2) 
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2) 
            
            """update parameter set"""
            if 'PP' not in self.kinetic_parameters:
                self.kinetic_parameters['PP'] = []
            self.kinetic_parameters['PP'].extend([ka,kb,kcat])
            if self.allosteric_terms:
                self.kinetic_parameters['PP'].extend(self.allosteric_terms)
            
            """the nominator and denominator terms for the flux"""
            PP_numerator   =  ' ( ( {1} * {2} ) * {3} * {4} ) {0} '.format(self.a_act,self.enzyme,kcat,S1,S2)
            PP_denominator =  ' ( ( {0} * {1} ) |+| ( {3} * {4} ) |+|  ( {5} * {6} ) ){2} '.format(kb,S1,self.a_inh,ka,S2,S1,S2)
        
            """flux term"""
            self.PP_flux = ' (' + PP_numerator + ')' + ' / ' + '(' + PP_denominator + ') '
            
            """update kinetoic"""
            self.kinetic['PP'] = self.PP_flux
            
            """update the equations"""
            self.equations['PP']  = construct_equation(self.reactants,self.products,self.PP_flux)

    def kinetic_RER(self,):
        """bi-reactant flux equations for rapid equilibrium random reactions"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
        
            """4 Reactions in the rapid random equilibrium reaction"""
            kai  = 'k_ai_({0}:{1}_{2})'.format(self.enzyme,S1,S2)
            ka   = 'k_a_({0}:{1})'.format(self.enzyme,S1)   
            kb   = 'k_b_({0}:{1})'.format(self.enzyme,S2)    
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2)  
            
            """Update parameter set"""
            if 'RER' not in self.kinetic_parameters:
                self.kinetic_parameters['RER'] = []
            self.kinetic_parameters['RER'].extend([kai,ka,kb,kcat])
            if self.allosteric_terms:
                self.kinetic_parameters['RER'].extend(self.allosteric_terms)
            
            """The nominator and denominator terms for the flux"""
            RER_numerator   = ' ( ( {1} * {2} ) * {3} * {4} ) {0} '.format(self.a_act,self.enzyme,kcat,S1,S2)
            RER_denominator =  ' ( ( {0} * {1} ) |+| ( {3} * {4} ) |+| ( {5} * {6} ) |+| ( {7} * {8} ) ){2} '.format(kai,kb,self.a_inh,kb,S1,ka,S2,S1,S2)
                  
            """Define flux term"""
            self.RER_flux = ' (' + RER_numerator + ')' + ' / ' + '(' + RER_denominator + ') '
            
            """Update kinetic"""
            self.kinetic['RER'] = self.RER_flux
            
            """update the equations"""
            self.equations['RER']  = construct_equation(self.reactants,self.products,self.RER_flux)
            
    def kinetic_BH(self,):
        """bi-reactant flux equations for rapid equilibrium random reactions"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
        
            """4 Reactions in the rapid random equilibrium reaction"""
            ka   = 'k_a_({0}:{1})'.format(self.enzyme,S1)   
            kb   = 'k_b_({0}:{1})'.format(self.enzyme,S2) 
            
            """The large K values"""
            k1f   = 'k_1f_({0}:{1})'.format(self.enzyme,S2)   
            k2f   = 'k_2f_({0}:{1})'.format(self.enzyme,S1)   
            k1r   = 'k_1r_({0}:{1})'.format(self.enzyme,S2)   
            k2r   = 'k_2r_({0}:{1})'.format(self.enzyme,S1)    

            
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2)  
            
            """Update parameter set"""
            if 'BH' not in self.kinetic_parameters:
                self.kinetic_parameters['BH'] = []
            self.kinetic_parameters['BH'].extend([ka,kb,kcat,k1f,k2f,k1r,k2r])
            if self.allosteric_terms:
                self.kinetic_parameters['BH'].extend(self.allosteric_terms)
            
            """The nominator and denominator terms for the flux"""
            BH_numerator   = ' ( ( {1} * {2} ) * {3} * {4} ) {0} '.format(self.a_act,self.enzyme,kcat,S1,S2)
            BH_denominator =  ' ( ( {0} * ( 1 |+| ( {1} / {2} ) ) ) |+| ( {3} * ( 1 |+| ( {4} / {5} ) ) ) |+| ( {6} * ( 1 |+| ( {1} / {2} ) ) ) |+| ( {7} * ( 1 |+| ( {4} / {5} ) ) ) ) {8} '.format(k1f,S2,kb,k2f,S1,ka,k1r,k2r,self.a_inh)
                  
            """Define flux term"""
            self.BH_flux = ' (' + BH_numerator + ')' + ' / ' + '(' + BH_denominator + ') '
            
            """Update kinetic"""
            self.kinetic['BH'] = self.BH_flux
            
            """update the equations"""
            self.equations['BH']  = construct_equation(self.reactants,self.products,self.BH_flux)
            
    def kinetic_Hill(self,):
        """bi-reactant flux equations for rapid equilibrium random reactions"""
        if len(self.reactants) > 1 and self.enzyme != False:  
            S1,S2 = tuple(sorted(self.reactants))
        
            """4 Reactions in the rapid random equilibrium reaction"""
            ka   = 'k_a_({0}:{1})'.format(self.enzyme,S1)   
            kb   = 'k_b_({0}:{1})'.format(self.enzyme,S2)    
            kcat = 'k_cat_({0}:{1}_{2})'.format(self.enzyme,S1,S2)  
            Na   = 'k_N_({0}:{1})'.format(self.enzyme,S1) 
            Nb    ='k_N_({0}:{1})'.format(self.enzyme,S2) 
            
            
            """Update parameter set"""
            if 'Hill' not in self.kinetic_parameters:
                self.kinetic_parameters['Hill'] = []
            self.kinetic_parameters['Hill'].extend([ka,kb,kcat,Na,Nb])
            if self.allosteric_terms:
                self.kinetic_parameters['Hill'].extend(self.allosteric_terms)
            
            """The nominator and denominator terms for the flux"""
            Hill_numerator   = ' ( ( {1} * {2} ) * {3} * {4} ) {0} '.format(self.a_act,self.enzyme,kcat,S1,S2)
            Hill_denominator =  ' ( ( {0} ** {1} * {2} ** {3} ) |+| ( {4} ** {1} * {5} ** {3} ) |+| ( {4} ** {1} * {2} ** {3} ) ) {6} '.format(ka,Na,S2,Nb,S1,kb,self.a_inh)
                  
            """Define flux term"""
            self.Hill_flux = ' (' + Hill_numerator + ')' + ' / ' + '(' + Hill_denominator + ') '
            
            """Update kinetic"""
            self.kinetic['Hill'] = self.Hill_flux
            
            """update the equations"""
            self.equations['Hill']  = construct_equation(self.reactants,self.products,self.Hill_flux)
            

class K_:
    def __init__(self,name,control = False,known = False):
        """name of the parameter and the reaction it occures in """
        self.name    = name
        
        """let the user know if it is a control value"""
        self.control = control
        
        """let the user know that this is a value that needs to be provided manually"""
        self.known   = known
    
class Km(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (0.1,20000)
        self.unit     = 'uM'
        
class Ka(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (0.1,20000)
        self.unit     = 'uM'
        
class Kb(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (0.1,20000)
        self.unit     = 'uM'
        
class Kai(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (0.1,20000)
        self.unit     = 'uM'
        
class Kfwd(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (0.1,20000)
        self.unit     = 'min'
        
class Kcat(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (1,45000)
        self.unit     = 'min'
        
class Kinh(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (150,40000)
        self.unit     = 'uM'
        
class Kact(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (0.1,20000)        
        self.unit     = 'uM'
        
class Kn(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (1.0,1.5)        
        self.unit     = 'Dimless'
   
class Kstock(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = True)
        
        """Boundary and unit of this parameter"""
        self.boundary = (1000,1001)        
        self.unit     = 'uM' 
        
class Kflow(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = True,known = False)
        
        """Boundary and unit of this parameter"""
        self.boundary = (6.75/60.,350/60.)        
        self.unit     = 'ul'    
        
class Kvolume(K_):
    def __init__(self,name):
        """super to inherit the other properties"""
        super().__init__(name,control = False,known = True)
        
        """Boundary and unit of this parameter"""
        self.boundary = (250,250)        
        self.unit     = 'ul'    
             
class base_Model:
    def __init__(self,model,
                      equations,
                      parameters,
                      # """compilation of the model"""
                      compilation = True,
                      name        = 'Test_1'):
        
        """Get the model object"""
        from Model import ModelObject
        
        """The model as defined by the reaction kinetics"""
        self.kinetic_model     = model
        r,k = zip(*self.kinetic_model)
        kinetics = {i:[] for i in k}
        for i in range(len(k)):
            kinetics[k[i]].append(i)
        identifier = ''
        for i in sorted(kinetics.keys()):
            identifier     += i 
            reaction_number = ''
            for j in kinetics[i]:
               reaction_number += str(j) 
            identifier += '({})'.format(reaction_number)
            
        """define model name"""
        self.name = name

        """Get the equations that are concatenated together"""
        self.equations = equations
        
        """Get the parameters and the control parameters"""
        self.parameters         = parameters

        """States of the model"""
        self.states = []
        
        """The stringmodel of the network for OED"""
        self.stringmodel = ' '
        for state,fluxterms in self.equations.items():
            self.states.append(state)
            for flux in list(set(fluxterms)):
                self.stringmodel += ' ' + flux + ' '
            self.stringmodel += ' \n '
            
        print(self.stringmodel)
        """Remove the last jump"""
        self.stringmodel = self.stringmodel[0:-4]

        library = {}
        """Next we need to classify the parameters""" 
        
        for i in self.parameters:
            if 'k_cat' in i:
                library[i] = Kcat(i,)
            if 'k_fwd_' in i:
                library[i] = Kfwd(i)
            if 'k_1r_' in i:
                library[i] = Ka(i)
            if 'k_2r_' in i:
                library[i] = Ka(i)
            if 'k_1f_' in i:
                library[i] = Ka(i)
            if 'k_2f_' in i:    
                library[i] = Ka(i)
                
            if 'k_m_'   in i:
                library[i] = Km(i)
            if 'k_a_'   in i:
                library[i] = Ka(i)
            if 'k_b_'   in i:
                library[i] = Kb(i)
            if 'k_ai_'  in i:
                library[i] = Kai(i)
            if 'k_inh_'   in i:
                library[i] = Kinh(i)
            if 'k_act_'   in i:
                library[i] = Kact(i)
            if 'k_N_'   in i:
                library[i] = Kn(i)
            if 'kflow'   in i:
                library[i] = Kflow(i)
            if 'volume'   in i:
                library[i] = Kvolume(i)
            if '(in)'     in i:
                library[i] = Kstock(i)
                                      
        self.boundaries = {}
        for ID,value in library.items():
            self.boundaries[ID] = value.boundary

        import random
        self.fixed = {}
        for ID,bound in self.boundaries.items():
            """get hte upper and lower bound"""
            l,u = bound
            """Get centroid of log scaled range"""
            self.fixed[ID] = 10**((math.log10(l) + math.log10(u))/2)*random.choice([0.25,0.5,0.75,1])
            
        """Get the control parameters"""
        self.control = [ID for ID,v in library.items() if v.control]

        """Get the predefined values"""
        self.known   = [ID for ID,v in library.items() if v.known]
        
        """Define the included rates"""
        exclude      = self.known+self.control
        self.include = [ID for ID,v in library.items()  if v not in exclude]

        """Plug the model into the Model Builder"""
        self.model = ModelObject(self.stringmodel,
                                 self.states,
                                 self.boundaries,
                                 self.fixed,
                                
                                 #"""the name of the model"""
                                 name               = self.name,
                                 control_parameters = self.control,
                                 include            = self.include)
        
        """build an antimony file and SBML file"""
        self.model.AntimonyConversion()
        self.model.SBMLconversion()
        
        """compile the model to C++"""
        if compilation:
            """File create on desktop if it exists there will be a model file and used during simulation"""
            self.model.PytoCompile()     

    def return_model(self,fixed = {}):
        for name,value in fixed.items():        
            self.model.fixed[name] = value
        self.include = [i for i in self.include if i not in fixed.keys()]
        self.model.include = self.include
        return self.model