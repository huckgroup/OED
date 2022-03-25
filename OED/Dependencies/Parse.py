# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:53:31 2019

@author: huckg
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May 06 14:11:45 2019

@author: huckg
"""



import numpy 
import os
import csv

from collections import defaultdict
from DataTransform import *
from Scoring import *
from pyDOE import *


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