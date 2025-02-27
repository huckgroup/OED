# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:01:40 2022

@author: bob van sluijs
"""

"""parse this model in the model object"""
from Model import ModelObject

def main():
    models,control = {},{}
    
    bounds = { '(in)p70_S19':(0.1,3.5),
               '(in)p19_deGFP':(0.1,3.5),
               'KcatP70':(0.1,10),
               'degmRNAS19':(0.01,0.3),
               'KcatP19':(0.1,10),
               'kd_S19':(1,100),
               'degmRNAdeGFP':(0.01,0.3),    
               'KcatmRNAS19':(0.1,10), 
               'KcatmRNAdeGFP':(0.1,10),
               'kmatdeGFPdark':(0.099,0.1),
               'dil':(0.026,0.043),
               }   
    
    maximum_likelihood = { '(in)p70_S19':1,
                           '(in)p19_deGFP':1,
                           'KcatP70':5,
                           'degmRNAS19':0.025,
                           'KcatP19':5,
                           'kd_S19':100,
                           'degmRNAdeGFP':0.025,    
                           'KcatmRNAS19':0.5, 
                           'KcatmRNAdeGFP':0.5,
                           'kmatdeGFPdark':0.088,
                           'dil':0.026
                           }   
    
    modelname = 'Toy_Model'  
    stringmodel = """ + (in)p70_S19 * dil - p70_S19 * dil 
                      + (in)p19_deGFP * dil - p19_deGFP * dil 
                      + p70_S19 * KcatP70 - dil * mRNAS19 - degmRNAS19 * mRNAS19 
                      + p19_deGFP * KcatP19 * ( S19 / kd_S19 ) / ( 1 |+| ( S19 / kd_S19 ) ) - dil * mRNAdeGFP - degmRNAdeGFP * mRNAdeGFP 
                      + KcatmRNAS19 * mRNAS19 - dil * S19 
                      + KcatmRNAdeGFP * mRNAdeGFP - kmatdeGFPdark * deGFPdark - dil * deGFPdark 
                      + kmatdeGFPdark * deGFPdark - dil * deGFP """
                      
    conditions         = []
    states             = ['p70_S19','p19_deGFP','mRNAS19','mRNAdeGFP','S19','deGFPdark','deGFP']
    observables        = ['deGFP']
    control_parameters = ['(in)p70_S19','(in)p19_deGFP','dil']
    models[len(models)]     = ModelObject(stringmodel,states,bounds,maximum_likelihood,observables = observables,name = modelname,control_parameters = control_parameters)   
    return models
    

    