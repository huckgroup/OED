# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:43:15 2021

@author: Bob
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:18:27 2019

@author: huckg
"""
from Parse import *
from Operations import *
from OptimizationOperations import *

import libsbml
import importlib
import amici
import os
import sys
import numpy as numpy
import copy
import seaborn as sns  
from matplotlib.pyplot import cm  
import matplotlib as mpl    
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance 
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy

from Scoring import __score__
from Model import ModelVariables
from SolveSystem import ModelSolver,AmiciDataObject,PythonDataObject

class IdentifiabilityAnalysis:
    def __init__(self,measurements,
                 #conditions and actual observables of the model
                 conditions   = [],
                 observables  = [],
                 include      = [],
                 time_dependent_parameters = {},
                 parameters                = {},
                 #simulation time
                 normalization = False,
                 show = True,
                 samples = 10,
                 time    = (0,1000),
                 dt      = 1,
                 #finite differences and forward problem
                 finite_difference = False,
                 forward_problem   = True):
        
        self.control = []
        """store simulated likelihoods"""
        self.simulated_likelihood = []
        """sample the likelihood space an simulate it, mxl as base first"""
        self.maximum_likelihoods  = []
        
        """the measurements"""
        self.measurements = measurements
#       _______________________________________________________________________        
        """loop thriugh the observables i.e. what needs to be plotted"""        
        self.observables = []
        for i in measurements.values():
            for j in i.observables:
                self.observables.append(j)
                
#       _______________________________________________________________________     
        """least squares error of the mean (maximum likelihood)"""
        self.lsq_sensitivities = {};
        """"now simulate the local devations from the respective means"""
        self.sensitivity_vector = {}
        """fixed parameters"""
        self.values = {}
        for i in range(len(measurements)):
            for ID,value in measurements[i].model.fixed.items():
                self.values[ID] = value
                if parameters != {} and ID in parameters:
                    self.values[ID] = parameters[ID]
                    measurements[i].model.fixed[ID] = parameters[ID]
                    
        sensitivity = {}
        for i in range(len(measurements)):
            for n in measurements[i].model.parameters:
                sensitivity[n] = []
            for n in measurements[i].model.control:
                self.control.append(n)
        self.control = list(set(self.control))   
            
        self.allstates = []
        for i in range(len(measurements)):
            for j in measurements[i].model.states:
                self.allstates.append(j)
        self.allstates     = list(set(self.allstates))
        
        self.allparameters = []
        for i in range(len(measurements)):
            for k,v in measurements[i].model.fixed.items():
                self.allparameters.append(k)
        self.allparameters = list(set(self.allparameters))
        
        """collection of all statevectors and all the parameters"""
        self.statevectors       = {i:[] for i in self.allstates}
        self.parametervectors   = {i:[] for i in self.allparameters}
        self.statesensitivities = {i:{j:[] for j in self.allstates} for i in self.allparameters}
        
        if forward_problem:        
            """"solve the system"""
            scores,simdata,forward_sensitivities = simulate_measurements(measurements,parameters = self.values,forward = True)

            """add states to statevector"""
            for i in range(len(simdata)):
                s = []
                for state,data in simdata[i].items():
                    s.append(state)
                    for point in data:
                        self.statevectors[state].append(point)
                        
                """padding of the array for experiments that do not contain observable or parameter with a zero in the dataset"""
                for ms in self.allstates:
                    if ms not in s:
                        for n in range(len(data)):
                            self.statevectors[ms].append(1)
                
                px = []
                for pID,pair in forward_sensitivities[i].items():
                     sx = []; px.append(pID)
                     for state,vector in pair.items():
                        sx.append(state)
                        for n in vector:
                            self.statesensitivities[pID][state].append(n)
                     for si in self.allstates:
                         if si not in sx:
                             for n in range(len(vector)):
                                     self.statesensitivities[pID][si].append(0)

                for p in self.allparameters:
                    if p not in px:
                        for s in self.allstates:
                            for n in range(len(vector)):
                                self.statesensitivities[p][s].append(0)
                      
            for p in self.allparameters:
                for n in range(len(self.statevectors[ms])):
                    self.parametervectors[p].append(self.values[p])
                    
        self.corr_s = {}
        for p in self.allparameters:
            t    = {}
            for s in self.allstates:
                norm = numpy.array(self.statesensitivities[p][s])  * (numpy.array(self.parametervectors[p])/numpy.array(self.statevectors[s]))
                t[s] = numpy.array([norm[i] for i in range(1,len(norm),1)])
            self.corr_s[p] = copy.copy(t) 
        self.observables = list(set(self.observables))
           
    def Correlation(self,observables = True,normalized = True,removed = [],name = '',additional_observation = []):
        import pandas
        removed += self.control

        """the sensitivities in the set"""
        vector = {p:pair for p,pair in self.corr_s.items() if p not in removed}
        """FIM"""
        self.FIM    = {}
        if observables == True:
            for p in list(vector.keys()):
                sx = []
                for o in self.observables + additional_observation:
                    for n in vector[p][o]:
                        # print(n)
                        if not math.isnan(n):
                            sx.append(n)
                """fisher information"""
                self.FIM[p] = copy.copy(sx)

    
        """DROP protential nan values in the right order rows then columns"""
        points = pandas.DataFrame(self.FIM).corr()
        points = points.dropna(how = 'all')
        points = points.dropna(axis='columns', how ='all')
        
        """the figure clustermap"""
        plt.figure(figsize = (10,10),dpi = 750)
        sns.set(font_scale=0.6)
        fig = sns.clustermap(
           points,
           cmap="vlag",
           metric='correlation', 
           linewidths=0.5,
           center=0,
           vmin=-1,
           vmax=1,
           xticklabels=True,
           yticklabels=True,
           annot_kws={"size": 6})

        plt.title("Correlations")  

        import os
        directory = os.path.join(os.path.join(os.path.expanduser('~'))) + '\\' + 'desktop'  + '\\' + 'Identification' + '\\' + name
        try:
            os.mkdir(directory)
        except:
            pass
        plt.savefig(directory + '\\correlation + {0}.png'.format(name),dpi = 500)
        plt.close()
 

    def SVD(self,observables = True,normalized = True,removed = [],name = '',additional_observation = []):
        import pandas
        import math
        removed += self.control

        """the sensitivities in the set"""
        vector = {p:pair for p,pair in self.corr_s.items() if p not in removed and 'in' not in p}
        """FIM"""
        self.FIM    = {}
        if observables == True:
            for p in list(vector.keys()):
                sx = []
                for o in self.observables + additional_observation:
                    for n in vector[p][o]:
                        if not math.isnan(n):
                            sx.append(n)
                """fisher information"""
                self.FIM[p] = copy.copy(sx)

        collinearity = {i:{} for i in self.FIM.keys()}
        absolute     = {i:{} for i in self.FIM.keys()}
        """calculate the svd of a pairwise sensititivity vector"""
        for combination in itertools.combinations(list(self.FIM.keys()),2):
            try:
                array = []           
                for i in combination:
                    # if i in include:
                        if i in self.FIM.keys():
                            arr = numpy.array(self.FIM[i])
                            array.append(arr / numpy.linalg.norm(arr))            
                array = numpy.array(array).T
                arr,sv,uni = numpy.linalg.svd(array)
                try:
                     c = 1./math.sqrt(float(abs(min(sv))))
                except ZeroDivisionError:
                     c = 0
                if math.isnan(c):
                     c = 10**-5
                if c == 0:
                    c = 10**-5
                if math.isinf(c):
                     c = 10**5       

                """create a full matrix"""
                if c != 0:
                    collinearity[combination[0]][combination[1]] = c
                    collinearity[combination[1]][combination[0]] = c
                    absolute[combination[0]][combination[1]] = c 
                    absolute[combination[1]][combination[0]] = c
                else:
                    collinearity[combination[0]][combination[1]] =  0
                    collinearity[combination[1]][combination[0]] =  0    
                    
                """autocorr to zero"""
                collinearity[combination[0]][combination[0]] = 0
                collinearity[combination[1]][combination[1]] = 0
            except:
                pass

        """absolute values"""
        absolute = pandas.DataFrame(absolute)
        absolute = absolute.dropna(how = 'all')
        absolute = absolute.dropna(axis='columns', how ='all')


        """hte points in a dataframe?"""
        points = pandas.DataFrame(collinearity)
        points = points.dropna(how = 'all')
        points = points.dropna(axis='columns', how ='all')
        
        """the figure clustermap"""
        plt.figure(figsize = (10,10),dpi = 750)
        sns.set(font_scale=0.6)  
        fig = sns.clustermap(
           points,
           metric='correlation', 
           linewidths=0.5,
           center=1,
           xticklabels=True,
           yticklabels=True,
           cmap="coolwarm",
            vmin = 0,
            vmax = 5,
           annot_kws={"size": 10})
        plt.title("Collinearities")          
        
        import os
        directory = os.path.join(os.path.join(os.path.expanduser('~'))) + '\\' + 'desktop'  + '\\' + 'Identification' + '\\' + name
        try:
            os.mkdir(directory)
        except:
            pass
        plt.savefig(directory + '\\SVD collinearity + {0}.png'.format(name),dpi = 500)  
        plt.close()
     
    def PCA(self,name = ''):
        from sklearn.preprocessing import StandardScaler
        import pandas as pd
        """3 plots that further doublecheck the clusters that were formed by doing
        a PCA analysis, in this instance we do the PCA, we subsequently plot the
        data of the scatter after which plot the kmeans cluster inertia reflexting
        the number of clusters there should be, more corroborating evidence"""
        plt.figure(figsize = (21,7))
        sns.set_style("whitegrid")
        ax = plt.subplot(1,3,1)
        parameters = list(self.FIM.keys())
        PCA_dataframe = pd.DataFrame(data = self.FIM)   
        pca = PCA(n_components=len(parameters))
        X_std = StandardScaler().fit_transform(PCA_dataframe)
        X_std += numpy.min(X_std)
        principalComponents = pca.fit_transform(X_std)
        features = range(pca.n_components_)
        ax.bar(features, pca.explained_variance_ratio_,color = 'black')
        ax.set_xlabel('PCA features',size = 7)
        ax.set_ylabel('variance %',size = 14)
        plt.xticks(features)        
        ax = plt.subplot(1,3,2)
        PCA_components = pd.DataFrame(principalComponents)
        ax.scatter(PCA_components[0], PCA_components[1], alpha=.1, color='black',s = 30)
        ax.set_xlabel('PCA 1',size = 14)
        ax.set_ylabel('PCA 2',size = 14)
        ax = plt.subplot(1,3,3)
        
        import os
        directory = os.path.join(os.path.join(os.path.expanduser('~'))) + '\\' + 'desktop'  + '\\' + 'Identification' + '\\' + name
        try:
            os.mkdir(directory)
        except:
            pass
        plt.savefig(directory + '\\PCA + {0}.png'.format(name),dpi = 500)
