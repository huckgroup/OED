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
from PlotData import *
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
                 split   = 1,
                 nm = '',
                 #finite differences and forward problem
                 finite_difference = False,
                 forward_problem   = True):
        
        self.control = []
        """store simulated likelihoods"""
        self.simulated_likelihood = []
        """sample the likelihood space an simulate it, mxl as base first"""
        self.maximum_likelihoods  = []
        self.nm = nm
        
        """the measurements"""
        self.measurements = measurements
#       ______________________________________________________________________        
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
            for i in range(len(simdata)):
                for k,v in simdata[i].items():
                    pass
                        # plt.plot(v,label = k)
                        # plt.legend(fancybox = True)
            # # print(forward_sensitivities)
            # plt.gca().set_facecolor('whitesmoke')
            # plt.xlabel('Time (min)')
            # plt.ylabel('Concentration')
            # plt.savefig('C:\\Users\\huckg\\OneDrive\\Desktop\\Revision Glycolysis (Fig)\\'+nm+'.png',dpi = 600)
            # plt.close()
            # plt.show()
                        # 
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
                """padding of the array for experiments that do not contain"""
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
            t = {}
            for s in self.allstates:
                t[s] = numpy.array(self.statesensitivities[p][s]) * numpy.array(self.parametervectors[p])/numpy.array(self.statevectors[s])                           
                # for i in range(len(numpy.array(self.statevectors[s]))):
                    # if self.statevectors[s] < 2.5:
                        # self.statesensitivities[p][s][i] = 0
                # a = numpy.array(self.statesensitivities[p][s]) * (numpy.array(self.parametervectors[p])/numpy.array(self.statevectors[s]))
                t[s] = numpy.array([t[s][i] for i in range(1,len(t[s]),split)])
            self.corr_s[p] = copy.copy(t)  

    def Correlation(self,observables = True,normalized = True,removed = [],name = '',paper = False):
        import pandas
        removed += self.control
        """the sensitivities in the set"""
        vector = {p:pair for p,pair in self.corr_s.items() if p not in removed}
        """FIM"""
        self.FIM    = {}
        if observables == True:
            for p in list(vector.keys()):
                sx = []
                for o in self.observables:
                    for n in vector[p][o]:
                        if not math.isnan(n):
                            sx.append(n)
                            
                """Fisher information"""
                self.FIM[p] = copy.copy(sx)
                
        """DROP protential nan values in the right order rows then columns"""
        self.FIMname = {}
        for p in self.FIM.keys():
            try:
                string = p
                s = copy.copy(string)
                partial = string.split('_')[0] + '_' + string.split('_')[1] + '_' + string.split('_')[2]
                final = s.replace(partial,'')
                final = final.replace('_',' ')
                if len(final.split(' ')) > 2:
                    final = final.split(' ')[0] + final.split(' ')[1]  + ':' + final.split(' ')[2]
                string = '$'+ string.split('_')[0] +'_{' + string.split('_')[1] + "}" + string.split('_')[2] + '_{' + final + '}' + '$'
                self.FIMname[string] = self.FIM[p]
            except:
                pass

        """the figure clustermap"""
        points = pandas.DataFrame(self.FIM).corr()
        points = points.dropna(how = 'all')
        points = points.dropna(axis='columns', how ='all')
        
        """the figure clustermap"""
        plt.figure(figsize = (6,6),dpi = 650)
        sns.set(font_scale=1)
        fig = sns.clustermap(
           points,
           cmap="RdBu",
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
        plt.tight_layout()
        plt.savefig(directory + '\\correlation + {0}.svg'.format(name),dpi = 500)
        plt.savefig(directory + '\\correlation + {0}.png'.format(name),dpi = 500)  
        plt.show()
        if paper == True:
            directory = os.path.join(os.path.join(os.path.expanduser('~'))) +  '\\' + 'Onedrive' +  '\\' + 'desktop'  + '\\'  + 'Figure Data' + '\\'
            plt.savefig(directory + '\\SVD correlation + {0}.svg'.format(name),dpi = 500)    

    def SVD(self,observables = True,normalized = True,removed = [],name = '',paper = False):
        import pandas
        import math

        """the sensitivities in the set"""
        vector = {p:pair for p,pair in self.corr_s.items() if p not in removed and 'in' not in p}

        
        removed = ["HK",'G6PDH',"6PGDH",'PRI','PRPPS']
        enzymes = removed
        
        """FIM"""
        self.FIM    = {}
        if observables == True:
            for p in list(vector.keys()):
                if p not in enzymes:
                    sx = []
                    for o in self.observables:
                        for n in vector[p][o]:
                            if not math.isnan(n):
                                sx.append(n)
                                
                    """fisher information"""
                    self.FIM[p] = copy.copy(sx)

        collinearity = {i:{} for i in self.FIM.keys()}
        absolute     = {i:{} for i in self.FIM.keys()}
        
        # """calculate the svd of a pairwise sensititivity vector"""
        # for combination in itertools.combinations(list(self.FIM.keys()),2):
        #     try:
        #         array = []           
        #         for i in combination:
        #             if i not in removed:
        #                 if i in self.FIM.keys():
        #                     arr = numpy.array(self.FIM[i])
        #                     array.append(arr / numpy.linalg.norm(arr))            
        #         array = numpy.array(array).T
        #         arr,sv,uni = numpy.linalg.svd(array)
        #         try:
        #              c = numpy.log10((1./math.sqrt(float(abs(min(sv))))-1))
        #              print(c)
        #         except ZeroDivisionError:
        #              c = 0

        #         """create a full matrix"""
        #         if c != 0:
        #             collinearity[combination[0]][combination[1]] = c
        #             collinearity[combination[1]][combination[0]] = c
        #             absolute[combination[0]][combination[1]] = c 
        #             absolute[combination[1]][combination[0]] = c
        #         else:
        #             collinearity[combination[0]][combination[1]] =  0
        #             collinearity[combination[1]][combination[0]] =  0    
                    
        #         """autocorr to zero"""
        #         collinearity[combination[0]][combination[0]] = 0
        #         collinearity[combination[1]][combination[1]] = 0
        #     except:
        #         pass

        # """absolute values"""
        # absolute = pandas.DataFrame(absolute)
        # absolute = absolute.dropna(how = 'all')
        # absolute = absolute.dropna(axis='columns', how ='all')

        # """hte points in a dataframe?"""
        # points = pandas.DataFrame(collinearity)
        # points = points.dropna(how = 'all')
        # points = points.dropna(axis='columns', how ='all')
        
        array = [] 

        for i in list(self.FIM.keys()):       
            if i in ['k_ai_(6PGDH:6PGA_NADP)',
                      'k_ai_(PRPPS:ATP_R5P)',
                      'k_b_(6PGDH:NADPH)',
                      'k_ai_(G6PDH:6PGL_NADPH)',
                      'k_cat_(PRPPS:ATP_R5P)',
                      'k_cat_(HK:ATP_Glucose)',
                      'k_ai_(PRPPS:AMP_PRPP)',
                      'k_m_(HK:ATP)',
                      'k_b_(PRPPS:AMP)',
                      'k_cat_(G6PDH:6PGL_NADPH)',
                      'k_cat_(PRI:RU5P)',
                      'k_b_(6PGDH:6PGA)',
                      'k_b_(G6PDH:6PGL)',
                      'k_b_(G6PDH:G6P)',
                      'k_b_(HK:ADP)',
                      'k_b_(PRPPS:ATP)',
                      'k_m_(PRI:RU5P)',
                      'k_cat_(PRPPS:AMP_PRPP)',
                      'k_m_(PRI:R5P)',
                      'k_cat_(HK:ADP)',
                      'k_ai_(HK:ATP_Glucose)',
                      'k_cat_(HK:ADP_G6P)',
                      'k_m_(HK:ADP)',
                      'k_ai_(HK:ADP_G6P)',
                      'k_ai_(6PGDH:NADPH_RU5P)',
                      'k_cat_(6PGDH:6PGA_NADP)',
                      'k_cat_(PRI:R5P)',
                      'k_ai_(G6PDH:G6P_NADP)',
                      'k_cat_(HK:ATP)',
                      'k_b_(HK:ATP)',
                      'k_fwd_6PGA',
                      'k_fwd_6PGL',
                      'k_cat_(G6PDH:G6P_NADP)',
                      'k_cat_(6PGDH:NADPH_RU5P)']:
            
                arr = []
                for s in self.FIM[i]:
                    arr.append(round(s,9))

                array.append(numpy.array(arr))
        
        array = numpy.array(array)
        print(array,array.shape)
        FIS =  numpy.matmul(array,array.T)
        
        
        """check the boolean"""
        j = 0
        indices = []
        for i in FIS:
            # data.append(i )
            testbool = False
            for n in i:
                if n != 0:
                    testbool = True
            if testbool == False:
                indices.append(j)
            j += 1

    
        def remove_indices(arr, indices):
            import numpy as np
            # Convert the list of indices to a set for faster lookup
            index_set = set(indices)
        
            # Get the number of rows and columns in the original array
            num_rows, num_cols = arr.shape
        
            # Create a mask to exclude the specified rows and columns
            row_mask = np.array([i not in index_set for i in range(num_rows)])
            col_mask = np.array([i not in index_set for i in range(num_cols)])
        
            # Use the masks to create a new array with rows and columns removed
            new_arr = arr[row_mask][:, col_mask]
        
            return new_arr


        # FIS = remove_indices(FIS,indices)
        det = numpy.linalg.slogdet(FIS)

        arr,sv,uni = numpy.linalg.svd(array)
        try:
              c = 1./math.sqrt(float(abs(min(sv))))-1
              print(c)
        except ZeroDivisionError:
              c = 0
        
        plt.bar([0],[det[1]])
        print(c,det[1])
        plt.title(c)
        plt.savefig("C:\\Users\\huckg\\OneDrive\\Desktop\\Revision Glycolysis (Fig)\\barplot "+self.nm+'.png',dpi=600)
        plt.close()
        # import os
        # directory = os.path.join(os.path.join(os.path.expanduser('~'))) + '\\' + 'desktop'  + '\\' + 'Identification' + '\\' + name
        # try:
        #     os.mkdir(directory)
        # except:
        #     pass
        # plt.title('Total Coll: ' + str(round(totalcoll,1)))
        # plt.savefig(directory + '\\SVD collinearity + {0}.svg'.format(name),dpi = 750)    
        # plt.savefig(directory + '\\SVD collinearity + {0}.png'.format(name),dpi = 750)  
        # plt.show()
        # plt.close()
    
        # if paper == True:
        #     directory = os.path.join(os.path.join(os.path.expanduser('~'))) + '\\' + 'Onedrive' +  '\\' + 'desktop'  + '\\' + 'Figure Data' + '\\'
        #     plt.savefig(directory + '\\SVD collinearity + {0}.svg'.format(name),dpi = 500)    
        #     import pickle
        #     with open(directory + 'data {}.pickle'.format(name), 'wb') as handle:
        #         pickle.dump(absolute, handle, protocol=pickle.HIGHEST_PROTOCOL)        
                
            
    def PCA(self,name = '',paper = True):
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
        ks,inertias = range(1, len(parameters)),[]
        import os
        directory = os.path.join(os.path.join(os.path.expanduser('~'))) + '\\' + 'desktop'  + '\\' + 'Identification' + '\\' + name
        try:
            os.mkdir(directory)
        except:
            pass
        plt.savefig(directory + '\\PCA + {0}.png'.format(name),dpi = 500)
        if paper == True:
            directory = os.path.join(os.path.join(os.path.expanduser('~'))) +'\\' +'Onedrive' +  '\\' + 'desktop'  + '\\' + 'Figure Data' + '\\'
            plt.savefig(directory + '\\SVD PCA + {0}.svg'.format(name),dpi = 500)    
            
    def analyse(self,observables = [],combination_size = 2,show = True,store = False,directory = "",name = '',absolute = False,removed = []):
        """the fixed parameters not part of the sensitivity"""
        fixed = {}
        for i in range(len(self.measurements)):
            for k,v in self.measurements[i].model.fixed.items():
                fixed[k] = v 

        """"remove the control paramters from the equation"""
        self.sensitivity_vector = {k:v for k,v in self.sensitivity_vector.items() if k not in list(set(self.control))}
        self.sensitivity_vector = {k:v for k,v in self.sensitivity_vector.items() if fixed[k] != 0}
        self.sensitivity_vector = {k:v for k,v in self.sensitivity_vector.items() if k not in removed}
        # print(self.sensitivity_vector.keys())
        ######################################################
        """add observables to the model"""
        if not observables:
            self.observables = self.observables      
        self.osv = {}  
        """observed collinearity"""
        for parameter,observed_sensitivity in self.sensitivity_vector.items():
            """stack the columns into a single column vector, first extract the relevant columns"""
            sensitivity = []
            for i in self.observables:
                sensitivity.append(observed_sensitivity[i])
            """create new sensitivity vector with concatenated observables"""
            self.osv[parameter] = numpy.concatenate(sensitivity)
        
            
        # self.observed_sensitivities,cumsense = single_value_decomposition(self.osv,self.observables,pair = True)
        ######################################################
        """absolute sensitivity"""
        self.absolute_sensitivities = {}
        """absolute collinearity"""
        for observable in self.allstates:
            absolute_sensitivity = {}
            for parameter,sensitivity in self.sensitivity_vector.items():
                """stack the columns into a single column vector, first extract the relevant columns"""
                absolute_sensitivity[parameter] = sensitivity[observable] 
            self.absolute_sensitivities[observable],cumsense = single_value_decomposition(absolute_sensitivity,observable, pair = True)
            
        """means of the colls"""
        means,medians = {},{}      
        for observable,c_obj in self.absolute_sensitivities.items():
            """"calculate the means and rankorder"""
            means[observable] = numpy.mean([i.collinearity for i in c_obj if i])
            medians[observable] = numpy.median([i.collinearity for i in c_obj if i])
        """unpack and rankorder"""
        obs,value = zip(*means.items())
        ordered = [i for i in numpy.argsort(value)]
     
        """store information addition"""
        self.seqinfadd = {}
        """forward and reversed introduction of observables"""
        ranked_observables = []
        for i in range(len(ordered)):
            ranked_observables.append(obs[ordered.index(i)])
            """sensitivity"""
            added_information = {}
            for parameter,observed_sensitivity in self.sensitivity_vector.items():
                """stack the columns into a single column vector, first extract the relevant columns"""
                sensitivity = []
                for j in ranked_observables:
                    sensitivity.append(observed_sensitivity[j])
                """create new sensitivity vector with concatenated observables"""
                added_information[parameter] = numpy.concatenate(sensitivity)      
            self.seqinfadd[tuple(ranked_observables)],cumsense = single_value_decomposition(added_information,ranked_observables)  

        """store information addition"""
        self.revseqinfadd = {}
        """forward and reversed introduction of observables"""
        ranked_observables = []
        for i in reversed(range(len(ordered))):
            ranked_observables.append(obs[ordered.index(i)])
            """sensitivity"""
            added_information = {}
            for parameter,observed_sensitivity in self.sensitivity_vector.items():
                """stack the columns into a single column vector, first extract the relevant columns"""
                sensitivity = []
                for j in ranked_observables:
                    sensitivity.append(observed_sensitivity[j])
                """create new sensitivity vector with concatenated observables"""
                added_information[parameter] = numpy.concatenate(sensitivity)      
            self.revseqinfadd[tuple(ranked_observables)],cumsense = single_value_decomposition(added_information,ranked_observables)  

        """Pairwise clustering between the parameters"""
        IDs = []
        for i in self.observed_sensitivities:
            for j in i.parameters:
                IDs.append(j)
        IDs = list(set(IDs))
        parameters   = {IDs[i]:i for i in range(len(IDs))}
        index        = {i:IDs[i] for i in range(len(IDs))}
        
        """collinearity between pairs"""
        distance_grid = numpy.zeros((len(parameters),len(parameters)))
        
        """allcoll"""
        coll = []
        
        """pairs of the parameters"""
        for i in self.observed_sensitivities:
            r,c = i.parameters
            
            """collinearities"""
            distance_grid[parameters[r],parameters[c]] = i.collinearity
            distance_grid[parameters[c],parameters[r]] = i.collinearity
            
            """windsorizing the data"""
            coll.append(i.collinearity)
        sort = coll 
        if max(coll) > 10**6:
            sort = [sorted(coll)[i] for i in range(int(0.9*len(coll)))] 
        maxsort = max(sort)
        for i in range(len(distance_grid)):
            for j in range(len(distance_grid)):
                if distance_grid[i,j] not in sort and distance_grid[i,j] != 0:
                    distance_grid[i,j] = 2*maxsort
        
        # Fails with "ValueError: similarities must be symmetric"
        prn,pcoll,pmean,pmedian,pci = {},{},[],[],[]
        for k,v in self.absolute_sensitivities.items():
            prn[k] = len([i.collinearity for i in v if i.collinearity])
            """collinaearity per parameter"""
            for i in v:
                if i.collinearity:
                    """loop through the parameter"""
                    for j in i.parameters:
                        if j not in pcoll.keys():
                            pcoll[j] = []
                        pcoll[j].append(i.collinearity)
                        
        """find mean, median and indices"""
        for p,v in pcoll.items():
                pmean.append(numpy.mean(v))
                pmedian.append(numpy.median(v))
                pci.append(p)
        indices = numpy.argsort(pmean)
        
        """"create the identification directory"""
        if store:
            if not os.path.exists(directory):
                os.makedirs(directory)  
                         
        parameter = {} 
        """rank calculation"""
        for p,obs in self.sensitivity_vector.items():
            vector = numpy.empty(1)
            for i in sorted(obs.keys()):
                vector = numpy.append(vector,[obs[i]])
                
            parameter[p] = copy.copy(vector.T)
        rankmatrix = numpy.array([v for k,v in parameter.items()])
        rank = numpy.linalg.matrix_rank(rankmatrix)
        
        """"parameter linkage, to what extend do parameters cluster together based 
        on the overall collinearity within the system i.e. you can find 90% of the 
        sensitivity within the parameters or you can find the observables best
        suited to find parameters the distance if proportional to the dissimilarity the bigger the larger 
        the distance we plot the dendogram and the correlation matrix in the same graph
        so see potential clusters and double check"""
        if store == True or show == True:
            sns.set_style("whitegrid")
            plt.figure(figsize = (14,7))
            ax = plt.subplot(1,2,1)
            dissimilarity = copy.copy(distance_grid)
            for i in range(len(distance_grid)):
                for j in range(len(distance_grid)):
                    if distance_grid[i,j] != 0:
                        v = 1./distance_grid[i,j]
                        if math.isnan(v):
                            dissimilarity[i,j] = 1
                        else:
                            dissimilarity[i,j] = v 
                            
            dst = squareform(dissimilarity)
            linkage_matrix = linkage(dst, "single")
            dendrogram(linkage_matrix, labels= [index[i] for i in range(len(parameters))],orientation = "right")
            ax = plt.subplot(1,2,2)
            plt.imshow(distance_grid, zorder=2, cmap='Blues', interpolation='nearest')
            plt.xticks(range(len(parameters)), [index[i] for i in range(len(parameters))],rotation=90)
            plt.yticks(range(len(parameters)), [index[i] for i in range(len(parameters))])
            plt.colorbar()
            plt.title("Pairwise Parameter Distance: Rank {}".format(str(rank)))
            plt.tight_layout()
            if store:
                plt.savefig(directory + "\\Dendogram {0} {1}.png".format(name,name),dpi = 600)
            if show:
                plt.show()
            plt.close()
                
            from sklearn.preprocessing import StandardScaler
            import pandas as pd
            """3 plots that further doublecheck the clusters that were formed by doing
            a PCA analysis, in this instance we do the PCA, we subsequently plot the
            data of the scatter after which plot the kmeans cluster inertia reflexting
            the number of clusters there should be, more corroborating evidence"""
            try:
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
                ax.set_xlabel('PCA features',size = 14)
                ax.set_ylabel('variance %',size = 14)
                plt.xticks(features)        
                ax = plt.subplot(1,3,2)
                PCA_components = pd.DataFrame(principalComponents)
                ax.scatter(PCA_components[0], PCA_components[1], alpha=.1, color='black',s = 30)
                ax.set_xlabel('PCA 1',size = 14)
                ax.set_ylabel('PCA 2',size = 14)
                ax = plt.subplot(1,3,3)
                ks,inertias = range(1, len(parameters)),[]
                for k in ks:
                    # Create a KMeans instance with k clusters: model
                    model = KMeans(n_clusters=k)
                    # Fit model to samples
                    model.fit(PCA_components.iloc[:,:3])
                    # Append the inertia to the list of inertias
                    inertias.append(model.inertia_)
                ax.plot(ks, inertias, '-o', color='black')
                ax.set_xlabel('number of clusters, k',size = 14)
                ax.set_ylabel('inertia',size = 14)
                plt.xticks(ks)
                plt.tight_layout()
                if store:
                    plt.savefig(directory + "\\PCA {0} {1}.png".format(name,self.model.name),dpi = 600)
                if show:
                    plt.show()
                plt.close()
            except:
                pass
                
            #get indices of sorted means, plot information in observables in barplot
            observed  =  [i for i in obs]
            sorted_observables = [observed[i] for i in numpy.argsort(value)]
            figure = plt.figure(figsize = (10,10))
            sns.set_style("whitegrid")  
            ax = figure.add_subplot(111)
            ax.bar(range(len(value)),sorted(value),color = "tab:blue",edgecolor = 'k')
            ax.set_ylabel("Mean Collinearity",size = 14,color = "tab:blue")
            ax.set_yscale("log")
            ax2 = ax.twinx()
            ax2.scatter(range(len(sorted_observables)),[prn[i] for i in sorted_observables],color = "DarkRed")
            ax2.set_ylabel("Parameter Number",size = 14,color = "DarkRed")
            ax2.set_ylim([0,max([prn[i] for i in sorted_observables])])
            plt.xticks(range(len(sorted_observables)),sorted_observables,rotation = 90)
            if store:
                plt.savefig(directory + "\\Information {0} {1}.png".format(name,name),dpi = 600)
            if show:
                plt.show()
            plt.close()                
            #parameter collinearity in the observable set
            figure = plt.figure(figsize = (10,10))
            plt.bar(range(len(pmean)),[pmean[i] for i in indices],width=0.4,color = "tab:blue",label = "Mean",edgecolor = 'k')
            plt.bar(numpy.array(range(len(pmean)))+0.4,[pmedian[i] for i in indices],width=0.4,color = "tab:red",label = "Median",edgecolor = 'k')
            plt.xticks(numpy.array(range(len(pmean)))+0.2,[pci[i] for i in indices],rotation = 90,size = 14)
            plt.legend(fancybox = True)
            plt.ylabel("Collinearity",size = 14)
            plt.yscale("log")
            if store:
                plt.savefig(directory + "\\Collinearity {0} {1}.png".format(name,name),dpi = 600)
            if show:
                plt.show()
            plt.close()            
            #Collinearity in terms of worst to best observable################# 
            data = {len(i):[] for i in self.seqinfadd.keys()}
            for k,v in self.seqinfadd.items():
                for i in v:
                    if i.collinearity!= 0:
                        data[len(k)].append(i.collinearity)
            sns.set_style("whitegrid")
            twixplot(data,y_axis = "Collinearity Index",x_axis = "Observable Number",show = False) 
            if store:
                plt.savefig(directory + "\\Collinearity LowHigh {0} {1}.png".format(name,name),dpi = 600)
            if show:
                plt.show()
            plt.close()                
            #Same but the reverse, best to worst###############################            
            data = {len(i):[] for i in self.revseqinfadd.keys()}
            for k,v in self.revseqinfadd.items():
                for i in v:
                    if i.collinearity!= 0:
                        data[len(k)].append(i.collinearity)
            sns.set_style("whitegrid")
            twixplot(data,y_axis = "Collinearity Index",x_axis = "Observable Number",show = False)
            if store:
                plt.savefig(directory + "\\Collinearity Highlow {0} {1}.png".format(name,name),dpi = 600)
            if show:
                plt.show()
            plt.close()   
            
