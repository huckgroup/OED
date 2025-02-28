# -*- coding: utf-8 -*-
"""
Created on Wed Feb 06 10:12:03 2019

@author: huckg
"""
import numpy
import copy
import math
from matplotlib.tri import Triangulation
from scipy.spatial import ConvexHull
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt 
from matplotlib import colors as mcolors
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import colors as mcolors
from Operations import *

def find_fixeddimension(triset,biset):
    p = [i for i in triset if i not in biset]
    return p[0] 

def build_array(indexinfo,dimensions,data):
    arrays  = {i:[] for i in dimensions}
    for i in indexinfo:
        for j in arrays.keys():
            arrays[j].append(data[j][i])
    return tuple([numpy.array(arrays[i]) for i in dimensions])

def heatmap(correlations,labels,criterion = 'Binary',store = False,show = True,path = ''):
    fig = plt.figure(figsize = (10,10))
    corr    = numpy.zeros((len(labels),len(labels)))
    indices = {labels[i]:i for i in range(len(labels))}
    for pair,value in correlations.items():
        row,column = pair
        corr[indices[row],indices[column]] = value
        corr[indices[column],indices[row]] = value
    mx = numpy.ma.masked_array(corr, mask=numpy.tri(corr.shape[0], k=-1))
    plt.imshow(mx, cmap='gist_earth', interpolation='nearest')
    plt.xticks([indices[i] for i in labels],labels,rotation = 90,size = 12)
    plt.yticks([indices[i] for i in labels],labels,size =12)
    plt.title("{}".format(criterion))
    plt.colorbar()
    plt.tight_layout()  
    if store:
        plt.savefig(path,dpi = 600)    
    if show:
        plt.show()  
    return 

def heatpanel(correlations,labels,criterion = 'Binary',ax = ''):
    corr  = numpy.zeros((len(labels),len(labels)))
    indices = {labels[i]:i for i in range(len(labels))}
    for pair,value in correlations.items():
        row,column = pair
        corr[indices[row],indices[column]] = value
        corr[indices[column],indices[row]] = value
    mx = numpy.ma.masked_array(corr, mask=numpy.tri(corr.shape[0], k=-1))
    im = ax.imshow(mx, cmap='gist_earth', interpolation='nearest')
    ax.set_xticks([indices[i] for i in labels])
    ax.set_yticks([indices[i] for i in labels])
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticklabels(labels)
    ax.set_title("Covariance: {}".format(criterion))   
    return ax,im

def barplot_1D(pearson,parameters,title = 'insert title',xl = 'insert label',yl = 'insert label',store = False,folder = '',path = '',ax = ''):
    xrn = range(len(parameters))
    sort_pearson = [0] * len(pearson)
    for i, x in enumerate(sorted(range(len(pearson)), key=lambda y: pearson[y])):
       sort_pearson[x] = i    
    tags = [None for i in range(len(parameters))]
    ordered_pearson = numpy.zeros((len(parameters)))
    
    cnt = 0
    for i in sort_pearson:
        tags[i] = parameters[cnt]
        ordered_pearson[i] = pearson[cnt]
        cnt += 1
    
    colors = []
    for i in ordered_pearson:
        if i > 0:
            colors.append('tab:blue')
        else:
            colors.append('tab:red')
    ax.bar(xrn, ordered_pearson, align='center',color = colors,edgecolor = 'k',alpha = 0.7)
    if yl != None:
        ax.set_ylabel(yl,fontsize=14)
    else:
        ax.set_ylabel(yl,fontsize=14)
    ax.set_title(title)
    ax.set_xticks(xrn)
    ax.set_xticklabels(tags,rotation=90)
    return ax



def twixplot(data,y_axis = '',x_axis = '',show = True):
    """"calculate the mean and median of the data array"""
    means   = [(index,numpy.mean(values)) for index,values in data.items()]
    medians  = [(index,numpy.median(values)) for index,values in data.items()]    
    """unzip the mean and medians """
    mean_idx,mean = zip(*means)
    median_idx,median = zip(*medians)
    """plot the figure """    
    fig = plt.figure(figsize =(10,10))
    ax = fig.add_subplot(111)
    ax.plot(mean_idx,mean,color = "DarkRed")
    ax.scatter(mean_idx,mean,color = "DarkRed")
    ax.set_yscale("log")
    ax2 = ax.twinx()
    ax2.plot(median_idx,median,color = "DarkBlue")  
    ax2.scatter(median_idx,median,color = "DarkBlue")  
    ax2.set_yscale("log")
    ax.legend(loc=0)
    ax.grid()
    ax.set_ylabel("Mean {}".format(y_axis),size = 14)
    ax2.set_ylabel("Median {}".format(y_axis),size = 14)
    ax.set_xlabel(x_axis)
    if show:
        plt.show()    

def analyse_distributions(parameters,scores,observables,criteria = ['mean','period','amplitude','oscillations'],store = False,folder = '',path = ''):   
    for parameter,information in scores.items():
        for criterion in criteria:
            executable = "{0} = {{i:[] for i in observables}}".format(criterion)
            exec(executable)
        values = []
        for i in information:
            scores,value = i
            values.append(value)
            for state in scores.keys():
                for k,v in scores[state].items():
                    eval(k)[state].append(v)
        clm = len(criteria)
        fig = plt.figure(figsize=(5*clm,5*len(observables)))  
        scr = ["DarkGreen","DarkBlue","DarkRed","black","teal"]
        row = 1; sc = 0;
        for state in scores.keys():
            for criterion in criteria:   
                output = eval(criterion)[state]
                ax = plt.subplot(len(criteria),clm,row)
                ax.plot(values,output,color = scr[sc])
                plt.ylabel(criterion + "\n" + state,color=scr[sc],fontsize = 14)
                plt.xlabel(parameter,fontsize = 12)  
                plt.tick_params(axis='both', which='minor', labelsize=14)
                row += 1 
            sc += 1
        fig.suptitle("Local Parameter Sensitivity: {}".format(parameter))
        plt.tight_layout()
        fig.subplots_adjust(top=0.88)
        if store:
            plt.savefig(path +  filecount(folder),dpi = 600)
        plt.show()
    return

def analyse_global_sensitivity(parameter,scores,observed,criteria = [],binary = [],store = False, folder = '', path = ''):
    for state in observed:
        if len(criteria) > 2 and len(criteria)%2 == 0:
            clm = len(criteria)/2
        else:
            clm = len(criteria)
        
        fig = plt.figure(figsize=(5*clm,8))  
        row = 1; sc = 0;
        for i in criteria:
            exec("{} = []".format(i))
        for i in range(len(scores)):
            for k,v in scores[i][state].items():
                exec("{0}.append({1})".format(k,v))    
        for i in criteria:
            if i not in binary:
                ax = plt.subplot(len(criteria),clm,row)
                binlength = max(eval(i)) - min(eval(i)) + 10
                if binlength >= 100:
                    binlength = 50
                ax.hist(eval(i),bins = int(binlength),color = 'tab:blue',edgecolor = 'k',alpha = 0.8)
                ax.set_title("{0} \n State: {1}".format(i,state))
                row += 1
        for i in criteria:
            if i in binary:
                ax = plt.subplot(len(criteria),clm,row)
                positive = [j for j in eval(i) if j == True]
                negative = [j for j in eval(i) if j == False]
                total    = len(positive) + len(negative)
                ax.bar([0.5,1.5], [100*float(len(positive))/float(total),100*float(len(negative))/float(total)], align='center',color = 'tab:blue',edgecolor = 'k',alpha = 0.8)
                ax.set_xlabel("   +                               -   ")   
                ax.set_ylabel("Percentage %")
                ax.set_title("{0} \n State: {1}".format(i,state))
                row += 1
        row = 0
        if store:
            plt.tight_layout()
            plt.savefig(path +  filecount(folder),dpi = 600)
        plt.show()
    return 

def analyse_phase_2D(triset,parameters,scores,observables,boundaries,criteria = [],binary = [],store = False,folder = '',path = ''):
    p_all,s_all,p_bin,s_bin = unpack_parameters(parameters,scores,observables,criteria = ['mean','period','amplitude','oscillations'],binary = binary)    
    p_pair = list(itertools.combinations(triset,2))
    
    for state in observables:
        clm = len(criteria) - len(binary)
        number = len(criteria) * len(p_pair)
        fig = plt.figure(figsize=(5*len(p_pair),5*len(criteria)))  
        row = 1; sc = 0;
        for biset in p_pair:
            for criterion in criteria:
                if criterion not in binary:
                    ax = plt.subplot(len(criteria),clm,row)
                    x,y = biset
                    p = find_fixeddimension(triset,biset) 
                    median = list(sorted(p_all[p]))[int(len(p_all[p])/2)]
                    indices=numpy.where(p_all[p] == median)[0]
                    x_arr,y_arr = build_array(indices,biset,p_all)
                    z_arr = build_array(indices,(criterion,),s_all[state])[0]
                    cols = numpy.unique(x_arr).shape[0]
                    
                    ax.set_title("Phase {0} \n {1}".format(state,criterion))
                    ax.tick_params(axis='both', which='minor', labelsize=12)
                    ax.set_xlabel(x,fontsize = 14)
                    ax.set_ylabel(y,fontsize = 14)
                    
                    X = x_arr.reshape(-1, cols)
                    Y = y_arr.reshape(-1, cols)
                    Z = z_arr.reshape(-1, cols)
                    cp = ax.contourf(X, Y, Z,cmap = 'gist_earth',antialiased = True)
                    fig.colorbar(cp,ax = ax)
                    row += 1
        plt.tight_layout()
        if store:
            plt.savefig(path + filecount(folder),dpi = 600)
        plt.show()      

def analyse_phase_3D(triset,parameters,scores,observables,boundaries,criteria = [],binary = [],store = False, folder ='',path = ''):
    p_all,s_all,p_bin,s_bin = unpack_parameters(parameters,scores,observables,criteria = criteria,binary = binary)    
    for state in observables:
        points = numpy.array([p_bin[state][i] for i in triset]).T
        cvx = ConvexHull(points)
        x, y, z = points.T
        
        tri = Triangulation(x, y, triangles=cvx.simplices)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(tri, z,cmap = "gist_earth")
        ax.set_xlabel(triset[0]); ax.set_ylabel(triset[1]); ax.set_zlabel(triset[2]);
        
        #set axis of the plot
        ax.set_xlim3d(boundaries[triset[0]])
        ax.set_ylim3d(boundaries[triset[1]])
        ax.set_zlim3d(boundaries[triset[2]]) 

        #compute centroid of the geometric structure           
        centroid = [numpy.mean(cvx.points[cvx.vertices,i]) for i in range(3)]   
        ax.tick_params(axis='both', which='minor', labelsize=12)
        ax.set_title("Phase {0} \n {1}: Volume = {2}".format(state,binary[0],str(round(cvx.volume,2))))
        plt.draw()  
        if store:
            plt.savefig(path +  filecount(folder),dpi = 600)
        plt.show()
    return

def Biplot(prs,obs,store = False, folder = '', path = ''):  
    fig = plt.figure(figsize=(5*len(obs),5))  
    clm,row = (len(obs),1)
    for o in obs:
        labels = prs[o].keys()
        dataset = numpy.array([prs[o][i] for i in labels])
        X = dataset.T
        
        """ We do PCA analysis and make a biplot of our data
        to assess to what extend our parameters are related to one another and how one
        might group them together for analysis"""
        scaler = StandardScaler()
        scaler.fit(X)
        X=scaler.transform(X)    
        pca = PCA()
        x_new = pca.fit_transform(X)
        ax = plt.subplot(len(obs),clm,row)
        def myplot(score,coeff,labels=None):
            """ plotting function to create a biplot that"""
            xs = score[:,0]
            ys = score[:,1]
            n = coeff.shape[0]
            scalex = 1.0/(xs.max() - xs.min())
            scaley = 1.0/(ys.max() - ys.min())
            ax.scatter(xs * scalex,ys * scaley,color = "DarkBlue", alpha = 0.2)
            labels = prs[o].keys()
            for i in range(n):
                ax.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
                ax.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')

        """ We show two dimensions of the biplot"""    
        plt.title("PCA analyisis {0}".format(o))
        ax.set_xlabel("PC{}".format(1))
        ax.set_ylabel("PC{}".format(2))
        ax.set_xlim([-1.5,1.5]);ax.set_ylim([-1.5,1.5])
        ax.grid()
        
        ''' Myplot is a function that creates the biplot'''
        myplot(x_new[:,0:2],numpy.transpose(pca.components_[0:2, :]))
        row += 1 
    plt.tight_layout()
    if store:
        plt.savefig(path + "phaseplot_" +  filecount(folder),dpi = 600)
    plt.show()
    return
    
    
def fill_plot(l_means,l_data,least_squares,percentage = 10,store = False, folder = '', path = ''):
    pct = int(100/percentage)
    """the least squares numbers"""
    fill_between = {i:[] for i in range(pct)}
    
    """per observable"""
    fill_observables = {i:copy.copy(fill_between) for i in l_means.keys()}
    
    """sort the data according to their lsq distance from the mean"""
    for observable,error in least_squares.items():
        ranks    = list(numpy.argsort(error))
        """"distance from the mean"""
        distance = 0        
        for i in range(len(ranks)):
            if i%len(ranks)/pct == 0 and i != 0:
                distance += 1          
            fill_observables[observable][distance].append(ranks.index(i))
            
    """"the level of the observables"""
    figure = plt.figure(figsize = (12,12))
    number = int(math.sqrt(len(l_means)))
    if number < 2 and len(l_means) > 1:
        number = 2
    if number == 0:
        number = 1

    """set plot window coordinates"""
    tln,clm,row = len(l_means),number,1    
    for observable,indices in fill_observables.items():
        ax = plt.subplot(tln,clm,row)
        
        means = []
        for levels in range(len(indices)):
            total = numpy.zeros(len(l_data[observable][-1]))
            for index in indices[levels]:
                total += l_data[observable][index]
            mean = total/float(len(indices[levels]))
            means.append(copy.copy(mean))
            

        for plotlevel in means:
            ax.fill_between(range(len(means[0])),means[0],plotlevel,alpha = 0.1,color = "tab:blue")      
        ax.plot(l_means[observable],color = 'k',label = observable)
        ax.legend(fancybox = True)
        ax.set_ylabel("Concentration",size = 14)
        ax.set_xlabel("Time",size = 14)
        ax.set_title("Likelihood Profile CI {}%".format(str(percentage)))
        row += 1
    plt.tight_layout()
    if store:
        plt.savefig(path + "likelihood_" +  filecount(folder),dpi = 600)
    plt.show()
            
    return