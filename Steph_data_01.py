'''
Created on May 18, 2014

@author: Keisuke
'''
from pylab import * 

class Dataset:    
    def __init__(self,plotcolor,markercolor,plotmarker,marker_size,V,Mu):
        self.plotcolor=plotcolor
        self.markercolor=markercolor
        self.plotmarker=plotmarker
        self.marker_size=marker_size
        self.V=V
        self.Mu=Mu

def dataset():    

    a=genfromtxt('Steph_data_01.csv', delimiter=',')
    nodata=zeros(a.shape[0])+nan
    #this nodata is for creating the column of nan
    
    plotcolor=['none','none','none','none','none']
    markercolor=['red','blue','cyan','#00FF00','black']
    plotmarker=['o','o','o','+','s']
    marker_size=[10,10,10,13,10]
  
    Diatom=Dataset(plotcolor[0],markercolor[0],plotmarker[0],marker_size[0],a[:,1],a[:,2])
    Calcifier=Dataset(plotcolor[1],markercolor[1],plotmarker[1],marker_size[1],a[:,4],a[:,5])
    Dino=Dataset(plotcolor[2],markercolor[2],plotmarker[2],marker_size[2],a[:,7],a[:,8])
    Small=Dataset(plotcolor[3],markercolor[3],plotmarker[3],marker_size[3],a[:,10],a[:,11])
    N2fixer=Dataset(plotcolor[4],markercolor[4],plotmarker[4],marker_size[4],a[:,13],a[:,14])

    return(Diatom,Calcifier,Dino,Small,N2fixer) 
