# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 20:14:42 2015

@author: Marcus
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import pandas as pd
from TransferMatrix import OpticalModeling as TM
#import time
#from scipy import signal
plt.style.use('ggplot')

# Device: (Name, thciness) ; thickness in nm 
# Names of layers of materials must match that in the library
# ex: layer name: Glass ; library should contains Glass_n & Glass_k
# starting from the side where light is incident from
# first layer is assumed to be a thick substrate whose thickness is irrelivant
# if a thin substrate is used, add "Air" as the first layer



#------------------------------ User input --------------------------------
#Device = [
#          ("Glass"  , 500), # layer 0
#          ("ITO"    , 145), # layer 1
#          ("ZnO"    ,  70), # layer 2
#          ("PbSTBAI", 300), 
#          ("PbSEDT" ,  45),
#          ("Au"     , 150),
##          ("ZnO"    ,  70),
##          ("PbSTBAI", 200),
##          ("PbSEDT" ,  45),
##          ("Au"     , 150)
#         ]


Device = [
          ("Glass"  , 500), # layer 0
          ("ITO"    , 145), # layer 1
          ("ZnO"    , 120), # layer 2
          ("PbSTBAI", 220),
          ("PbSEDT" , 45),
          ("Au"     , 150)
         ]


##############  vary the thickness of one layer     ##############
VaryOneLayer = True # vary the thickness of one layer or two layers(False)
ToVary = 3 # the layer to vary
#t_range = np.arange(100, 601, 10) # start, end (not included), step
t_range = np.arange(50, 501, 50)
#t_range = [ 50,75,125,150,250,300,350] # manually input range

# target: layer of interest (layer index), usually the light absorber. 
# Will calculate the max Jsc in this layer (assuming 100% IQE)
# alternatively, can use 'R' for reflection or 'T' for transmission
target = 3
#target = 'R'


##############  vary the thickness of two layers     ##############
VaryTwoLayer = not VaryOneLayer # vary the thickness of two layers

ToVary2= 2
t2_range = np.arange(20, 201, 10)
target2 = None # for tandem only, calculate and plot the Jsc of the tandem
               # cell with absorber target1 and target 2 (min of these)
               # i.e. the current limiting case. Use None for non-tandem device
#########################################################################


libname = "Index_of_Refraction_library.csv"
Solarfile = "SolarAM15.csv" # Wavelength vs  mW*cm-2*nm-1

posstep = 1.0 # thickness step size
wavelength = (300, 1200) # wavelength range (nm)
WLstep = 1.0 # wavelength step size (nm)

SaveName = "Result"

cbarlegend = True # colorbar as legend for the "thickness vs Abs" plot
                   # but if there are >25 curves, colorbar will be used

interp_countour = True # True : Contour plot (interplate data)
                        # False: Heatmap (no interp)
# ----------------------------END of user input---------------------------

class VaryThickness(TM):

    def set_t_update(self, layerind, newthickness):
        """
        set the thickness of the layer "layerind" and then update
        relevant parameters
        """
        self.t[layerind] = newthickness
        self.t_cumsum = np.cumsum(self.t)
        self.x_pos = np.arange(self.WLstep/2.0, sum(self.t) , self.posstep)
        self.x_ind = self.x_indice()
        return None

    def VaryOne(self, L_vary, t_range, target, toPrint=False, 
                PlotJsc=True, PlotAbs=True, cbarlegend=False):
        """
        vary the thickness of the layer with index L_vary and then run optical 
        modeling. 
        L_vary: index of the layer to vary
        t_range: thickness to vary, an iterable
        self.varyJsc: a list of "Jsc in each layer" w.r.t the varying thickness
        self.varyAbsorption: a list of "absorption in each layer"
        return None
        """
        self.varyJsc = []
        self.varyAbs = []
        self.varyT = []
        self.varyR = []
        self.L_vary = L_vary
        self.t_range = t_range
        if not isinstance(target, str):
            target = int(target) # in case user input is a float
        else:
            target = target.upper() # in case user input is a lowercase letter
            
        for ti in t_range:
            self.set_t_update(L_vary, ti) # update t & related variables
            self.RunSim(plotE=False, plotAbs=False, plotGen=False)
            self.varyJsc.append(self.Jsc)
            self.varyAbs.append(self.Absorption)
            self.varyT.append(self.Transmission)
            self.varyR.append(self.Reflection)
            if toPrint:
                print "calculating: ", self.layers[L_vary], "=", ti, "nm,"
                #print self.layers[L_vary], "=", ti, "nm,", "Max Jsc in", 
                #print self.layers[target], "=", np.round(self.Jsc[target-1],2)
        if PlotJsc and not isinstance(target, str):
            self.PlotVaryJsc(target)
        if PlotAbs:
            self.PlotVaryAbs(target, cbarlegend)
        return None

    def PlotVaryJsc(self, target):
        ftitle = 'Max Jsc in L' + str(target) + ' ' + self.layers[target]
        figJsc = plt.figure(ftitle)
        figJsc.clf()
        axJsc = figJsc.add_subplot(111)
        xlabel = 'Thickness of ' + self.layers[self.L_vary] + ' (nm)'
        ylabel = 'Jsc' + " (mA/cm$^2$)"
        axJsc.set_xlabel(xlabel, size=16)
        axJsc.set_ylabel(ylabel, size=16)
        axJsc.plot(self.t_range, [Jsc[target-1] for Jsc in self.varyJsc], 
                   '-o', linewidth=2, color='r', markersize=8)
        axJsc.tick_params(labelsize=14)
        figJsc.suptitle(ftitle, fontsize=14)
        return None
        
    def PlotVaryAbs(self, target, cbarlegend=False):

        #figAbs = plt.figure('absorption', figsize=(16*0.8, 9*0.8))
        figAbs = plt.figure(figsize=(16*0.8, 9*0.8))
        figAbs.clf()
        axAbs = figAbs.add_subplot(111)
        axAbs.set_xlabel('Wavelength (nm)', size=24)
        
        cmap = plt.get_cmap('rainbow')
        num_color = len(self.t_range)
        vmin, vmax = t_range[0], t_range[-1]
        normalize = mcolors.Normalize(vmin, vmax)

        if target=='R':
            ftitle = 'Reflection'
            axAbs.set_ylabel('Reflection (%)', size=24)
            for ind, t in enumerate(self.t_range):
                axAbs.plot(self.WL, 100.0*self.varyR[ind], 
                           linewidth=2, label = str(t) + " nm",
                           color=cmap(normalize(t)))

        elif target=='A':
            ftitle = '1-R-T'
            axAbs.set_ylabel('1-R-T (%)', size=24)
            for ind, t in enumerate(self.t_range):
                axAbs.plot(self.WL, 100.0*(1-self.varyR[ind]-self.varyT[ind]), 
                           linewidth=2, label = str(t) + " nm",
                           color=cmap(normalize(t)))

            
        elif target=='T':
            ftitle = 'Transmission'
            axAbs.set_ylabel('Transmission (%)', size=24)  
            for ind, t in enumerate(self.t_range):
                axAbs.plot(self.WL, 100.0*self.varyT[ind], 
                           linewidth=2, label = str(t) + " nm",
                           color=cmap(normalize(t)))
        else:
            #ftitle = 'Modeled Absorption in '+self.layers[target]
            ftitle = 'Absorption in L'+str(target)+' '+self.layers[target]
            targethead = self.Absorption.columns[target-1] 
            axAbs.set_ylabel('Modeled Absorption (%)', size=24)
        
            for ind, t in enumerate(self.t_range):
                axAbs.plot(self.WL, 100.0*self.varyAbs[ind][targethead], 
                           linewidth=2, label = str(t) + " nm",
                           color=cmap(normalize(t)))
                       
        axAbs.tick_params(labelsize=18)
        figAbs.suptitle(ftitle, fontsize=14)
        
        # use normal legend
        if num_color <= 20 and not cbarlegend: 
            axAbs.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), 
                         numpoints=1, fontsize=14,
                         title='Thickness of\n ' + self.layers[self.L_vary],
                         borderaxespad=0).draggable()
            figAbs.subplots_adjust(right=0.82)
        # use colorbar legend when specified or >20 lines 
        else:
            scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=cmap)
            scalarmappaple.set_array(t_range)
            cblegend = plt.colorbar(scalarmappaple)
            cblabel = 'Thickness of ' + \
                       self.layers[self.L_vary] + " (nm)"

            cblegend.set_label(cblabel, fontsize=14, labelpad=16)
            cblegend.ax.tick_params(labelsize=14)

        return None


    def VaryTwo(self, L1, t1_range, L2, t2_range, target1, target2=None, 
                toPlot=True, print1=False, print2=False, 
                interp_countour=False):
        """
        vary the thickness of the two layers L1 and L2. 
        L1, L2: indice of the layers to vary
        t1_range, t2_range: thickness to vary, an iterable
        target1: the layer used to calculate Jsc (@100%IQE)
        target2: Only used for tandem cells, where target1 and target2
                 are the two absorber. Calculate the min of (Jsc1, Jsc2),
                 i.e. the current limiting Jsc.
                 Default is None (single junction cells)
        The result will be stored in self.v2Jsc, a 3-D np array
        
        return None
        """
        v2Jsc=[]
        self.L1, self.t1range = L1, t1_range
        self.L2, self.t2range = L2, t2_range
        for t1 in self.t1range:
            self.set_t_update(self.L1, t1)
            if print1:
                print "Calculating: ", self.layers[L1], "=", t1, "nm,",
                print "Varying ", self.layers[L2]
            self.VaryOne(self.L2, self.t2range, target1, toPrint=print2,
                         PlotJsc=False, PlotAbs=False)
            v2Jsc.append(self.varyJsc)
        self.v2Jsc = np.array(v2Jsc)
        if toPlot:
            self.PlotTwo(target1, target2, interp_countour)
        return None

    def PlotTwo(self, target1, target2=None, interp_contour=False):
        
        figV2 = plt.figure()
        figV2.clf()
        axV2  = figV2.add_subplot(111)
        xlabel = 'Thickness of L' + str(self.L1)+ ' ' + self.layers[self.L1]
        ylabel = 'Thickness of L' + str(self.L2)+ ' ' + self.layers[self.L2]
        axV2.set_xlabel(xlabel + ' (nm)', size=14)
        axV2.set_ylabel(ylabel + ' (nm)', size=14)
        X, Y = np.meshgrid(self.t1range, self.t2range)

        v2J1 = self.v2Jsc.take(target1-1, axis=2)
        if target2!=None: # for tandem
            v2J2 = self.v2Jsc.take(target2-1, axis=2)
            p, q = v2J1.shape[0], v2J1.shape[1]
            J = np.array([ [min(v2J1[i][j], v2J2[i][j]) for j in range(q)] 
                               for i in range(p)])
            V2title = ('Max Jsc in the device\n' +
                       '(min of L'+ str(target1) + ' ' + self.layers[target1]+
                       ' and ' +
                       'L'+ str(target2) + ' ' + self.layers[target1] +")" )
        else:
            J = v2J1.T
            V2title = ('Max Jsc in L' + str(target1) + 
                        ' ' + self.layers[target1])                          
        ## heat map, no interpolation
        if not interp_contour:
            CS=axV2.pcolormesh(X, Y, J)
            axV2.set_xlim(self.t1range[0], self.t1range[-1])
            axV2.set_ylim(self.t2range[0], self.t2range[-1])
        else: ## contourf, interpolate data
            axV2.contourf(X, Y, J, 100, lw=0.1)
            CS = axV2.contourf(X, Y, J, 100)
        axV2.tick_params(labelsize=14)
        figV2.colorbar(CS)
        figV2.suptitle(V2title, fontsize=14)       
        
        return None

## To do:
    def SaveVary(SaveName):
        return 


if __name__=="__main__":

    VT = VaryThickness(Device, libname=libname, wavelength=wavelength,
                       WLstep = WLstep, posstep = posstep,
                       Solarfile = Solarfile)
    if VaryOneLayer:
        VT.VaryOne(ToVary, t_range, target, toPrint=True,cbarlegend=cbarlegend)

    if VaryTwoLayer:
        VT.VaryTwo(L1=ToVary, t1_range = t_range, 
                   L2=ToVary2, t2_range = t2_range, 
                   target1 = target, target2=target2, toPlot=True,
                   print1=True, print2=False, interp_countour=interp_countour)
                   

