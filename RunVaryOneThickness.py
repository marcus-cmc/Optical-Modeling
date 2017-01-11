# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from VaryThickness import OMVaryThickness

plt.style.use('ggplot')


# Device: (Name, thciness) ; thickness in nm
# Names of layers of materials must match that in the library
# ex: layer name: Glass ; library should contains Glass_n & Glass_k
# starting from the side where light is incident from
# first layer is assumed to be a thick substrate whose thickness is irrelivant
# if a thin substrate is used, add "Air" as the first layer


#------------------------------ User input --------------------------------


Device = [
          ("Glass"  ,   0), # layer 0
          ("ITO"    , 145), # layer 1
          ("ZnO"    , 120), # layer 2
          ("PbS"    , 220),
          ("Au"     , 100)
         ]


##############  vary the thickness of one layer     ##############
#VaryOneLayer = False # vary the thickness of one layer or two layers(False)
ToVary = 2 # the layer to vary
#t_range = np.arange(100, 601, 10) # start, end (not included), step
t_range = np.arange(20, 301, 10)
#t_range = [50, 75, 125, 150, 250, 300, 350] # manually input range

# target: layer of interest (layer index), usually the light absorber.
# Will calculate the max Jsc in this layer (assuming 100% IQE)
# alternatively, can use 'R' for reflection or 'T' for transmission
target = 3
#target = 'R'


##############  vary the thickness of two layers     ##############
#VaryTwoLayer = not VaryOneLayer # vary the thickness of two layers
#
#ToVary2= 1
#t2_range = np.arange(50, 300, 50)
#target2 = None # for tandem only, calculate and plot the Jsc of the tandem
#               # cell with absorber target1 and target 2 (min of these)
#               # i.e. the current limiting case. Use None for non-tandem device
#########################################################################


libname = "Index_of_Refraction_library_Demo.csv"
Solarfile = "SolarAM15.csv" # Wavelength vs  mW*cm-2*nm-1

posstep = 1.0 # thickness step size
WLrange = [350, 1200] # wavelength range (nm)
WLstep = 2.0 # wavelength step size (nm)

SaveName = "Result"

cbarlegend = True # colorbar as legend for the "thickness vs Abs" plot
                   # but if there are >25 curves, colorbar will be used

interp_countour = True # True : Contour plot (interplate data)
                        # False: Heatmap (no interp)
#########################################################################


if __name__=="__main__":
    VT = OMVaryThickness(Device, libname=libname, WLrange=WLrange,
                         plotWL = None, WLstep = WLstep, posstep = posstep,
                         Solarfile = Solarfile)

    VT.VaryOne(ToVary, t_range, target, toPrint=False,
               cbarlegend=cbarlegend)
    pass

#    if VaryTwoLayer:
#        VT.VaryTwo(L1=ToVary, t1_range = t_range,
#                   L2=ToVary2, t2_range = t2_range,
#                   target1 = target, target2=target2, toPlot=True,
#                   print1=True, print2=False, interp_countour=interp_countour)

