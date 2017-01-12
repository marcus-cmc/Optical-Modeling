# -*- coding: utf-8 -*-
"""
@author: C. Marcus Chuang
"""

"""
Instruction:
You can use this file to enter user inputs and then run optical modeling,
or treat this file (all the default values) as an example of how to run it.

To run your own simulation, you should have a file of the library of refraction
indices of the materials in the device you want to model and then replace those
variables in the 'mandatary user input' section to run.

"""

from TransferMatrix import OpticalModeling
import matplotlib.pyplot as plt
plt.style.use('ggplot') # or use 'classic' or any in plt.style.available

#------------------------ User input ------------------------------------

#------------------- Mandatary user input -------------------------------
"""
Device: (Name, thciness) ; thickness in nm
*** Names of layers of materials must match that in the library ***
For example, if the layer is called "Glass",
the library should contain a columns called "Glass_n" and "Glass_k"
where n, k are the real and imaginary part of the refraction index.

Device stack starts from the side where light is incident from.
The first layer is assumed to be a thick substrate whose thickness is
irrelivant. If a thin substrate is used, add "Air" as the first layer
indice are 0-based, i.e, the index of the first layer is 0
"""
Device = [
          ("Glass"  , 100), # layer 0
          ("ITO"    , 145), # layer 1
          ("ZnO"    , 120),
          ("PbS"    , 250),
          ("Au"     , 150)
         ]

"""
file name of the  refraction index library, must be in the same folder
The first row in the file is the header, for a material named "Mater1",
there sould be a "Mater1_n" column, and a "Mater1_k" column.
The header for the wavelength column ***MUST be*** : "Wavelength (nm)"
"""

libname = "Index_of_Refraction_library_Demo.csv"

# file name of the AM1.5 solar spectra, already in the folder
# change it ONLY if you want to use your own file
Solarfile = "SolarAM15.csv" # Wavelength vs  mW*cm-2*nm-1

wavelength_range = [350, 1200] # wavelength range (nm) to model [min, max]

#------------------- End of Mandatary input -------------------------------




# Below are the optioanl input ------

posstep = 0.5 # step size for thickness, must be smaller than the thinnest layer
WLstep = 2.0 # wavelength step size (nm)

# selected wavelengths for "E vs position plot",
# can be several values or a single value,
# if set to None or not provided, it would automatically calculate 1 or 3 values
# with gap of a multiple of 50 nm
plotWL = [450, 600, 700, 950]

plotE = True   # whehter to plot E-field vs wavelength
plotAbs = True # whether to plot absorption vs wavelength
plotGen = True # whether to plot generation rate and spectral absorption rate

SaveName = "Result" # prefix of the file names
# whether to save the data as csv files
saveDataE, saveDataAbs, saveDataGen = False, False, False
# wherther to save the figures
# default format is vector graphic 'pdf' (with non-transparent background)
# can also use 'png', 'jpg' or someother format matplotlib supports
figformat = 'pdf'
saveFigE , saveFigAbs , saveFigGen  = False, False, False

#------------------- End of user input -------------------------------


if __name__=="__main__":

    #### --------------- run with Mandatary input only --------
    ## initialize an OpticalModeling obj OM
    #OM = OpticalModeling(Device, libname=libname, WLrange=wavelength_range)
    ## do all the caculation
    #OM.RunSim()
    #### ---------------------------------------------------------------------


    #### --------------- run with all the user input  --------
    # initialize an OpticalModeling obj OM
    OM = OpticalModeling(Device, libname=libname, WLrange=wavelength_range,
                         plotWL=plotWL, WLstep=WLstep, posstep=posstep)
    OM.RunSim(plotE=plotE, plotAbs=plotAbs, plotGen=plotGen,
              saveFigE=saveFigE, saveFigAbs=saveFigAbs, saveFigGen=saveFigGen,
              figformat='pdf', savename=SaveName)
    ####----------------------------------------------------------------------

    plt.show()

    summary = OM.JscReport() # print and return a summary report
    ### save data as csv or not
    OM.SaveData(savename=SaveName,
                saveE=saveDataE, saveAbs=saveDataAbs, saveGen=saveDataGen)

    """
    Note:
    After calling OM.RunSim(), the simulation would have been done, and all the
    data can be re-accessed and figures can be re-plotted without running the
    simulation again.

    To get some data, access these attributes:
        OM.Absorption  : wavelength vs absorption for each material
        OM.E   : electric field (row: position; column: wavelength)
        OM.Gx  : generation rate (row: position; column: wavelength)
        OM.Jsc : Jsc for each layer (not including the substrate)

    You can re-plot the figures again and decide whether to save them,
    what name and fileformat of the figures you want to save by calling the
    following methods and replace the default variables with what you prefer:

        OM.PlotAbs(savename="Result", savefig=False, figformat='pdf')
        OM.PlotE(  savename="Result", savefig=False, figformat='pdf')
        OM.PlotGen(savename="Result", savefig=False, figformat='pdf')

    To save the simulation data as .csv files, call the OM.SaveData() method and
    replace the default with what you prefer:
         OM.SaveData(savename="Result",
                     saveAbs=True, saveE=False, saveGen=False)
    """


