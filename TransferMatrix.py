# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 12:40:50 2015

@author: Marcus
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
 

#------------------------------ User input --------------------------------


# Device: (Name, thciness) ; thickness in nm 
# Names of layers of materials must match that in the library
# ex: if the layer is called Glass, 
# the library should contains Glass_n and Glass_k
# starting from the side where light is incident from
# The first layer is assumed to be a thick substrate whose thickness is 
# irrelivant, if a thin substrate is used, add "Air" as the first layer

# ----- Mandatary input ------
Device = [
          ("Glass"  , 1000), 
          ("ITO"    , 145), 
          ("ZnO"    , 120), 
          ("PbSTBAI", 220), 
          ("PbSEDT" ,  0),
          ("Au"     , 150)
         ]

libname = "Index_of_Refraction_library.csv"
Solarfile = "SolarAM15.csv" # Wavelength vs  mW*cm-2*nm-1

wavelength = (300, 1200) # wavelength range (nm) to model
plotWL = [450, 600, 700, 950] # selective wavelengths for "E vs position" 



# ----- Optional Input -----

plotE = True   # plot E-field vs wavelength
plotAbs = True # plot absorption vs wavelength
plotGen = True # plot generation rate and spectral absorption rate
saveDataE, saveDataAbs, saveDataGen = False, True, False
saveFigE , saveFigAbs , saveFigGen  = False, True, False
SaveName = "Result"

posstep = 1.0 # step size for thickness
WLstep = 1.0 # wavelength step size (nm)


# ---------------------------- End of  user input ----------------------------

# Constants
h = 6.626e-34 # Js Planck's constant
c = 2.998e8   # m/s speed of light
q = 1.602e-19 # C electric charge

class OpticalModeling(object):
    def __init__(self, Device=Device, libname=libname, wavelength=wavelength, 
                 WLstep  = 2.0, posstep = 0.5, Solarfile = "SolarAM15.csv"):
        # layer (materials)
        self.layers = [Device[0][0]] + [mat[0] for mat in Device[1:] 
                       if float(mat[1]) > 0]
        # thickness; set the thickness of the substrate (first layer) to 0
        self.t = [0] + [mat[1] for mat in Device[1:] if float(mat[1]) > 0]
        self.t_cumsum = np.cumsum(self.t)
        # wavelength
        self.WL = np.arange(wavelength[0], wavelength[1] + WLstep, WLstep)
        self.WLstep = WLstep
        self.posstep = posstep
        #positions to evaluate field
        self.x_pos = np.arange(self.WLstep/2.0, sum(self.t) , self.posstep)
        # material i is in x_pos[x_ind[i-1]:xind[i]]
        #self.x_ind = self.x_indice(self.x_pos, self.t, self.t_cumsum)
        self.x_ind = self.x_indice()
        self.AM15 = self.LoadSolar(Solarfile) # load/reshape into desired WL range
        self.nk = self.Load_nk(libname)
        
        self.E = None
        self.Reflection = None
        self.Transmission = None
        self.AbsRate = None
        self.Absorption = None
        self.Gx = None
        self.Jsc = None   

    def RunSim(self, plotE=True, plotAbs=True, plotGen=True,
                     saveFigE=False, saveFigAbs=False, saveFigGen=False):
        self.CalE()
        self.CalAbs()
        self.CalGen()
        if plotE:
            self.PlotE(savefig=saveFigE)
        if plotAbs:
            self.PlotAbs(savefig=saveFigAbs)
        if plotGen:
            self.PlotGen(savefig=saveFigGen)
        return None

    def SaveData(self, savename="Resulttest_",
                 saveDataE=False, saveDataAbs=False, saveDataGen=False):

        if saveDataE:
            df_E = pd.DataFrame(abs(self.E**2))
            df_E.index = self.x_pos
            df_E.to_csv(savename+"_E.csv", header=self.WL) 
        if saveDataGen:
            df_Gen = pd.DataFrame(self.Gx)
            df_Gen.index = self.x_pos
            df_Gen.to_csv(savename+"_Gen.csv", header=self.WL)
        if saveDataAbs:
            df_Abs = pd.DataFrame(self.Absorption)
            df_Abs["Transmission"] = self.Transmission
            df_Abs["Reflection"] = self.Reflection
            df_Abs.index = self.WL
            df_Abs.to_csv(savename+"_Absorption.csv")

        return None
        
        
    def LoadSolar(self, Solarfile):
        Solar = pd.read_csv(Solarfile, header = 0) 
        AM15 = np.interp(self.WL, Solar.iloc[:,0], Solar.iloc[:,1])
        return AM15 # mW/cm2 nm

    def Load_nk(self, libname):
        # load data
        data=pd.read_csv(libname, header = 0)
        # initialize nk, a dataframe for complex indices of each layer
        # the column names in "nk" are the same as in "layers"  
        nk = pd.DataFrame()
        nk["WL"] = self.WL
        nk["Air"] = [1]*len(self.WL) 
        # interp n,k to desired WL range and store them in nk 
        for mater in self.layers:
            if mater not in nk:
                n=np.interp(self.WL, data["Wavelength (nm)"], data[mater+"_n"])
                k=np.interp(self.WL, data["Wavelength (nm)"], data[mater+"_k"])
                nk[mater] = [complex(n[i],k[i]) for i in xrange(len(self.WL))]
        return nk

    def x_indice(self):
        """
        return a list of indice "x_ind" for use in x_pos
        material i corresponds to the indices range  
        [x_ind[i-1], xind[i]) in x_pos
        Note: the first layer is glass and is excluded in x_pos 
        """
        x_ind = [len(self.x_pos) for _ in xrange(len(self.t)) ]
        j = 0
        for i, x in enumerate(self.x_pos):
            if x > self.t_cumsum[j]:
                x_ind[j] = i
                j += 1
        return x_ind         
        
    def CalE(self):
        # Calculate Incoherent power transmission through substrate
        # T = |4*n1*n2 / (n1+n2)^2| , R = |((n1-n2)/(n1+n2))^2|
        T_glass = abs(4*1*self.nk[self.layers[0]] / 
                      (1+self.nk[self.layers[0]])**2 )
        R_glass = abs(((1-self.nk[self.layers[0]])/
                      (1+self.nk[self.layers[0]]))**2 )
        
        # Calculate transfer matrices and field at each wavelength and position
        
        self.E = np.zeros([len(self.x_pos), len(self.WL)], dtype=complex)
        R = self.WL*0.0
        T = self.WL*0.0

        Imats, Lmats = {}, {}
        for matind in xrange(len(self.layers)-1):
            mater, nex = self.layers[matind], self.layers[matind+1]
            Lmats[matind] = self.L_mat(matind)
            if (mater, nex) not in Imats:
                Imats[(mater, nex)] = self.I_mat(mater, nex)
        Imats[(nex, "Air")] =  self.I_mat(nex, "Air")
        Lmats[len(self.layers)-1] = self.L_mat(len(self.layers)-1)
        
        for i, w in enumerate(self.WL):
            # Calculate transfer matrices for incoherent reflection and 
            # transmission at the first interface
            S = Imats[(self.layers[0], self.layers[1])][i]
            for L in xrange(1, len(self.layers)-1):
                mater, nex = self.layers[L], self.layers[L+1] 
                S = S.dot(Lmats[L][i])
                S = S.dot(Imats[(mater,nex)][i])

            S = S.dot(Lmats[len(self.layers)-1][i])
            S = S.dot(Imats[(nex,'Air')][i])

            # JAP Vol 86 p.487 Eq 9 
            # Power Reflection from layers other than substrate
            R[i] = (abs(S[1,0]/S[0,0]))**2 
            
            #Transmission of field through glass substrate 
            #Griffiths 9.85 + multiple reflection geometric series
            T[i] = ( abs(2.0/(1 + self.nk[self.layers[0]][i])) / 
                    (1-R_glass[i]*R[i] )**0.5  )


            # Calculate all other transfer matrices
            for L in xrange(1, len(self.layers)):
                mater = self.layers[L]
                xi = 2.0*np.pi*self.nk[mater][i] / w
                dj = self.t[L]
                # x: distance from previous layer
                x = (self.x_pos[self.x_ind[L-1]:self.x_ind[L]] - 
                     self.t_cumsum[L-1] ) 
        
                # Calculate S matrices (JAP Vol 86 p.487 Eq 12 and 13)
                S_prime = Imats[(self.layers[0], self.layers[1])][i]
                for matind in xrange(2, L+1):
                    mater, pre = self.layers[matind], self.layers[matind-1]
                    S_prime = S_prime.dot(Lmats[matind-1][i])
                    S_prime = S_prime.dot(Imats[(pre, mater)][i])

                S_dprime = np.eye(2, dtype=complex)                 

                for matind in xrange(L, len(self.layers)):
                    if matind < len(self.layers)-1:
                        mater, nex = self.layers[matind], self.layers[matind+1]
                        S_dprime = S_dprime.dot(Imats[(mater, nex)][i])
                        S_dprime = S_dprime.dot(Lmats[matind+1][i])
                    else:
                        mater, nex = self.layers[matind], "Air"
                        S_dprime = S_dprime.dot(Imats[(mater, nex)][i])
                        
                # Normalized Field profile (JAP Vol 86 p.487 Eq 22)
                numerator = (S_dprime[0,0] * np.exp(complex(0,-1.0)*xi*(dj-x))+
                             S_dprime[1,0] * np.exp(complex(0, 1.0)*xi*(dj-x)))                   
                denom = (S_prime[0,0]*S_dprime[0,0]*
                         np.exp(complex(0,-1.0)*xi*dj) ) 
                denom += (S_prime[0,1] * S_dprime[1,0]* 
                          np.exp(complex(0, 1.0)*xi*dj) ) 
                        
                self.E[self.x_ind[L-1]:self.x_ind[L], i]=T[i]*numerator/denom        
        
        # Overall Reflection from device with incoherent reflections 
        # at first interface
        self.Reflection = R_glass + T_glass**2*R / (1-R_glass*R)
        
        return None


    def CalAbs(self):
        """
        Calculate normalized intensity absorbed /cm3-nm at each position and
        wavelength as well as the total reflection expected from the device
        """
        # Absorption coefficient in cm^-1 (JAP Vol 86 p.487 Eq 23)
        a = pd.DataFrame()
        for mater in self.layers[1:]:
            if mater not in a:
                a[mater] = 4*np.pi*self.nk[mater].imag/(self.WL*1e-7)
                
        # initialize Absrate with E^2, multiply nk later
        self.AbsRate = abs(self.E)**2 
        self.Absorption = pd.DataFrame() # initialize Absorption (% of light)
        for matind in xrange(1, len(self.t)):
            mater = self.layers[matind]
            posind = self.x_ind[matind-1], self.x_ind[matind]
            mlabel = "L" + str(matind) + "_" + mater
            self.AbsRate[posind[0]:posind[1]]*=a[mater]*np.real(self.nk[mater])
            self.Absorption[mlabel] = (np.sum(self.AbsRate[
                                              posind[0]:posind[1]],0)
                                       * self.posstep * 1e-7 )        
        self.Transmission = 1 - np.sum(self.Absorption,1) - self.Reflection

        return None


    def CalGen(self):
        """
        Calculate generation rate as a function of position in the device 
        and calculates Jsc (in mA/cm^2)
        """    
        # Energy dissipation mW/cm3-nm at each position and wavelength 
        # (JAP Vol 86 p.487 Eq 22)  
        if self.AbsRate is None:
            self.CalAbs()
        Q = self.AbsRate * self.AM15
        self.Gx = Q*1e-12/(h*c)*self.WL
    
        Gx_x = [np.sum(self.Gx[self.x_ind[i-1]:self.x_ind[i]]) 
                for i in xrange(1,len(self.layers))]
        self.Jsc = np.array(Gx_x) * self.WLstep * self.posstep * q * 1e-4
        return None
        
    def PlotE(self, plotWL=plotWL, savename=SaveName+"_Efiled", savefig=False):
        """
        Plots electric field intensity |E|^2 vs position in device for
        wavelengths specified in the initial array, plotWavelengths. 
        """
            
        # E-field, selected wavelength
        fig1 = plt.figure(1)
        plt.clf()
        ax1 = fig1.add_subplot(111)
         
        ax1.set_ylabel('Normalized |E|$^2$Intensity', size=20)
        ax1.set_xlabel('Position in Device (nm)', size=20)
        ax1.tick_params(labelsize=18)
     
        E2 = abs(self.E**2)
        for i, w in enumerate(plotWL):
            label = "%s nm" % w
            #xind = min(enumerate(WL), key= lambda x: abs(x[1]-w))[0]
            # find the index closest to the desired wavelength
            xind = min(xrange(len(self.WL)), key= lambda x: abs(self.WL[x]-w))
            #E2 = abs(E[:, xind])**2
            ax1.plot(self.x_pos, E2[:,xind], label=label, linewidth=2)
        ax1.set_ylim(ymin=0)
         
        # E-field, contour
        fig2 = plt.figure("E-filed")
        plt.clf()
        ax2  = fig2.add_subplot(111)
        ax2.set_ylabel('Wavelength (nm)', size=20)
        ax2.set_xlabel('Position (nm)', size=20)
        X, Y = np.meshgrid(self.x_pos, self.WL)
        #ax2.contourf(X,Y,E2.T, 50, lw=0.1)
        CS = ax2.contourf(X,Y,E2.T, 50)
        for c in CS.collections: # avoid white gaps when converting to pdf
            c.set_edgecolor("face")
        ax2.tick_params(labelsize=18)
        fig2.colorbar(CS)
        fig2.suptitle('Normalized E-field Intensity', 
                       fontsize=20)
    
        # layer bar
        for matind in xrange(2, len(self.layers)+1):
            ax1.axvline(self.t_cumsum[matind-1], color="black")
            ax2.axvline(self.t_cumsum[matind-1], color="black")
            x_text = (self.t_cumsum[matind-2] + self.t_cumsum[matind-1])/2.0
            ax1.text(x_text, ax1.get_ylim()[1]+0.01, self.layers[matind-1], 
                     size=16, va="bottom", ha="center")
            ax2.text(x_text, ax2.get_ylim()[1]+0.01, self.layers[matind-1], 
                     size=14, va="bottom", ha="center")                 
        ax1.set_xlim(0, max(self.x_pos))
        ax1.legend(loc='upper right', fontsize=18).draggable()
        #plt.tight_layout()
        if savefig:
            fig1.savefig(savename+"_selectedWL.pdf", transparent=True)
            fig2.savefig(savename+".pdf", transparent=True)            
        return None
 
        
    def PlotAbs(self, savename=SaveName+"_Absorption", savefig=False):
        """
        Plots normalized intensity absorbed /cm3-nm at each position and
        wavelength as well as the reflection expected from the device
        """
    
        fig3 = plt.figure("Absorption")
        plt.clf()
        ax3  = fig3.add_subplot(111)
        ax3.set_ylabel('Fraction of Light', size=20)
        ax3.set_xlabel('Wavelength (nm)', size=20)
        ax3.tick_params(labelsize=18)
        
        for matind in xrange(1, len(self.t)):
            mater = self.layers[matind]
            mlabel = "L" + str(matind) + "_" + mater
            ax3.plot(self.WL, self.Absorption[mlabel], 
                     label=mlabel, linewidth=2)
        ax3.plot(self.WL, self.Transmission, label="Transmission", linewidth=2)    
        ax3.plot(self.WL, self.Reflection, label="Reflection", linewidth=2)
        ax3.legend(loc='upper right', fontsize=14).draggable()
        ax3.set_ylim(ymax=1.0,ymin=0)
        fig3.show()

        if savefig:
            fig3.savefig(savename+".pdf", transparent=True)
            
        return None

    
    def PlotGen(self, savename=SaveName, savefig=False):
        """
        Plots generation rate as a function of position in the device
        """    
  
        Gx_pos = np.sum(self.Gx,1)
        fig4 = plt.figure("Generation Rate")
        fig4.clf()
        ax4  = fig4.add_subplot(111)
        ax4.set_xlabel('Position (nm)', size=20)
        ax4.set_ylabel('Generation Rate (1/sec$\cdot$cm$^3$)', 
                       size=20)
        ax4.plot(self.x_pos, Gx_pos, linewidth=2, color="r")
        
        fig5 = plt.figure("Absorption Rate")
        fig5.clf()
        ax5  = fig5.add_subplot(111)
        ax5.set_ylabel('Wavelength (nm)', size=20)
        ax5.set_xlabel('Position (nm)', size=20)
        X, Y = np.meshgrid(self.x_pos, self.WL)
        
        ax4.tick_params(labelsize=18)
        ax5.tick_params(labelsize=18)
        #ax5.contourf(X, Y, self.Gx.T, 50) 
        CS = ax5.contourf(X, Y, self.Gx.T, 50, vmax=6e19)
        for c in CS.collections: # avoid white gaps when converting to pdf
            c.set_edgecolor("face")
        ax5.tick_params(labelsize=18)
        fig5.colorbar(CS)
        fig5.suptitle('Photon Absorption Rate (1/sec$\cdot$nm$\cdot$cm$^3$)', 
                       fontsize=16, fontweight='bold')
    
        for matind in xrange(2, len(self.layers)+1):
            ax4.axvline(self.t_cumsum[matind-1], color="black")
            ax5.axvline(self.t_cumsum[matind-1], color="black")
            x_text = (self.t_cumsum[matind-2] + self.t_cumsum[matind-1])/2.0
            ax4.text(x_text, ax4.get_ylim()[1]+0.01, self.layers[matind-1], 
                     size=14, va="bottom", ha="center")
            ax5.text(x_text, ax5.get_ylim()[1]+0.01, self.layers[matind-1], 
                     size=14, va="bottom", ha="center")
        if savefig:
            fig4.savefig(savename+"_Gen_position_.pdf", transparent=True)
            fig5.savefig(savename+"_AbsorptionRate.pdf", transparent=True)
        return None

    def JscReport(self):
        '''
        OM: Optical Modeling object
        print Jsc report: layer/materials/thickness/Jsc_max(100% IQE)
        return Jsc report (a pd dataframe)
        '''
        self.JscData = pd.DataFrame()
        self.JscData["Layer No."] = range(1,len(self.layers))
        self.JscData["Material"] = self.layers[1:]
        self.JscData["Thickness (nm)"] = self.t[1:]
        self.JscData["Jsc_Max (mA/cm^2)"] = np.round(self.Jsc, 2)
        print self.JscData.to_string(index=False)
        return None
    
    def I_mat(self, mat1, mat2):
        """
        Calculate the transfer matrix I for Reflection and Transmission
        at an interface between two materials.
        mat1, mat2: name of the materials
        return I, a  numpy array with shape len(self.WL)x2x2 
        I[i] is the transfer matrix at wavelength self.WL[i]
        
        I[i] = 1/T[i] *  [ [    1,  R[i] ]
                           [ R[i],     1 ] ]
        """
        n1s, n2s = self.nk[mat1], self.nk[mat2] # complex dielectric constants
        R=(n1s-n2s)/(n1s+n2s)
        T=2.0*n1s/(n1s+n2s)
        I = np.array([ [[1.0, R[i]],[R[i], 1.0]]/T[i] 
                       for i in xrange(R.shape[0])])
        return I
        
    def L_mat(self, matind): 
        """
        Calculate the propagation matrix L, through a material
        matind: index of the material
        material name : mat = self.layers[matind]
        thickness     : d   = self.t[matind] 
        complex dielectric constants:  self.nk[mat]
        
        return L, a num[y array with shape len(self.WL)x2x2 array
        L[i] is the propogation matrix at wavelength self.WL[i]
        
        L[i] = [ [ exp(-x*d),        0 ]
                 [         0, exp(-x*d)] ]
        where x = n*cos(phi)* 2*(pi)/(lambda),  
        (n:complex, phi:incident angle, here phi= 0
        """
        
        mat, d = self.layers[matind], self.t[matind]  # d: thickness
        x = 2.0*(np.pi)*self.nk[mat]/self.WL
        L = np.array([ [[np.exp((-1j)*x[i]*d), 0], [0 , np.exp(1j*x[i]*d)]]
                        for i in xrange(x.shape[0])])
        return L
        

if __name__=="__main__": 
    
    plt.clf()
    OM = OpticalModeling(Device, WLstep = WLstep, posstep = posstep)
    #OM.RunSim()
    OM.RunSim(plotE, plotAbs, plotGen,
              saveFigE=saveFigE, saveFigAbs=saveFigAbs, saveFigGen=saveFigGen)
    Jsc= OM.JscReport()
    OM.SaveData(SaveName, saveDataE, saveDataAbs, saveDataGen)

    plt.show()