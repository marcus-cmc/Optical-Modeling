# Optical-Modeling :
Modeling the light propogation in a multi-layer thin-film stack using a Transfer Matrix method.

##
## OpticalModeling class
The `OpticalModeling` object is designed to model the light propogation in a stack of thin films with different materials at normal incidence. It can be used to calculate the following properties for thin films:
* **Light absorption** 
* **Transmission**
* **Reflection** 

And the following properties in solar cells under standard AM 1.5 solar irradiation at normal incidence:
* **Electric field profile**
* **Charge carrier generation rates** (equivalent to **photon absorption rate**)
* __*Jsc*__ (short circuit current density, assuming 100 % IQE, i.e. all absorbed photons are converted into charge carriers)

### Examples
Below are some example output figures for a device stack consisting of these materials at the given thickness (their refraction indices are included in the `Index_of_Refraction_library_Demo.csv` file):

Layer No. | Material| Thickness| |note|
|---|---|---|---|---|
|0| Glass| N/A | |substrate|
|1| **ITO** | 145 nm | |transparent electrode|
|2| **ZnO** | 120 nm | |window layer|
|3| **PbS** | 250 nm | |main light absorbing layer|
|4| **Au**  | 150 nm | |gold electrode|


#
#### Absorption, Transmission, Reflection of each material in the device stack
<img src="/Example_OpticalModeling_Figures/Fig_Absorption.png" width="600" "Absorption">


#
####Carrier generation (photon absorption) rate at different wavelengths and position.


<img src="/Example_OpticalModeling_Figures/Fig_AbsorptionRate.png" width="600" "Carrier generation (photon absorption) rate">


#
####Position vs carrier generation (photon absorption) rate in the device. 
This is basically the summation of the generation rates over different wavelengths for each position in the figure shown above)

<img src="/Example_OpticalModeling_Figures/Fig_Gen_position_.png" width="600" "Generation rave vs position">


#
####Electric field profile in each material at different wavelengths and position.
<img src="/Example_OpticalModeling_Figures/Fig_Efield.png" width="600" "E-field map">


#
####Position vs electric field in the device at selected wavelengths 
(selected slices of the figure above)

<img src="/Example_OpticalModeling_Figures/Fig_Efield_selectedWL.png" width="600" "E-field selected wl">



Finally, the code would print out a brief summary that looks like this in the console :

```python
"""
Summary of the modeled results between 350.0 and 1200.0 nm

Layer No. Material  Thickness (nm)  Jsc_Max (mA/cm^2)
        1      ITO             145               2.61
        2      ZnO             120               0.97
        3      PbS             250              24.07
        4       Au             150               1.47
"""
```
#
===

#

##
## OMVaryThickness

The `OMVaryThickness` object is a subclass of `OpticalModeling`. It adds some features so that you can vary the thickness of **one** or **two** layers in the thin film stack to see how the properties of interest – either **absorption** (and thus __*Jsc*__) in '*any*' layer, **transmission**, or **reflection**) – respond to the change of thickness.


For example, you can change the thickness of `layer2`, and watch how the following properties changes with respect to it:
* **Transmission** of the whole stack
* **Reflection** of the whole stack
* **Absorption** in `layer1`, `layer2`, `layer3`...

This feature would be very helpful for design of experiments. Below are some example outputs by using the `OMVaryThickness` object.



#
## Examples for OMVaryThickness
These examples use the same device stack to that in the example for `OpticalModeling`.
First create a `VaryThickness` object and then call the method `VaryOne()` with the desired parameters to run the simulation


#### OMVaryThickness Example 1
In this example, we vary the thickness of the **PbS** layer (`ToVary=3`) and set the target to **PbS** itself (`target=3`). 
```python
VT1 = OMVaryThickness(Device, libname="Index_of_Refraction_library_Demo.csv", WLrange=[350, 1200])
VT1.VaryOne(ToVary = 3, t_range = range(50, 601, 10), target = 3)
```

The top figure shows how the absorption in the PbS changes with its own thickness and the bottom one shows how the Jsc changes with it.


<img src="/Example_VaryThickness_Figures/VaryPbS_AbsPbS.png" width="800" >
<img src="/Example_VaryThickness_Figures/VaryPbS_JscPbS.png" width="480" >

# 
#### OMVaryThickness Example 2
In the second example, we vary the thickness of the **ZnO** layer (`ToVary = 2`)to see how the properties of the **PbS** layer (`target = 3`) change with it. 
```python
VT2 = OMVaryThickness(Device, libname="Index_of_Refraction_library_Demo.csv", WLrange=[350, 1200])
VT2.VaryOne(ToVary = 2, t_range = range(20, 551, 10), target = 3)
```
As you can see, because of the interference effects in thin films, the absorption of the **PbS** layer is not a monotonic function of the thickness of **ZnO**. It also shows strong wavelength dependence: every wavelength shows a different behavior. This example demostrates the usefulness of the `OMVaryThickness` object for the design of experiments.

<img src="/Example_VaryThickness_Figures/VaryZnO_AbsPbS.png" width="800" >
<img src="/Example_VaryThickness_Figures/VaryZnO_JscPbS.png" width="480" >

Once the simulation has been done after calling `VaryOne()`, we can just call the `PlotVaryAbs(target)` and `PlotVaryJsc(target)` methods to generate this plot again or generate other plots (different `target` layers) without running the simulation again. We can call 
```python
VT2.PlotVaryAbs(target = 2) # ZnO layer
```
to get the absorption in the ZnO layer, which increases with the thickness of itself:  


<img src="/Example_VaryThickness_Figures/VaryZnO_AbsZnO.png" width="800">


####
#### 
 
 
If, instead, we use target 1 (ITO layer)
```python
VT2.PlotVaryAbs(target = 1) # ITO layer
```
we get the absorption of the ITO layer with respect to the thickness of the ZnO layer. (ITO is very transparent to the visible light but it could show strong absorption in the near-infrared (>750nm) ). This figure probably does not provide too much useful information in practical, but I found the inteference pattern very beautiful (science!) -- I intentionally simulated a lot of data points to make the overlap of all curves look more interesting. That's why I decided to show it here.


<img src="/Example_VaryThickness_Figures/VaryZnO_AbsITO.png" width="800">
#### 
#### 
#### 
In addition to the absorption in each layer, we can also use `"T"` to plot the transmision,`"A"` for the total absorption (sum over all layers), and `"R"` for reflection in the device stack with respect to the thickness of the `ToVary` layer.
```python
VT2.PlotVaryAbs(target = "T") # Transmission
VT2.PlotVaryAbs(target = "A") # absorption
VT2.PlotVaryAbs(target = "R") # Reflection
```
The transmission is very low (< 1%) because there is an opaque electrode (150 nm of gold).
<img src="/Example_VaryThickness_Figures/VaryZnO_T.png" width="800" >

The overall absorption in the device stack
<img src="/Example_VaryThickness_Figures/VaryZnO_A.png" width="800" >

The reflection of the device: everything not absorbed is reflected back!
<img src="/Example_VaryThickness_Figures/VaryZnO_R.png" width="800">


### How to Run OpticalModeling
There is a **`RunModeling.py`** file desgined to take user inputs and then run the optical modeling based on your inputs. To run the simulation, simply provide a library of the refraction indices of the materials of interest and specify the materials and thickness in the thin film stack in the **`RunModeling.py`** file and the run it. You can get the some output figures with options to save the data as .csv files and figure in your desired format (such as .pdf vector graphics or .png raster graphics). More information on how to run it and how to specify inputs, and how to save data/figures can be found in the comments in the **`RunModeling.py`** file.

### How to Run VaryThickness
There is a **`RunVaryOneThickness.py`** file desgined to take user inputs and then run the optical modeling based on your inputs.







