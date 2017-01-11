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

### How to Run OpticalModeling
There is a **`RunModeling.py`** file desgined to take user inputs and then run the optical modeling based on your inputs. To run the simulation, simply provide a library of the refraction indices of the materials of interest and specify the materials and thickness in the thin film stack in the **`RunModeling.py`** file and the run it. You can get the some output figures with options to save the data as .csv files and figure in your desired format (such as .pdf vector graphics or .png raster graphics). More information on how to run it and how to specify inputs, and how to save data/figures can be found in the comments in the **`RunModeling.py`** file.


### Examples
Below are some example output figures for a device stack consisting of these materials at the given thickness (their refraction indices are included in the `Index_of_Refraction_library_Demo.csv` file):
* Glass (substrate)
* ITO (145 nm)
* ZnO (120 nm)
* PbS (250 nm)
* Au (150nm)

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
## VaryThickness

The `VaryThickness` object is a subclass of `OpticalModeling`. It adds some features so that you can vary the thickness of **one** or **two** layers in the thin film stack to see how your desired results – either **absorption** (and thus __*Jsc*__) in __'any'__ layer, **transmission, or **reflection**) – respond to the change of thickness.


For example, you can change the thickness of `layer2`, and watch how the following properties changes with respect to it:
* **Transmission** of the whole stack
* **Reflection** of the whole stack
* **Absorption** in `layer1`, `layer2`, `layer3`...

This feature would be very helpful for design of experiments. Below are some example outputs by using the `VaryThickness` object.


#
#### VaryThickness Example 1
In this example, we vary the thickness of the `PbS` layer. The top figure shows how the absorption in the `PbS` changes with its own thickness and the bottom one shows how the Jsc changes with it.

<img src="/Example_VaryThickness_Figures/VaryOneLayer_Abs2.png" width="800" "VaryPbs - PbsAbs">
<img src="/Example_VaryThickness_Figures/VaryOneLayer_Jsc2.png" width="480" "VaryPbs - PbsJsc">

# 
#### VaryThickness Example 2
In the second example, we vary the thickness of the `ZnO` layer and see how the properties of the `PbS` layer change with it.

<img src="/Example_VaryThickness_Figures/VaryZnO_targetPbSAbs.png" width="800" "VaryZnO - PbSAbs">
<img src="/Example_VaryThickness_Figures/VaryZnO_targetPbS_Jsc.png" width="480" "VaryZnO - PbSJsc">









