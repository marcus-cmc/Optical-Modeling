# Optical-Modeling (Transfer Matrix Method):
Modeling light propogation in a multi-layer thin-film stack

## OpticalModeling class
The `OpticalModeling` object is designed to model the light propogation in a stack of thin films with different materials at normal incidence. It can be used to calculate the following properties for thin films:
* **Light absorption** 
* **Transmission**
* **Reflection** 

And the following properties in solar cells under standard AM 1.5 solar irradiation at normal incidence:
* **Electric field profile**
* **Charge carrier generation rates** (equivalent to **photon absorption rate**)
* __*Jsc*__ (short circuit current density, assuming 100 % IQE, i.e. all absorbed photons are converted into charge carriers)


Below are some examples output figures for a device stack consisting of these materials at the given thickness:
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
#
### VaryThickness

The `VaryThickness` object is a subclass of `OpticalModeling`. It adds some features so that you can vary the thickness of **one** or **two** layers in the thin film stack to see how your desired results – either **absorption** (and thus __*Jsc*__) in __'any'__ layer, **transmission, or **reflection**) – respond to the change of thickness.


For example, you can change the thickness of `layer2`, and watch how the following properties changes with respect to it:
* **Transmission** of the whole stack
* **Reflection** of the whole stack
* **Absorption** in `layer1`, `layer2`, `layer3`...
This feature would be very helpful for design of experiments. Below are some example outputs by using the `VaryThickness` object.







