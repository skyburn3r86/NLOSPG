# NLOSPG
Script for controlling COMSOL 5.4 Simulation of Integrated Single Photon Source. Comsol provides the orthogonal basis for the 
Quantum State Simulations.
## Main
Loads the main parameters for the Simulation. Everything can be controlled from here. For sucessful simulation, start with the 
Function Setup, which receives the dimensions of the device in varargin. The default settings are already given. Next, load 
the parameters for the Dispersion relation into the model. additionally, you should run the geometry function to build the 
structure from the given parameters. After the structure is built, Materials need to be assined, finally, the mesh should be built, 
the Physics transferred to the model and computed. 
## Quantum Evaluation
### Classes:
**Class for Classical Field**:  
Fields normalized such that they fit into the other Quantum Mechanical framework. Instead of Photon numbers, 
takes the limit for very large N and takes input power and relates the mean Photon number to the Photon flux [Power] (<img src="https://latex.codecogs.com/gif.latex?P=dN/dt\hbar\omega "/> ) by assuming 
a linearized relationship <img src="https://latex.codecogs.com/gif.latex?F=\Delta N/\Delta t "/>  
**Class for Quantum Field**:  
Initializes the Quantum field including number states. No coherent states initialized
### Functions:
**overlap**:  
Calculates the g0 coupling rate  
**Perturbation Losses**:  
Calculates the losses in a perturbative approach. Very ineffective formulism, O(NPhotons^Order) speed, with Order being the perturbative order.  
**Propagate**:  
Propagates the Quantum and Classical fields by along the vector z  
**QuantumTest**:  
Test for the Propagation and evaluation  
**test_ExtractField**:  
Extracting the and Selecting Modes from the Matlab Simulation
## Plotting:
Some functions for Plotting and extracting results.

## Simulation:
Skeleton for Generation of the COMSOL model.

## Evalution:
Function for Evaluation of the Data.