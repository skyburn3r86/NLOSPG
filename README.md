# NLOSPG
Script for controlling COMSOL 5.4 Simulation of Integrated Single Photon Source. Comsol provides the orthogonal basis for the 
Quantum State Simulations.
## Main
Loads the main parameters for the Simulation. Everything can be controlled from here. For sucessful simulation, start with the 
Function Setup, which receives the dimensions of the device in varargin. The default settings are already given. Next, load 
the parameters for the Dispersion relation into the model. additionally, you should run the geometry function to build the 
structure from the given parameters. After the structure is built, Materials need to be assined, finally, the mesh should be built, 
the Physics transferred to the model and computed. 
