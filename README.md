# GBTStraightPipeLinear
FEM implementation of GBT for analysis of thin-walled straight pipe members 


GBTStraightPipeLinear is a C++ code currently under development for stress and deformation analysis of straight thin-walled circular pipes based on the Generalized Beam Theory (GBT). This code is based on the PhD thesis: Generalized Beam Theory (GBT) for the analysis of thin-walled circular pipe members, https://e-pub.uni-weimar.de/opus4/frontdoor/index/index/docId/4572 , chapter 2



## Getting started

We recomend to use Visual Studio IDE and vcpkg C++ library manager.  Before you run the program, please update the path of the vcpkg.cmake file in the CMakeLists.txt.


## Dependencies
 * C++ 20
 * Libraries: Eigen, OpenMP, matplot++ and gnuplot

## Code structure
 src\
 * main.cpp -> preprocessing, solver and postprocessing
 
 include\
 
 preprocessing headers 
 
 * boundary_conditions.h	
 * LM_connectivity_matrix.h
 * ik_tensor.h
 * mode_coupling_ik.h	
 * element_stiffness_matrix.h
 * loading.h	
 
 postprocessing headers
 
 * GBT_function.h	
 * amplitude_function.h	
 * displacements.h	
 * mesh.h
## Numerical Example
The numerical example is developed to validate and illustrate the application and capabilities of the linear GBT formulation and its numerical implementation. Here, a short cantilever pipe is considered as a numerical example with the physical properties and boundary conditions shown below.

![example](https://github.com/AbinetKH/GBTStraightPipeLinear/blob/master/doc/example.png)



### The generalized modal amplitude vector
![The generalized modal amplitude vector](https://github.com/AbinetKH/GBTStraightPipeLinear/blob/master/doc/dispVector.png)

### Force vector
![Force vector](https://github.com/AbinetKH/GBTStraightPipeLinear/blob/master/doc/externalForceVector.png)

### Element stiffness matrix
![Element stiffness matrix](https://github.com/AbinetKH/GBTStraightPipeLinear/blob/master/doc/stiffnessmatrix.png)

### Deformation shape of a short cantilever pipe
![Deformation shape of a short cantilever pipe](https://github.com/AbinetKH/GBTStraightPipeLinear/blob/master/doc/plot.png)

Citation
--------

If you use any of the GBTStraightPipeLinear packages to produce scientific articles, please cite: https://e-pub.uni-weimar.de/opus4/frontdoor/index/index/docId/4572
