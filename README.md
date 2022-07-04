# GBTStraightPipeLinear
FEM implementation of GBT for analysis of thin-walled straight pipe members 


GBTStraightPipeLinear is a C++ code currently under development for stress and deformation analysis of straight thin-walled circular pipes based on the Generalized Beam Theory (GBT). This code is based on the PhD thesis: Generalized Beam Theory (GBT) for the analysis of thin-walled circular pipe members, https://e-pub.uni-weimar.de/opus4/frontdoor/index/index/docId/4572 , chapter 2

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

