# SKETCH-N

## Description

### Program Name and Title 
SKETCH-N: A nodal code for solving neutron diffusion equations of steady-state and kinetics problems.

### List of authors
Zimin V.G., Semenov A.A., Romanenko V.I., Glazkov O.V.

### Simulation scope 
The SKETCH-N code solves neutron diffusion equations in x-y-z and hexagonal-z geometries for steady-state and neutron kinetics problems. The code can treat an arbitrary number of neutron energy groups and delayed neutron precursors.

### Methods used 
The polynomial, semi-analytic and analytic nodal methods based on the nonlinear iteration procedure can be used for spatial discretization of diffusion equations. The time integration of the neutron kinetics problem is performed by a fully implicit scheme with an analytical treatment of the delayed neutron precursors. The steady-state eigenvalue problems are solved by inverse iterations with Wielandt shift, the Chebyshev adaptive acceleration procedure is used for the neutron kinetics problems. The block symmetric Gauss-Seidel preconditioner is applied in the both iterative methods. The flux-weighting homogenization procedure is used for partially-rodded nodes to minimize a rod cusping effect. Simple one-phase model of the thermal-hydraulics of fuel assembly is included in the code. The code also has an interface module for a coupling with transient analysis codes , such as TRAC. The interface module performs a data exchange between the codes, synchronizes a time stepping and maps the neutronics data onto thermal-hydraulics spatial mesh and vice versa. The interface module is based on the message passing library PVM (Parallel Virtual Machine).

### Limitations 
The code can treat the neutron diffusion problems in Cartesian geometry. Few-group macro cross sections and their dependencies are provided by a code user. The code does not have fuel burn-up modelling capabilities. An external thermal-hydraulics code is generally required for the calculation of the “real-life” problems.

### Typical Running Time 
The running time of the full-core case C1 of the PWR NEACRP rod ejection benchmark (2 neutron energy groups, 6 groups of the delayed neutron precursors, 884x18 neutronics nodes, 910 time steps) is 68 minutes on Sun UltraSPARC I (143 MHz) with an internal thermal hydraulics model.

### Code features 
Dimensions of a problem are specified as parameters in the include files, the code should be recompiled when the problem dimensions are changed. The code has PVM- based interface module developed for a coupling with transient thermal-hydraulics codes. The interface model has been used for a coupling of the SKETCH-N code with the J-TRAC (TRAC-PF1) and TRAC-BF1 codes.

### Related and Auxiliary Programs 
PVM library is used for the interface module of the code.

### Status
The SKETCH-N code has been verified by solving the steady-state and neutron kinetics benchmark problems. The coupled J-TRAC/SKETCH-N code system has been verified against NEACRP PWR rod ejection and rod withdrawal benchmarks. NEACRP BWR cold water injection benchmark has been used for verification of the TRAC-BF1/SKETCH-N system.

### References
* Zimin V. G. “Nodal Neutron Kinetics Models Based on Nonlinear Iteration Procedure for LWR Analysis”, PhD Thesis, Research Laboratory for Nuclear Reactors, Tokyo Institute of Technology, August 1997.
* Zimin, V.G., and H. Ninokata, “Nodal Neutron Kinetics Model Based on Nonlinear Iteration Procedure for LWR Analysis,” Ann. Nucl. Energy, 25, 507-528, 1998
* Zimin, V.G., H. Ninokata, and L. R. Pogosbekyan “Polynomial and Semi-Analytic Nodal Methods for Nonlinear Iteration Procedure,” Proc. of the Int. Conf. on the Physics of Nuclear Science and Technology, October 5-8, 1998, Long Island, New York, American Nuclear Society, vol. 2, pp. 994- 1002, 1998
* Zimin V. G., H. Asaka, Y. Anoda, and M. Enomoto, “Verification of the J-TRAC code with 3D Neutron Kinetics Model SKETCH-N for PWR Rod Ejection Analysis”, Proc. of the 9 International Topical Meeting on Nuclear Reactor Thermal Hydraulics (NURETH 9), San Francisco, California, October 3- 8, CD-ROM, 1999.
* Asaka, H., V. G. Zimin, T. Iguchi, and Y. Anoda, “Coupling of the Thermal-Hydraulic TRAC Codes with 3D Neutron Kinetics Code SKETCH-N”, Proc. of OECD/CSNI Workshop on Advanced Thermal- Hydraulic and Neutronics Codes: Current and Future Applications, Barcelona, Spain, 10-13 April, 2000
* Zimin, V. G., H. Asaka, Y. Anoda, E. Kaloinen and R. Kyrki-Rajamaki, “Analysis of NEACRP 3D BWR Core Transient Benchmark”, Proc. of the 4 Intl. Conf. on Supercomputing in Nuclear Application SNA 2000, September 4-7, 2000, Tokyo, Japan.

### Machine Requirements 
A workstation under UNIX.

### Program Language Used 
Fortran 77 and Fortran 90 

### Other Programming or Operating Information or Restrictions 
The interface module requires [PVM] (http://www.epm.ornl.gov/pvm/pvm\_home.html) installed on a computer. The PVM is a public domain software available from NETLIB.

### Material Available
Source Code
Sample LWR Benchmark Problems Input and Output Files 
SKETCH-N Manual, vol. I. Model Description
SKETCH-N Manual, vol. II User’s Guide

### Keywords
kinetics, three-dimensional, neutron diffusion, nodal methods, nonlinear iteration procedure, reactor transient analysis.

### Contacts
	Vyacheslav G. Zimin
	Division 836 “Laboratory of Simulator Systems”, National Research Nuclear University “MEPhI”, Kashirskoe shosse, 31, Moscow, 115409, Russia
	Tel.: +7-963-648-87-81
	E-mail: vgzimin@mail.ru

