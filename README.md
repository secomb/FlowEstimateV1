FlowEstimateV1

The purpose of this program is to estimate blood flows in a network of microvessels with known geometry, but with incomplete boundary conditions. If a complete set of boundary conditions is given, then NetFlowV1 can be used: see https://physiology.arizona.edu/people/secomb/netflow.

The essential idea of the method is to solve for a set of flows that satisfies conservation of flow at all internal nodes and also satisfies any known boundary conditions, while minimizing the sum of squared deviations of pressures and flows from target values based on typical hemodynamic properties of such networks. The method is described in:
Fry BC, Lee J, Smith NP, Secomb TW. 2012 Estimation of blood flow rates in large microvascular networks. Microcirculation 19, 530-538. (doi:10.1111/j.1549-8719.2012.00184.x).

An objective function is set up in terms of the unknown nodal pressures. The constraints of flow conservation at internal and known boundary nodes are introduced using Lagrange multipliers. Differentiating with respect to the pressures and the Lagrange multipliers gives a system of linear equations:

                 N                 IUB'
     __                                    __ __  __          __               __
     |                       |              | | p1 |          |                 |
     |                                      | | p2 |          |                 |
     |   H’ = 1 + ktH/kpw    |     -KT      | | p3 |          |p+kt*LMtau/kpw   |
     |         or 1                         | | .  |          |   or fixed p    |  N
     |                       |              | | .  |          |                 |
     |                                      | | pn |          |                 |
     |-----------------------+--------------| |----|  ____    |-----------------|
     |                                      | | ?1 |  ____    |      q01        |
     |                       |              | | ?2 |          |      q02        |
     |                                      | | ?3 |          |      q03        |
     |          -K           |        0     | | .  |          |  also includes  | IUB'
     |                                      | | .  |          |   terms from    | 
     |                       |              | | .  |          | fixed p values  |
     |                                      | | .  |          |        .        |
     |_                      |            __| |_ __|          |_              __|   

where N = nnod = total number of nodes;
Nkp = numberknownpress = number of boundary nodes with known pressures;
Nu  = numberunknown = number of boundary nodes with unknown pressures and flows;
IUB' = number of internal nodes + number of boundary nodes with known flow
      = N - Nu - Nkp;
IUB' + N = matrixdim = 2*N - Nkp - Nu is the dimension of the matrix.

The target shear stress refers only to the magnitude and the flow directions are unknown a priori. This introduces a combinatorial aspect into the optimization problem. A heuristic approach is used. The parameter ktau, which gives the weight of the squared shear stress deviations in the objective function, is increased stepwise, with the target flow direction set in each segment set equal to the computed flow direction at the previous step. The outline of the procedure is:
(1) Start with ktau = 0.
(2) Determine flows.
(3) If a segment flow is negative, set the corresponding flow direction to -1, so the target shear stress is set to be negative.
(4) Repeat until flow directions no longer change for the given ktau value.
(5) Successively double ktau, and repeat (2), (3) and (4).
Flows are given in nl/min, pressures in mmHg, viscosities in cP, lengths and diameters in microns and times in s. In these units, the procedure is repeated until the ratio of ktau to kpress (the pressure deviation weight) is 4.096, as found by Fry et al. (2012). 

The original implementations of the method used a direct method (LU decomposition) to solve the linear system. For large networks (thousands of nodes) this method is slow and may run into memory limitations. The current version provides a choice of a direct solver or an iterative method that takes advantage of the sparseness of the system. The system is arranged so as to make it strictly symmetric and a preconditioner is used to put ones on the diagonals of the H’ and K matrices. Then the conjugate gradient method is used to solve the system.

Notes

1. We supply several sample problems in the sub-directories. For input file formatting, please see the files in these directories. These files contain further explanatory comments.

2. Input files must be placed in the same working directory where the C code resides. They can be copied from the folders containing sample problems. The following input files are used:
Network.dat This file specifies the network structure, and the known boundary conditions. It uses the same structure as Network.dat files for NetFlow and Greens, but with some significant differences. (1) The data on lines 1 to 6 of the file are ignored. (2) The flow rates and hematocrits listed for each segment are ignored. However, the sign of the flow rates in this section are used to compute the number of segments with reversed flows (relative to the given data). (3) In the final section, boundary conditions are included only for boundary nodes where the boundary pressure (type 0) or flow (type 2) is given. 
ContourParams.dat Although no contour plots are generated, the coordinates of three corners are used to define the plane on which the network is projected in picturenetwork.cpp.
RheolParams.dat This gives parameters needed for estimating the apparent viscosity in each segment. Set varyviscosity = 1 to get diameter-dependent viscosity. Set phaseseparation = 1 to compute phase separation in diverging bifurcations. Warning: including this effect leads to unstable flow directions and lack of convergence of the method. It is preferable to assume fixed viscosity.
FlowEstParams.dat These values are specific to this program: target pressure, target shear stress, number of ktau steps, maximum number of iterations for each tau value.

3. The program generates several output files:
Run_summary.out tracks the variation in RMS deviations of pressure and shear stress from target values
Results.out gives statistics of the resulting distributions of variables, and a list of estimated variables for each segment.
NetworkNew.dat recreates the Network.dat file but with the computed flow values (and possibly hematocrits). The list of boundary nodes is now complete, and is given in terms of pressure boundary conditions. 
histogram-logflows.out.ps, histogram-logstress.out.ps and histogram-pressures.out.ps give histograms of the data in the corresponding .out files, for log10(flow), log10(shear stress) and pressure.
network.exelem, network.exenode and cmgui.com.txt can be used to obtain 3D visualizations of the network using CMGUI. When CMGUI is started, use File - Open - com file, select greens.com.txt and hit “All” for the visualization. For details of cmgui, see: http://sourceforge.net/projects/cmiss/files/cmgui/cmgui-wx-2.8.0/

4. We have tested this package using Microsoft Visual C++ 2017 under Windows 7. For error reporting and suggestions please contact Dr. Timothy W. Secomb, (520) 626-4513, email secomb@u.arizona.edu. We welcome your comments and suggestions.

5. Brendan Fry, Jack Lee, Nic Smith, Amy Smith, Bohan Li and Tim Secomb contributed to the development of this method.

6. This program is freely available for non-commercial use, provided appropriate acknowledgement is given. Commercial users please contact us before using this program. No assurance is given that it is free of errors and any use is at the user’s risk.

Updated August 2018
