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

Brendan Fry, Amy Smith, Bohan Li and Tim Secomb contributed to the development of this method.

Updated August 2018
