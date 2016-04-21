# Euler1D

C++ implementations of the Godunov, WAF, MUSCL, FORCE, FLIC, and SLIC methods for solving the 1D Euler equations, as described in [1]. Each of these methods can be found in the files bearing their names. Possible inputs are given below:

* vector<state> X: The initial condition of the cells of the system. The state class is defined in variables.h and holds the cell's variables.
* int BCL: The left boundary condition (0 for transitive, 1 for reflexive)
* int BCR: The right boundary condition, similarly to the left
* double dx: The spatial distance between neighbouring cells
* double C: The CFL number
* double t: The final time
* double w: Parameter omega used in the MUSCL and SLIC solvers
* int l: The type of limiter used (if appropriate). 0 - SUPERBEE, 1 - VANLEER, 2 - VANALBADA, 3 - MINBEE
 
Note that implementations of both the approximate HLLC and exact Riemann solvers are included.

[1] http://www.springer.com/gb/book/9783540252023
