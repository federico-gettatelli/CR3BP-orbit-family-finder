# CR3BP-orbit-family-finder
Computes the initial states of periodic orbits in the circular restricted three body problem. Works with lyapunov and halo orbits. Still wip.

# HOW TO USE
In main.m you can find all the parameters to change to find the initial conditions of the orbits of a determined family.
The algorithm used to find the initial condition is a correction method that starts from an initial state and corrects it with the informations from the state transition matrix after half a period. In particular the variable "CHANGE" define the x/z distance between successive orbits. If "CHANGE" is to high the algorithm is not able to converge to a solution and the code stops. So "CHANGE" should be a small value, for halo < 1000-5000 km, for lyapunov < 100-500 km. 
