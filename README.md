# three-way_bulging
these codes are used to analyze a system of membrane tubules under different conditions. 
The most basic system is a simple tube or a cylinder made of a lipid membrane. The more complex component of the tubular network is the tubular junction, modeled by three tubules connected via a three-way junction.
Here, we wish to examine different conditions applied to the system, or, in other words, different ensembles constraining or controlling different parameters, such as tension/surface area, and pressure/volume.
in order to numerically simulate the structure we are using Ken Brakke's "Surface Evolver", thus, we minimize the energy and find the specific relations of the different parameters in each specific system. 
in order to process and utilize the energy minimization results, we use MATLAB. By this analysis, we are analyzing the different relations, interpolate the governing rules of the system, that emerge in the master curves, and lastly, use the master curve to simulate any ensemble we wish to generate.

Here, you can find files used to simulate different scenarios, for example, one can use these files to simulate an experiment with fixed volume, and varying tension, or, alternatively, varying pressure, and fixed tension. Different ensembles with different boundary conditions can be simulated. Then, one can use the MATLAB codes in order to analyze the results. One would find out that under the normalization we offer, all the ensembles, with all the different boundary conditions, fit on a single curve we call- "master curve".
*Note- the simulation must not violate a simple condition which is that the ends of the tubules must be long enough to allow full relaxation and convergence into a simple cylindrical shape.
