# K-to-clust

Using explicit generation of integer partitions to find equilibrium distributions of clusters for small N without relying on the Law of Mass Action. See J. T. Kindt, "Accounting for Finite-Number Effects on Cluster Size Distributions in Simulations of Equilibrium Aggregation", J. Chem. Theor. Comp. v. 9 p. 147-152 (2013).

Fortran 90 code:  K-to-clust.f90

Input file format: (see 60in for sample) First line should contain Nread (integer), N (integer), concentration(real)

(Nread is the number of values for association constants that will be provided. N is the total number of monomers in the system of interest. concentration is N/V in the system of interest.)

Nread further lines will give a cluster size s starting with 1 and an association constant, i.e. the equilibrium constant for the formation of the cluster from free monomers (should equal 1 for s=1)

Final output (see 60out) gives list of cluster sizes, concentration (number/volume) of each cluster size in the finite-N system, concentration in a bulk system at that same total concentration (excluding any clusters bigger than N),association constant (same as input), apparent association constant ( cluster concentration over s-th power of monomer concentration) and number of that cluster in the most probable partition of monomers into clusters.

Input file includes equilibrium association constants for a hypothetical micelle-forming system, same as used for figure 1 of the JCTC (2013) paper.  The output concentrations show a bimodal distribution of cluster sizes, arising from an equilibrium between structures containing 1 large cluster or 2 smaller clusters.
