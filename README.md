# autoDarc
This code parses the output of Nigel Badnell's `AUTOSTRUCTURE` atomic code (https://amdpp.phys.strath.ac.uk/autos/) and outputs a `TARGET.INP` that can be used in the `DARC` R-matrix codes of Norrington et al (https://connorb.freeshell.org).

The `AUTOSTRUCTURE` atomic code is a 'kappa-averaged' relativistic code - outputting non-relativistic orbitals (see Cowan's book). In order to connect these with a relativistic code, the package held in this repository conversts them to their relativistic counterparts via the relations,

$$\begin{array}{c}
P_{n\kappa}(r) = P_{nl}(r),\\
Q_{n\kappa}(r) = \frac{\alpha}{2} \Big(\frac{d}{dr} + \frac{\kappa}{r}\Big) P_{nl} + O(\alpha^2),
\end{array}$$

which are then re-normalised. To do this conversion, the code makes a selection of points from `AUTOSTRUCTURE`'s grid, and then interpolates them on to the default grid in the `DARC` and `GRASP` codes,

$$
r_j = R \big( e^{(j-1)h} -1 \big),
$$

using routines from the standard numerical `FORTRAN` recipes. The orbitals are then able to be used in `DARC` in a CI sense - using these semi-relativistic orbitals for an R-matrix calculation. Alternatively, these can be used as a starting point in `GRASP` calculations. The code outputs an `ORBIN.DAT` file which can be used for `GRASP0` calculations (https://connorb.freeshell.org), and additionally a `rwfn.out` file that could be used for `GRASP2018` calculations (https://github.com/compas/grasp?tab=readme-ov-file). This is advantageous as the orbitals of `AUTOSTRUCTURE` can be scaled to suit the needs of the user.

The code is compiled simply with your compiler of choice, e.g 
`gfortran autodarc.f90 -o autodarc.x`. The code is a simple utility one, and as such does not warrant anything more complex.
