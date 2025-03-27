# autoDarc
This code parses the output of Nigel Badnell's `AUTOSTRUCTURE` atomic code (https://amdpp.phys.strath.ac.uk/autos/) and outputs a `TARGET.INP` that can be used in the `DARC` R-matrix codes of Norrington et al (https://connorb.freeshell.org).

The `AUTOSTRUCTURE` atomic code is a 'kappa-averaged' relativistic code - outputting non-relativistic orbitals. These radial orbitals are then converted to their relativistic counterparts via the relations,

$$\begin{array}{ccc}
P_{n\kappa}(r) = &P_{nl}(r),\\
Q_{n\kappa}(r) = &\frac{\alpha}{2} \Big(\frac{d}{dr} + \frac{\kappa}{r}\Big) P_{nl} + O(\alpha^2),
\end{array}$$
which are then normalised. To do this, the code makes a selection of points from `AUTOSTRUCTURE`'s grid, and then interpolates them on to the default grid in the `DARC` and `GRASP` codes,
