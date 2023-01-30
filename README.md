# Hydrogel-Diffusion-Coefficient-Estimation

Basic code to read data from diffusion experiments using hydrogel beads, then fit a diffusion- and partition coefficient to the data. The partial differential equation describing the diffusion of the various compounds out of the hydrogel bead is discretized by finite differencing. The solution concentration is allowed to vary with time and influences the boundary condition at the surface of the bead. The resulting system of equations is linear and solved using a simple matrix exponential method. The predicted solution concentration is compared to the measured values to enable regression fo the diffusion coefficient and the partition coefficient.

The bootstrap method is applied to estimate a probability distribution describing the parameters. An identifiability analysis is also conducted, which shows both parameters are identifiable.
