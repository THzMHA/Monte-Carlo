Two examples of simulating the propagation of polarized THz wave in turbid medium
1. Run plot_mie.m to display single-particle propeties based on Mie theory in absorbing medium.
2. Run mie_polar.m to obtain the diffuse reflectance and degree of polarization of exiting THz wave.

mie_polar.m uses a executable file iquv.exe, which is modified and compiled from the [Meridian Plane Monte Carlo](https://omlc.org/software/polarization/index.html) codes written by J. Ramella-Roman in C. The matlab code generates necessary inputs, fetches them to iquv.exe, and analyzes the outputs in the end.

Original version of iquv.exe uses four polarization inputs (H, V, P, R) to obtain the Mueller matrix. In this example, we only simulate P and R since the results of all linear polarizations are virtually the same. Thus, a typical output from iquv.exe looks like:
>
launch P
>
R= 0.02711577	 -0.00005200	 -0.02343410	 0.00000701
>
launch R
>
R= 0.02691129	 0.00000948	 0.00001508	 -0.02034429`
