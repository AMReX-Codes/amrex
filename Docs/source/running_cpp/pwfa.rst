In-depth explanation of a PWFA example
--------------------------------------

This example illustrates the simulation of a PWFA with realistic parameters in the bubble regime. The simulation is specified by 4 particle species, namely, the drive beam (driver), the witness beam (beam), the plasma electron (plasma_e), and the plasma ion (plasma_p). The species parameters are summarized in the following table.

======== ====================================================
Species  Parameters
======== ====================================================
driver   γ = 48923; N = 30x109; σz = 30 μm; σx = σy = 3.45 μm
beam     γ = 48923; N = 10x109; σz = 10 mm; σx = σy = 0.1 mm
plasma_e n = 1x1023 m-3
plasma_p n = 1x1023 m-3
======== ====================================================

The separation between the driver and witness beams is set to 115 μm.

The simulation can be done in the lab frame as well as in a Lorentz-boosted frame, where the computational costs can be substantially reduced. In the lab frame simulation, there is no need to include the plasma ions since they are stationary during the time scale of concern. In a boosted frame, this is no longer valid as they have finite velocities. Therefore plasma_p is also defined in the example without loss of generality.

The simulation parameters are defined in the lab frame, which includes the longitudinal and transverse dimensions of the simulation box, and the diagnostic time snapshots for back-transformed data to the lab frame from a boosted-frame simulation. Thus, when one has defined the grid size in the lab frame, the longitudinal resolution remains the same but the transverse grid sizes need to be adjusted approximately in the boosted frame with the following relation

The time step in the boosted frame is increased as

Here γ is the Lorentz factor of the boosted frame. In the boosted frame with β close to 1 in the forward direction of the beam propagation, the beam length and plasma length change, respectively, according to

Define the total run time of a simulation by the full transit time of the beam through the plasma, and they are given by, respectively in the lab and boosted frame



assuming the plasma moving at c opposite to the beam direction. Thus the number of time steps in the lab and boosted frame are

It should be pointed out that this example is performed in 2D x-y geometry, which is not equivalent to the realistic simulation. However, the fast turnaround time in 2D simulation helps determine the numerical requirements and the optimized boosted frame, which can then be used in 3D simulations.

