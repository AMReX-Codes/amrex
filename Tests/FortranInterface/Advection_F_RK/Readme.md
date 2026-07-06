# Advection_F_RK

This test advects a single scalar field with a face-centered velocity field
using the Fortran interfaces of AMReX.

It is based on `Tests/FortranInterface/Advection_F`, but removes particles and
adds explicit Runge-Kutta time integration options:

- `myamr.time_integrator=rk2`
- `myamr.time_integrator=rk3`
- `myamr.time_integrator=rk4`

The RK3 and RK4 paths exercise the Fortran interface to
`amrex::FillPatcher<MultiFab>` for coarse/fine boundary filling during RK
stages. The RK3 case is registered in CTest.

The directory `Exec/SingleVortex` includes a GNUmakefile and sample inputs
files. Plotfiles are generated when enabled by the inputs and can be viewed
with AMReX-compatible visualization tools.
