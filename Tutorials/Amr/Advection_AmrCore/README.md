 use_math: true

##  Multi-Level Scalar Advection

### What Features Are We Using

* Mesh data 
* Dynamic AMR with and without subcycling

### The Problem

Consider a drop of dye (we'll define $$\phi$$ to be the concentration of dye) 
in a thin incompressible fluid that is spinning 
clock-wise then counter-clockwise with a prescribed motion.  We consider the dye to be a 
passive tracer that is advected by the fluid velocity.  The fluid is thin enough that we can model
this as two-dimensional motion; here we have the option of solving in a 2D or 3D computational domain.

In other words, we want to solve for $$\phi(x,y,t)$$ by evolving 

$$\frac{\partial \phi}{\partial t} + \nabla \cdot (\bf{u^{spec}} \phi)  = 0$$

in time ($$t$$), where the velocity $${\bf{u^{spec}}} = (u,v)$$ is a divergence-free field computed by defining

$$\psi(i,j) = \sin^2(\pi x) \sin^2(\pi y)  \cos (\pi t / 2) / \pi $$

and defining

$$u = -\frac{\partial \psi}{\partial y},  v = \frac{\partial \psi}{\partial x}.$$

Note that because $${\bf{u^{spec}}}$$ is defined as the curl of a scalar field, it is analytically divergence-free

In this example we'll be using AMR to resolve the scalar field since the location of the dye is
what we care most about.

### The AMR Algorithm

To update the solution in a patch at a given level, we compute fluxes ($${\bf u^{spec}} \phi$$)
on each face, and difference the fluxes to create the update to phi.   The update routine
in the code looks like

```cplusplus
  // Do a conservative update
  {
    phi_out(i,j,k) = phi_in(i,j,k) +
                ( AMREX_D_TERM( (flxx(i,j,k) - flxx(i+1,j,k)) * dtdx[0],
                              + (flxy(i,j,k) - flxy(i,j+1,k)) * dtdx[1],
                              + (flxz(i,j,k) - flxz(i,j,k+1)) * dtdx[2] ) );
  }
```

In this routine we use the macro AMREX_D_TERM so that we can write dimension-independent code; 
in 3D this returns the flux differences in all three directions, but in 2D it does not include
the z-fluxes.

Knowing how to synchronize the solution at coarse/fine boundaries is essential in an AMR algorithm;
here having the algorithm written in flux form allows us to either make the fluxes consistent between
coarse and fine levels in a no-subcycling algorithm, or "reflux" after the update in a subcycling algorithm.

The subcycling algorithm can be written as follows
```C++
void
AmrCoreAdv::timeStepWithSubcycling (int lev, Real time, int iteration)
{

    // Advance a single level for a single time step, and update flux registers
    Real t_nph = 0.5 * (t_old[lev] + t_new[lev]);
    DefineVelocityAtLevel(lev, t_nph);
    AdvancePhiAtLevel(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStepWithSubcycling(lev+1, time+(i-1)*dt[lev+1], i);
        }

        if (do_reflux)
        {
            // update lev based on coarse-fine flux mismatch
            flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);
        }

        AverageDownTo(lev); // average lev+1 down to lev
    }

}
```

while the no-subcycling algorithm looks like
```C++
void
AmrCoreAdv::timeStepNoSubcycling (Real time, int iteration)
{
    DefineVelocityAllLevels(time);
    AdvancePhiAllLevels (time, dt[0], iteration);

    // Make sure the coarser levels are consistent with the finer levels
    AverageDown ();

    for (int lev = 0; lev <= finest_level; lev++)
        ++istep[lev];
}
```

### Running the Code

First git clone amrex.  Then from within the amrex repo,

```
cd Tutorials/Advection_AmrCore/Exec
```
Note that you can choose to work entirely in 2D or in 3D ... whichever you prefer.
The instructions below will be written for 3D but you can substitute the 2D executable.

To build in 2d, type
```
make DIM = 2
```

To build in 3d, type
```
make DIM = 3
```

In this directory you'll see 


```
inputs -- an inputs file for both 2D and 3D
```

To run in serial, 

```
./main3d.gnu.MPI.ex inputs
```

To run in parallel, for example on 4 ranks:

```
mpiexec -n 4 ./main3d.gnu.MPI.ex inputs
```

The following parameters can be set at run-time -- these are currently set in the inputs
file but you can also set them on the command line.  

```
stop_time          =  2.0                # the final time (if we have not exceeded number of steps)
max_step           = 1000000             # the maximum number of steps (if we have not exceeded stop_time)

amr.n_cell         =  64  64   8         # number of cells at the coarsest AMR level in each coordinate direction

amr.max_grid_size  = 16                  # the maximum number of cells in any direction in a single grid

amr.plot_int       = 10                  # frequency of writing plotfiles

adv.cfl            = 0.9                 # cfl number to be used for computing the time step

adv.phierr = 1.01  1.1  1.5              # regridding criteria  at each level

```

The base grid here is a square of 64 x 64 x 8 cells, made up of 16 subgrids each of size 16x16x8 cells.  
The problem is periodic in all directions.

We have hard-wired the code here to refine based on the magnitude of $$\phi$$.    Here we set the 
threshold level by level.  If $$\phi > 1.01$$ then we want to refine at least once; if $$\phi > 1.1$$ we
want to resolve $$\phi$$ with two levels of refinement, and if $$\phi > 1.5$$ we want even more refinement.

Note that you can see the total runtime by looking at the line at the end of your run that says

```
Total Time: 
```

and you can check conservation of $$\phi$$ by checking the line that prints, e.g. 

```
Coarse STEP 8 ends. TIME = 0.007031485953 DT = 0.0008789650903 Sum(Phi) = 540755.0014
```

Here Sum(Phi) is the sum of $$\phi$$ over all the cells at the coarsest level.

Questions to answer:

```
1. How do the subcycling vs no-subycling calculations compare?
    a.   How many steps did each take at the finest level? Why might this not be the same?
    b.   How many cells were at the finest level in each case? Why might this number not be the same?

2  What was the total run time for each calculation?  Was this what you expected?

3. Was phi conserved (over time) in each case?
      a.  If you set do_refluxing = 0 for the subcycling case, was phi still conserved?
      b.  How in the algorithm is conservation enforced differently between subcycling and not?

4. How did the runtimes vary with 1 vs. 4 MPI processes?  
   We suggest you use a big enough problem here -- try running 

   mpiexec -n 1 ./main3d.gnu.MPI.ex inputs_for_scaling

   mpiexec -n 4 ./main3d.gnu.MPI.ex inputs_for_scaling

5. Why could we check conservation by just adding up the values at the coarsest level?
```

### Visualizing the Results

Here is a sample slice through a 3D run with 64x64x8 cells at the coarsest level and three finer levels (4 total levels).

![Sample solution](amr101_3D.gif)

After you run the code you will have a series of plotfiles.  To visualize these
we will use ParaView 5.8, which has native support for AMReX Grid, Particle,
and Embedded Boundary data (in the AMR 101 exercise we only have grid data).

#### Make a Movie with the ParaView 5.8 Script

To use the ParaView 5.8 python script, simply do the following to generate `amr101_3D.gif`:

```
$ make movie3D
```

If you run the 2D executable, make the 2D movie using:

```
$ make movie2D
```

Notes:

- To delete old plotfiles before a new run, do `rm -rf plt*`

- You will need `+ffmpeg`.

#### Using ParaView 5.8 Manually

To do the same thing with ParaView 5.8 manually,

```
1. Start Paraview 5.8
2. File --> Open ... and select the collection of directories named "plt.." --> [OK]
3. From the "Open Data With..." dialog that pops up, select "AMReX/BoxLib Grid Reader" --> [OK]
4. Check the "phi" box in the "Cell Array Status" menu that appears
5. Click green Apply button
6. Click on the "slice" icon -- three to the right of the calculator
   This will create "Slice 1" in the Pipeline Browser which will be highlighted.
7. Click on "Z Normal"
8. Unclick the "Show Plane" button
9. Click green Apply button
10. Change the drop-down menu option (above the calculator row) from "vtkBlockColors" to "phi"
```

You are now ready to play the movie!  See the "VCR-like" controls at the top. Click the play button.

### Additional Topics to Explore

* What happens as you change the max grid size for decomposition?

* What happens as you change the refinement criteria (i.e. use different values of $$\phi$$)?
  (You can edit these in inputs)  
