This tutorial initializes particles, and steps through "nsteps" time steps where in each time step 
we 
 * compute the timestep = cfl * "cutoff" (particle radius) / max_particle_vel
 * compute or update the grid neighbors 
 * calculate particle neighbor lists
 *  compute forces due to particle-particle collisions
 * update the particle velocities then particle positions.   

At every time step we print out dt and the number of particles.
