
#include <BLProfiler.H>

#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::Init(MultiFab& dummy_mf)
{
    BL_PROFILE("MyPC::Init()");

    Real x_off = 0.5;
    Real y_off = 0.5;
    Real z_off = 0.5;

    Real charge = 1.;

    // Initialize one particle per cell with local position (x_off,y_off,z_off) relative to nondimensional cell size [0:1]
    InitOnePerCell(x_off,y_off,z_off,charge,dummy_mf);

    // Randomly initialize "num_particles" number of particles, each with charge "charge"
    // bool serialize = false;
    // int iseed   = 10;
    // InitRandom(num_particles,iseed,charge,serialize);
}
