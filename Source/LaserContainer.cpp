
#include <limits>

#include <WarpXConst.H>
#include <ParticleContainer.H>

LaserContainer::LaserContainer (AmrCore* amr_core, int ispecies)
    : WarpXParticleContainer(amr_core, ispecies)
{
    charge = PhysConst::q_e; // note that q_e is defined to be positive.
    mass = std::numeric_limits<Real>::max();
}

void
LaserContainer::AllocData ()
{
    // have to resize here, not in the constructor because GDB was not
    // ready in constructor.
    m_particles.resize(GDB().finestLevel()+1);
}

void
LaserContainer::InitData ()
{
    BoxLib::Abort("LaserContainer::InitData: Not implemented");
}

void
LaserContainer::Evolve (int lev,
			const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
			const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
			MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BoxLib::Abort("LaserContainer::Evolve: Not implemented");
}
