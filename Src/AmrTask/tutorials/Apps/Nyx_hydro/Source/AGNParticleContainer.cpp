#include "AGNParticleContainer.H"
#include "AMReX_RealVect.H"
#include "agn_F.H"

using namespace amrex;

using std::cout;
using std::endl;

void
AGNParticleContainer::moveKickDrift (amrex::MultiFab&       acceleration,
		                     int                    lev,
                    		     amrex::Real            dt,
		                     amrex::Real            a_old,
				     amrex::Real            a_half,
				     int                    where_width)
{
    BL_PROFILE("AGNParticleContainer::moveKickDrift()");

    //If there are no particles at this level
    if (lev >= this->GetParticles().size())
        return;

    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    amrex::MultiFab* ac_ptr;
    if (this->OnSameGrids(lev, acceleration))
    {
        ac_ptr = &acceleration;
    }
    else
    {
        ac_ptr = new amrex::MultiFab(this->m_gdb->ParticleBoxArray(lev),
					 this->m_gdb->ParticleDistributionMap(lev),
					 acceleration.nComp(),acceleration.nGrow());
        for (amrex::MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary(periodic);
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 1;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_agn_particles(&Np, particles.data(),
                                (*ac_ptr)[pti].dataPtr(),
                                ac_box.loVect(), ac_box.hiVect(),
                                plo,dx,dt,a_old,a_half,&do_move);
        }
    }

    if (ac_ptr != &acceleration) delete ac_ptr;
    
    ParticleLevel&    pmap          = this->GetParticles(lev);
    if (lev > 0 && sub_cycle)
    {
        amrex::ParticleLocData pld; 
        for (auto& kv : pmap) {
            AoS&  pbox       = kv.second.GetArrayOfStructs();
            const int   n    = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                if (p.id() <= 0) continue;

                // Move the particle to the proper ghost cell. 
                //      and remove any *ghost* particles that have gone too far
                // Note that this should only negate ghost particles, not real particles.
                if (!this->Where(p, pld, lev, lev, where_width))
                {
                    // Assert that the particle being removed is a ghost particle;
                    // the ghost particle is no longer in relevant ghost cells for this grid.
                    if (p.id() == amrex::GhostParticleID)
                    {
                        p.id() = -1;
                    }
                    else
                    {
                        std::cout << "Oops -- removing particle " << p.id() << std::endl;
                        amrex::Error("Trying to get rid of a non-ghost particle in moveKickDrift");
                    }
                }
            }
        }
    }
}

void
AGNParticleContainer::moveKick (MultiFab&       acceleration,
                                int             lev,
                                Real            dt,
                                Real            a_new,
                                Real            a_half) 
{
    BL_PROFILE("AGNParticleContainer::moveKick()");

    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    MultiFab* ac_ptr;
    if (OnSameGrids(lev,acceleration))
    {
        ac_ptr = &acceleration;
    }
    else 
    {
        ac_ptr = new MultiFab(ParticleBoxArray(lev),
				  ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_ptr); mfi.isValid(); ++mfi)
            ac_ptr->setVal(0.);
        ac_ptr->copy(acceleration,0,0,acceleration.nComp());
        ac_ptr->FillBoundary(periodic);
    }

    const Real* plo = Geom(lev).ProbLo();

    int do_move = 0;

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        if (Np > 0)
        {
           const Box& ac_box = (*ac_ptr)[pti].box();

           update_agn_particles(&Np, particles.data(),
                                (*ac_ptr)[pti].dataPtr(),
                                ac_box.loVect(), ac_box.hiVect(),
                                plo,dx,dt,a_half,a_new,&do_move);
        }
    }
    
    if (ac_ptr != &acceleration) delete ac_ptr;
}

void AGNParticleContainer::ComputeOverlap(int lev)
{
    BL_PROFILE("AGNParticleContainer::ComputeOverlap()");
    Vector<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = neighbors[index].size() / pdata_size;

        nyx_compute_overlap(&Np, particles.data(), 
                            &Ng, neighbors[index].dataPtr(), dx);

    }
}

void AGNParticleContainer::Merge(int lev)
{
    BL_PROFILE("AGNParticleContainer::Merge()");
    Vector<int> my_id;

    const Real* dx = Geom(lev).CellSize();

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        PairIndex index(pti.index(), pti.LocalTileIndex());
        int Ng = neighbors[index].size() / pdata_size;

        agn_merge_particles(&Np, particles.data(), 
                            &Ng, neighbors[index].dataPtr(), dx);
    }
}

void AGNParticleContainer::ComputeParticleVelocity(int lev,
                                                   amrex::MultiFab& state_old, 
                                                   amrex::MultiFab& state_new,
                                                   int add_energy)
{
    BL_PROFILE("AGNParticleContainer::ComputeParticleVelocity()");
    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    state_old.FillBoundary(periodic);
    state_new.FillBoundary(periodic);

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const Box& soldbox = state_old[pti].box();
        const Box& snewbox = state_new[pti].box();

        agn_particle_velocity(&Np, particles.data(), 
                              state_old[pti].dataPtr(), 
                              soldbox.loVect(), soldbox.hiVect(),
                              state_new[pti].dataPtr(), 
                              snewbox.loVect(), snewbox.hiVect(),
                              dx, &add_energy);
    }
}

void AGNParticleContainer::AccreteMass(int lev,
                                       amrex::MultiFab& state,
                                       amrex::MultiFab& density_lost,
                                       amrex::Real dt)
{
    BL_PROFILE("AGNParticleContainer::AccreteMass()");
    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    state.FillBoundary(periodic);

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const Box& sbox = state[pti].box();

        agn_accrete_mass(&Np, particles.data(),
                         state[pti].dataPtr(),
                         density_lost[pti].dataPtr(),
                         sbox.loVect(), sbox.hiVect(),
                         &dt, dx);
    }
}

void AGNParticleContainer::ReleaseEnergy(int lev, amrex::MultiFab& state, amrex::MultiFab& D_new, amrex::Real a)
{
    BL_PROFILE("AGNParticleContainer::ReleaseEnergy()");
    const Real* dx = Geom(lev).CellSize();
    const Periodicity& periodic = Geom(lev).periodicity();

    state.FillBoundary(periodic);
    D_new.FillBoundary(periodic);

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        AoS& particles = pti.GetArrayOfStructs();
        int Np = particles.size();

        const Box& sbox = state[pti].box();
        const Box& Dbox = D_new[pti].box();
        agn_release_energy(&Np, particles.data(), 
                           state[pti].dataPtr(), 
                           sbox.loVect(), sbox.hiVect(),
                           D_new[pti].dataPtr(),
                           Dbox.loVect(), Dbox.hiVect(),
                           &a, dx); 
    }
}

void AGNParticleContainer::writeAllAtLevel(int lev)
{
  BL_PROFILE("AGNParticleContainer::writeAllAtLevel()");
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& particles = pti.GetArrayOfStructs();
      size_t Np = pti.numParticles();
      Print() << "There are " << Np  << " AGN particles at level " << lev
              << " in grid " << pti.index() << std::endl;
      for (unsigned i = 0; i < Np; ++i)
        {
          const ParticleType& p = particles[i];
          const IntVect& iv = Index(p, lev);

          int id = p.id();
          int cpu = p.cpu();
          RealVect xyz(p.pos(0), p.pos(1), p.pos(2));
          Real mass = p.rdata(0);
          RealVect uvw(p.rdata(1), p.rdata(2), p.rdata(3));
          Real energy = p.rdata(4);
          Real mdot = p.rdata(5);

          Print() << "[" << i << "]: id " << id << " cpu " << cpu
                  << " mass " << mass
                  << " index " << iv
                  << " position " << xyz
                  << " velocity " << uvw
                  << " energy " << energy
                  << " mdot " << mdot
                  << endl;
        }
    }
}
