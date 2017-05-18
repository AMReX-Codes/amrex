#include "ElectrostaticParticleContainer.H"

#include "electrostatic_pic_F.H"

using namespace amrex;

void ElectrostaticParticleContainer::InitParticles() {

    ParticleType p;

    p.id()   = ParticleType::NextID();
    p.cpu()  = ParallelDescriptor::MyProc(); 
    
    p.pos(0) = -3.0e-7; 
    p.pos(1) =  0.0;
    p.pos(2) =  0.0;

    std::array<Real,PIdx::nattribs> attribs;
    attribs[0] = 1.0;
    attribs[1] = 0.0;
    attribs[2] = 0.0;
    attribs[3] = 0.0;

    // Add to level 0, grid 0, and tile 0
    // Redistribute() will move it to the proper place.
    std::pair<int,int> key {0,0};
    auto& particle_tile = GetParticles(0)[key];

    particle_tile.push_back(p);
    particle_tile.push_back_real(attribs);

    Redistribute();
}

void
ElectrostaticParticleContainer::DepositCharge(ScalarMeshData& rho) {
    
    const int lev = 0;

    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;
    nba.surroundingNodes();
    
    const Real* dx  = gm.CellSize();
    const Real* plo = gm.ProbLo();
    const int ng = 1;
    
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const Box& box = nba[pti];
        
        auto& wp = pti.GetAttribs(PIdx::w);
        const auto& particles = pti.GetArrayOfStructs();
        int nstride = particles.dataShape().first;
        const long np  = pti.numParticles();
        
        FArrayBox& rhofab = (*rho[0])[pti];
        
        deposit_cic(particles.data(), nstride, np,
                    wp.data(), &this->charge,
                    rhofab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);
    }
    
    rho[0]->SumBoundary(gm.periodicity());
}

void
ElectrostaticParticleContainer::
FieldGather(const VectorMeshData& E) {
    const int lev = 0;
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);

    BoxArray nba = ba;
    nba.surroundingNodes();

    const Real* dx  = gm.CellSize();
    const Real* plo = gm.ProbLo();
    const int ng = 1;

    BL_ASSERT(OnSameGrids(lev, *E[lev][0]));

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        const Box& box = nba[pti];

        const auto& particles = pti.GetArrayOfStructs();
        int nstride = particles.dataShape().first;
        const long np  = pti.numParticles();

        auto& attribs = pti.GetAttribs();
        auto& Exp = attribs[PIdx::Ex];
        auto& Eyp = attribs[PIdx::Ey];
        auto& Ezp = attribs[PIdx::Ez];

        const FArrayBox& exfab = (*E[lev][0])[pti];
        const FArrayBox& eyfab = (*E[lev][1])[pti];
        const FArrayBox& ezfab = (*E[lev][2])[pti];

        Exp.assign(np,0.0);
        Eyp.assign(np,0.0);
        Ezp.assign(np,0.0);

        interpolate_cic(particles.data(), nstride, np,
                        Exp.data(), Eyp.data(), Ezp.data(),
                        exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
                        box.loVect(), box.hiVect(), plo, dx);
    }
};

void 
ElectrostaticParticleContainer::
Evolve(const VectorMeshData& E, ScalarMeshData& rho, const Real& dt) {
    
    const int lev = 0;

    const auto& gm = m_gdb->Geom(lev);
    const RealBox& prob_domain = gm.ProbDomain();
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const Real* dx  = gm.CellSize();
    const Real* plo = gm.ProbLo();
    const int ng = 1;
    
    BoxArray nba = ba;
    nba.surroundingNodes();
    
    BL_ASSERT(OnSameGrids(lev, *rho[lev]));
    {
	for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
	{
	    const Box& box = nba[pti];
            
            // Particle structs
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;           
            const long np  = pti.numParticles();

            // Particle attributes
            auto& attribs = pti.GetAttribs();
            auto&  wp = attribs[PIdx::w];
            auto& vxp = attribs[PIdx::vx];
            auto& vyp = attribs[PIdx::vy];
            auto& vzp = attribs[PIdx::vz];
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
            auto& Ezp = attribs[PIdx::Ez];

	    // Data on the grid
	    const FArrayBox& exfab  = (*E[lev][0])[pti];
	    const FArrayBox& eyfab  = (*E[lev][1])[pti];
	    const FArrayBox& ezfab  = (*E[lev][2])[pti];
	    FArrayBox&       rhofab = (*rho[lev])[pti];

	    Exp.assign(np,0.0);
	    Eyp.assign(np,0.0);
	    Ezp.assign(np,0.0);

	    //
	    // Field Gather
	    //
            interpolate_cic(particles.data(), nstride, np, 
                            Exp.data(), Eyp.data(), Ezp.data(),
                            exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
                            box.loVect(), box.hiVect(), plo, dx);

	    //
	    // Particle Push
	    //
            push_leapfrog(particles.data(), nstride, np,
                          vxp.data(), vyp.data(), vzp.data(),
                          Exp.data(), Eyp.data(), Ezp.data(),
                          &this->charge, &this->mass, &dt,
                          prob_domain.lo(), prob_domain.hi());

	    //
	    // Charge Deposition
	    // xxxxx this part needs to be thread safe if we have OpenMP over tiles
	    //
            deposit_cic(particles.data(), nstride, np,
                        wp.data(), &this->charge,
                        rhofab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);           
	}
    }
}

void ElectrostaticParticleContainer::pushX(const Real& dt) {
    const int lev = 0;
    BL_PROFILE("WPC::PushXES()");
    
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        int nstride = particles.dataShape().first;
        const long np  = pti.numParticles();
        
        auto& attribs = pti.GetAttribs();
        auto& vxp = attribs[PIdx::vx];
        auto& vyp = attribs[PIdx::vy];
        auto& vzp = attribs[PIdx::vz];
        
        push_leapfrog_positions(particles.data(), nstride, np,
                                vxp.data(), vyp.data(), vzp.data(), &dt);
    }
}

void ElectrostaticParticleContainer::writeParticles(int n) {
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}
