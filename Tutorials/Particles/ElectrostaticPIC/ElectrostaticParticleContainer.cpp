#include "ElectrostaticParticleContainer.H"

#include "electrostatic_pic_F.H"

using namespace amrex;

void ElectrostaticParticleContainer::InitParticles() {

    if ( ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber() ) {

        ParticleType p;

        p.id()   = ParticleType::NextID();
        p.cpu()  = ParallelDescriptor::MyProc(); 
        
        p.pos(0) = -5.0e-6; 
        p.pos(1) =  0.0;
        p.pos(2) =  0.0;
        
        std::array<Real,PIdx::nattribs> attribs;
        attribs[PIdx::w] = 1.0;
        attribs[PIdx::vx] = 0.0;
        attribs[PIdx::vy] = 0.0;
        attribs[PIdx::vz] = 0.0;
        attribs[PIdx::Ex] = 0.0;
        attribs[PIdx::Ey] = 0.0;
        attribs[PIdx::Ez] = 0.0;
        
        // Add to level 0, grid 0, and tile 0
        // Redistribute() will move it to the proper place.
        std::pair<int,int> key {0,0};
        auto& particle_tile = GetParticles(0)[key];
        
        particle_tile.push_back(p);
        particle_tile.push_back_real(attribs);

    }

    Redistribute();
}

void
ElectrostaticParticleContainer::DepositCharge(ScalarMeshData& rho) {
    
    int num_levels = rho.size();
    int finest_level = num_levels - 1;

    // info for coarse/fine interpolation
    PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    NodeBilinear mapper;
    
    // temporary MF with zero ghost cells
    Array<std::unique_ptr<MultiFab> > tmp(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        const BoxArray& ba = rho[lev]->boxArray();
        const DistributionMapping& dm = rho[lev]->DistributionMap();
        tmp[lev].reset(new MultiFab(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }

    const int ng = rho[0]->nGrow();
    for (int lev = 0; lev < num_levels; ++lev) {       

        rho[lev]->setVal(0.0, ng);

        // first deposit on this level
        const auto& gm = m_gdb->Geom(lev);
        const auto& ba = m_gdb->ParticleBoxArray(lev);
        const auto& dm = m_gdb->DistributionMap(lev);
    
        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();
        BoxArray nba = ba;
        nba.surroundingNodes();
    
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Box& box = nba[pti];
            
            auto& wp = pti.GetAttribs(PIdx::w);
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();
            
            FArrayBox& rhofab = (*rho[lev])[pti];
            
            deposit_cic(particles.data(), nstride, np,
                        wp.data(), &this->charge,
                        rhofab.dataPtr(), box.loVect(), box.hiVect(), 
                        plo, dx, &ng);
        }
        
        rho[lev]->SumBoundary(gm.periodicity());
        
        // handle coarse particles that deposit some of their mass onto fine
        if (lev < finest_level) {
            amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *rho[lev], 0, 0, 1, 
                                         m_gdb->Geom(lev), m_gdb->Geom(lev+1),
                                         cphysbc, fphysbc,
                                         m_gdb->refRatio(lev), &mapper, bcs);
        }

        // handle fine particles that deposit some of their mass onto coarse
        // Note - this will double count the mass on the coarse level in 
        // regions covered by the fine level, but this will be corrected
        // below in the call to average_down_nodal.
        //        if (lev > 0) {
        //            amrex::sum_fine_to_coarse(*rho[lev], *rho[lev-1], 0, 1, 
        //                                      m_gdb->refRatio(lev-1), m_gdb->Geom(lev-1), m_gdb->Geom(lev));
        //        }

        rho[lev]->plus(*tmp[lev], 0, 1, 0);     
    }

    std::unique_ptr<MultiFab> crse;
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        BoxArray cba = rho[lev+1]->boxArray();
        const DistributionMapping& fdm = rho[lev+1]->DistributionMap();
        cba.coarsen(m_gdb->refRatio(lev));
        crse.reset(new MultiFab(cba, fdm, 1, 0));
        amrex::average_down_nodal(*rho[lev+1], *crse, m_gdb->refRatio(lev));
        rho[lev]->copy(*crse, m_gdb->Geom(lev).periodicity());
    }

    for (int lev = 0; lev < num_levels; ++lev) {
        rho[lev]->mult(-1.0/PhysConst::ep0, 1);
    }
}

void
ElectrostaticParticleContainer::
FieldGather(const VectorMeshData& E) {

    const int num_levels = E.size();
    const int ng = E[0][0]->nGrow();

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = m_gdb->Geom(lev);
        const auto& ba = m_gdb->ParticleBoxArray(lev);

        BoxArray nba = ba;
        nba.surroundingNodes();

        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();

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
                            box.loVect(), box.hiVect(), plo, dx, &ng);
        }
    }
};

void 
ElectrostaticParticleContainer::
Evolve(const VectorMeshData& E, ScalarMeshData& rho, const Real& dt) {
    
    const int num_levels = E.size();
    const int ng = E[0][0]->nGrow();

    for (int lev = 0; lev < num_levels; ++lev) {
        
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
        const auto& ba = m_gdb->ParticleBoxArray(lev);
        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();
    
        BoxArray nba = ba;
        nba.surroundingNodes();
    
        const int ng = 1;
        BL_ASSERT(OnSameGrids(lev, *rho[lev]));
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
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
            
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
            Ezp.assign(np,0.0);
            
            //
            // Field Gather
            //
            interpolate_cic(particles.data(), nstride, np, 
                            Exp.data(), Eyp.data(), Ezp.data(),
                            exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
                            box.loVect(), box.hiVect(), plo, dx, &ng);
            
            //
            // Particle Push
            //
            push_leapfrog(particles.data(), nstride, np,
                          vxp.data(), vyp.data(), vzp.data(),
                          Exp.data(), Eyp.data(), Ezp.data(),
                          &this->charge, &this->mass, &dt,
                          prob_domain.lo(), prob_domain.hi());
        }
    }
}

void ElectrostaticParticleContainer::pushX(const Real& dt) {

    int finest_level = finestLevel();
    for (int lev = 0; lev <= finestLevel(); ++lev) {    
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np  = pti.numParticles();
            
            auto& attribs = pti.GetAttribs();
            auto& vxp = attribs[PIdx::vx];
            auto& vyp = attribs[PIdx::vy];
            auto& vzp = attribs[PIdx::vz];
            
            push_leapfrog_positions(particles.data(), nstride, np,
                                    vxp.data(), vyp.data(), vzp.data(), 
                                    &dt, prob_domain.lo(), prob_domain.hi());

        }
    }
}

void ElectrostaticParticleContainer::writeParticles(int n) {
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}
