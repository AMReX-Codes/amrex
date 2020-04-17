#include <iomanip>

#include "ElectrostaticParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

#include "electrostatic_pic_F.H"

using namespace amrex;

void ElectrostaticParticleContainer::InitParticles() {

    if ( ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber() ) {

        ParticleType p;

        p.id()   = ParticleType::NextID();
        p.cpu()  = ParallelDescriptor::MyProc(); 
        
        p.pos(0) = -2.5e-6; 
        p.pos(1) =  0.0;
#if BL_SPACEDIM == 3
        p.pos(2) =  0.0;
#endif
        
        std::array<Real,PIdx::nattribs> attribs;
        attribs[PIdx::w]  = 1.0;
        attribs[PIdx::vx] = 0.0;
        attribs[PIdx::vy] = 0.0;
#if BL_SPACEDIM == 3
        attribs[PIdx::vz] = 0.0;
#endif
        attribs[PIdx::Ex] = 0.0;
        attribs[PIdx::Ey] = 0.0;
#if BL_SPACEDIM == 3
        attribs[PIdx::Ez] = 0.0;
#endif
        
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

    // each level deposits it's own particles
    const int ng = rho[0]->nGrow();
    for (int lev = 0; lev < num_levels; ++lev) {       

        rho[lev]->setVal(0.0, ng);

        const auto& gm = m_gdb->Geom(lev);
        const auto& ba = m_gdb->ParticleBoxArray(lev);
    
        const Real* dx  = gm.CellSize();
        const Real* plo = gm.ProbLo();
        BoxArray nba = ba;
        nba.surroundingNodes();
    
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Box& box = nba[pti];
            
            auto& wp = pti.GetAttribs(PIdx::w);
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const Long np  = pti.numParticles();
            
            FArrayBox& rhofab = (*rho[lev])[pti];
            
            deposit_cic(particles.data(), nstride, np,
                        wp.data(), &this->charge,
                        rhofab.dataPtr(), box.loVect(), box.hiVect(), 
                        plo, dx, &ng);
        }

        rho[lev]->SumBoundary(gm.periodicity());
    }

    // now we average down fine to crse
    std::unique_ptr<MultiFab> crse;
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        const BoxArray& fine_BA = rho[lev+1]->boxArray();
        const DistributionMapping& fine_dm = rho[lev+1]->DistributionMap();
        BoxArray coarsened_fine_BA = fine_BA;
        coarsened_fine_BA.coarsen(m_gdb->refRatio(lev));
        
        MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, 1, 0);
        coarsened_fine_data.setVal(0.0);
        
        IntVect ratio(D_DECL(2, 2, 2));  // FIXME
        
        for (MFIter mfi(coarsened_fine_data); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const Box& crse_box = coarsened_fine_data[mfi].box();
            const Box& fine_box = (*rho[lev+1])[mfi].box();
            sum_fine_to_crse_nodal(bx.loVect(), bx.hiVect(), ratio.getVect(),
                                   coarsened_fine_data[mfi].dataPtr(), crse_box.loVect(), crse_box.hiVect(),
                                   (*rho[lev+1])[mfi].dataPtr(), fine_box.loVect(), fine_box.hiVect());
        }
        
        rho[lev]->copy(coarsened_fine_data, m_gdb->Geom(lev).periodicity(), FabArrayBase::ADD);
    }
    
    for (int lev = 0; lev < num_levels; ++lev) {
        rho[lev]->mult(-1.0/PhysConst::ep0, ng);
    }
}

void
ElectrostaticParticleContainer::
FieldGather(const VectorMeshData& E,
            const Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks) {

    const int num_levels = E.size();
    const int ng = E[0][0]->nGrow();

    if (num_levels == 1) {
        const int lev = 0;
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
            const Long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
#if BL_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
#if BL_SPACEDIM == 3
            Ezp.assign(np,0.0);
#endif

            const FArrayBox& exfab = (*E[lev][0])[pti];
            const FArrayBox& eyfab = (*E[lev][1])[pti];
#if BL_SPACEDIM == 3
            const FArrayBox& ezfab = (*E[lev][2])[pti];
#endif

            interpolate_cic(particles.data(), nstride, np,
                            Exp.data(), Eyp.data(), 
#if BL_SPACEDIM == 3                
                            Ezp.data(),
#endif
                            exfab.dataPtr(), eyfab.dataPtr(), 
#if BL_SPACEDIM == 3
                            ezfab.dataPtr(),
#endif
                            box.loVect(), box.hiVect(), plo, dx, &ng);
        }

        return;
    }

    const BoxArray& fine_BA = E[1][0]->boxArray();
    const DistributionMapping& fine_dm = E[1][0]->DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(IntVect(D_DECL(2,2,2)));

    MultiFab coarse_Ex(coarsened_fine_BA, fine_dm, 1, 1);
    MultiFab coarse_Ey(coarsened_fine_BA, fine_dm, 1, 1);
#if BL_SPACEDIM == 3
    MultiFab coarse_Ez(coarsened_fine_BA, fine_dm, 1, 1);
#endif
    
    coarse_Ex.copy(*E[0][0], 0, 0, 1, 1, 1);
    coarse_Ey.copy(*E[0][1], 0, 0, 1, 1, 1);
#if BL_SPACEDIM == 3
    coarse_Ez.copy(*E[0][2], 0, 0, 1, 1, 1);
#endif

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
            const Long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
#if BL_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            Exp.assign(np,0.0);
            Eyp.assign(np,0.0);
#if BL_SPACEDIM == 3
            Ezp.assign(np,0.0);
#endif

            const FArrayBox& exfab = (*E[lev][0])[pti];
            const FArrayBox& eyfab = (*E[lev][1])[pti];
#if BL_SPACEDIM == 3
            const FArrayBox& ezfab = (*E[lev][2])[pti];
#endif

            if (lev == 0) {
                interpolate_cic(particles.data(), nstride, np,
                                Exp.data(), Eyp.data(), 
#if BL_SPACEDIM == 3                
                Ezp.data(),
#endif
                                exfab.dataPtr(), eyfab.dataPtr(), 
#if BL_SPACEDIM == 3
                                ezfab.dataPtr(),
#endif
                                box.loVect(), box.hiVect(), plo, dx, &ng);                
            } else {
                
                const FArrayBox& exfab_coarse = coarse_Ex[pti];
                const FArrayBox& eyfab_coarse = coarse_Ey[pti];
#if BL_SPACEDIM == 3
                const FArrayBox& ezfab_coarse = coarse_Ez[pti];
#endif                
                const Box& coarse_box = coarsened_fine_BA[pti];
                const Real* coarse_dx = Geom(0).CellSize();
                
                interpolate_cic_two_levels(particles.data(), nstride, np,
                                           Exp.data(), Eyp.data(), 
#if BL_SPACEDIM == 3                    
                                           Ezp.data(),
#endif
                                           exfab.dataPtr(), eyfab.dataPtr(), 
#if BL_SPACEDIM == 3
                                           ezfab.dataPtr(),
#endif
                                           box.loVect(), box.hiVect(), dx, 
                                           exfab_coarse.dataPtr(), eyfab_coarse.dataPtr(),
#if BL_SPACEDIM == 3
                                           ezfab_coarse.dataPtr(),
#endif
                                           (*masks[1])[pti].dataPtr(),
                                           coarse_box.loVect(), coarse_box.hiVect(), coarse_dx,
                                           plo, &ng, &lev);
            }
        }
    }
};

void 
ElectrostaticParticleContainer::
Evolve(const VectorMeshData& E, ScalarMeshData& rho, const Real& dt) {
    
    const int num_levels = E.size();

    for (int lev = 0; lev < num_levels; ++lev) {
        
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
    
        BL_ASSERT(OnSameGrids(lev, *rho[lev]));
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            
            // Particle structs
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;           
            const Long np  = pti.numParticles();
            
            // Particle attributes
            auto& attribs = pti.GetAttribs();
            auto& vxp = attribs[PIdx::vx];
            auto& vyp = attribs[PIdx::vy];

#if BL_SPACEDIM == 3
            auto& vzp = attribs[PIdx::vz];
#endif

            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];

#if BL_SPACEDIM == 3
            auto& Ezp = attribs[PIdx::Ez];
#endif
            
            //
            // Particle Push
            //
            push_leapfrog(particles.data(), nstride, np,
                          vxp.data(), vyp.data(), 
#if BL_SPACEDIM == 3
                vzp.data(),
#endif
                          Exp.data(), Eyp.data(), 
#if BL_SPACEDIM == 3
                          Ezp.data(),
#endif
                          &this->charge, &this->mass, &dt,
                          prob_domain.lo(), prob_domain.hi());
        }
    }
}

void ElectrostaticParticleContainer::pushX(const Real& dt) {
    for (int lev = 0; lev <= finestLevel(); ++lev) {    
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const Long np  = pti.numParticles();
            
            auto& attribs = pti.GetAttribs();
            auto& vxp = attribs[PIdx::vx];
            auto& vyp = attribs[PIdx::vy];
#if BL_SPACEDIM == 3
            auto& vzp = attribs[PIdx::vz];
#endif            
            push_leapfrog_positions(particles.data(), nstride, np,
                                    vxp.data(), vyp.data(), 
#if BL_SPACEDIM == 3
                                    vzp.data(), 
#endif
                                    &dt, prob_domain.lo(), prob_domain.hi());

        }
    }
}

void ElectrostaticParticleContainer::writeParticles(int n) {
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}
