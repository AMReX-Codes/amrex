
#include <TracerParticles.H>

//
// Uses midpoint method to advance particles using umac.
//
void
TracerParticleContainer::AdvectWithUmac (MultiFab* umac, int lev, Real dt)
{
    BL_PROFILE("ParticleContainer<NR,NI>::AdvectWithUmac()");
    BL_ASSERT(OK(true, lev, umac[0].nGrow()-1));
    BL_ASSERT(lev >= 0 && lev < m_particles.size());

    D_TERM(BL_ASSERT(umac[0].nGrow() >= 1);,
           BL_ASSERT(umac[1].nGrow() >= 1);,
           BL_ASSERT(umac[2].nGrow() >= 1););

    D_TERM(BL_ASSERT(!umac[0].contains_nan());,
           BL_ASSERT(!umac[1].contains_nan());,
           BL_ASSERT(!umac[2].contains_nan()););

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();
    const Real*     plo      = geom.ProbLo();

    PArray<MultiFab> umac_pointer;
    // We assume that if umac[0]'s boxArray matches then the others will too...
    if (OnSameGrids(lev, umac[0]))
    {
	umac_pointer.resize(BL_SPACEDIM, PArrayNoManage);
        for (int i = 0; i < BL_SPACEDIM; i++)
	    umac_pointer.set(i, &umac[i]);
    }
    else
    {
	umac_pointer.resize(BL_SPACEDIM, PArrayManage);
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
	    int ng = umac[i].nGrow();
	    
	    umac_pointer.set(i, new MultiFab(m_gdb->ParticleBoxArray(lev),
					     umac[i].nComp(),
					     ng,
					     m_gdb->ParticleDistributionMap(lev),
					     Fab_allocate,
					     IntVect::TheDimensionVector(i)));
	    umac_pointer[i].copy(umac[i],0,0,umac[i].nComp(),ng,ng);
        }
    }

    for (int ipass = 0; ipass < 2; ipass++)
    {
        PMap& pmap = m_particles[lev];

        for (auto& kv : pmap)
        {
            const int grid = kv.first;
            PBox&     pbox = kv.second;
            const int n    = pbox.size();

            FArrayBox* fab[BL_SPACEDIM] = { D_DECL(&umac_pointer[0][grid],
                                                   &umac_pointer[1][grid],
                                                   &umac_pointer[2][grid]) };

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                
                if (p.m_id <= 0) continue;

                BL_ASSERT(p.m_grid == grid);

                const Real len[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0])/dx[0] + Real(0.5),
                                                       (p.m_pos[1]-plo[1])/dx[1] + Real(0.5),
                                                       (p.m_pos[2]-plo[2])/dx[2] + Real(0.5)) };

                const IntVect cell(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

                const Real frac[BL_SPACEDIM] = { D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    IntVect ecell = cell;

                    ecell[d] = p.m_cell[d] + 1;

                    Real efrac[BL_SPACEDIM] = { D_DECL(frac[0], frac[1], frac[2]) };

                    efrac[d] = (p.m_pos[d]-plo[d])/dx[d] - p.m_cell[d];

                    for (int j = 0; j < BL_SPACEDIM; j++)
                    {
                        if (efrac[j] > 1) efrac[j] = 1;
                        if (efrac[j] < 0) efrac[j] = 0;
                    }

                    const Real vel = ParticleBase::InterpDoit(*fab[d], ecell, efrac, 0);

                    if (ipass == 0)
                    {
                        //
                        // Save old position and the vel & predict location at dt/2.
                        //
                        p.m_data[d] = p.m_pos[d];
                        p.m_pos[d] += 0.5*dt*vel;
                    }
                    else
                    {
                        //
                        // Update to final time using the orig position and the vel at dt/2.
                        //
                        p.m_pos[d]  = p.m_data[d] + dt*vel;
                        // Save the velocity for use in Timestamp().
			p.m_data[d] = vel;
                    }
                }
                
                ParticleBase::RestrictedWhere(p,m_gdb, umac[0].nGrow()); 
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "TracerParticleContainer::AdvectWithUmac() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}

//
// Uses midpoint method to advance particles using cell-centered velocity
//
void
TracerParticleContainer::AdvectWithUcc (const MultiFab& Ucc, int lev, Real dt)
{
    BL_ASSERT(Ucc.nGrow() > 0);
    BL_ASSERT(OK(true, lev, Ucc.nGrow()-1));
     BL_ASSERT(lev >= 0 && lev < m_particles.size());

    BL_ASSERT(!Ucc.contains_nan());

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    BL_ASSERT(OnSameGrids(lev,Ucc));

    int idx[BL_SPACEDIM] = {D_DECL(0,1,2)};

    for (int ipass = 0; ipass < 2; ipass++)
    {
        PMap& pmap = m_particles[lev];

        for (auto& kv : pmap)
        {
            const int grid = kv.first;
            PBox&     pbox = kv.second;
            const int n    = pbox.size();

	    const FArrayBox& fab = Ucc[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];
                
                if (p.m_id <= 0) continue;

                BL_ASSERT(p.m_grid == grid);

		Real v[BL_SPACEDIM];

		ParticleBase::Interp(p, geom, fab, idx, v, BL_SPACEDIM);

		if (ipass == 0) {
		    //
		    // Save old position and the vel & predict location at dt/2.
		    //
		    for (int d = 0; d < BL_SPACEDIM; d++)
		    {
			p.m_data[d] = p.m_pos[d];
                        p.m_pos[d] += 0.5*dt*v[d];
                    }
		} else {
		    //
		    // Update to final time using the orig position and the vel at dt/2.
		    //
		    for (int d = 0; d < BL_SPACEDIM; d++)
		    {
                        p.m_pos[d]  = p.m_data[d] + dt*v[d];
                        // Save the velocity for use in Timestamp().
			p.m_data[d] = v[d];
                    }
                }
                
                ParticleBase::RestrictedWhere(p,m_gdb, Ucc.nGrow()); 
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "TracerParticleContainer::AdvectWithUcc() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}
