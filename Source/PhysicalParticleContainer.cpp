#include <limits>
#include <sstream>

#include <ParticleContainer.H>
#include <WarpX_f.H>
#include <WarpX.H>

using namespace amrex;

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core, int ispecies,
                                                      const std::string& name)
    : WarpXParticleContainer(amr_core, ispecies),
      species_name(name)
{
    plasma_injector.reset(new PlasmaInjector(species_id, species_name));
    charge = plasma_injector->getCharge();
    mass = plasma_injector->getMass();
}

void PhysicalParticleContainer::InitData() {
    AddParticles(0); // Note - only one level right now
}

void
PhysicalParticleContainer::AddParticles (int lev, Box part_box) {
    BL_PROFILE("PhysicalParticleContainer::AddParticles()");

    if ( not plasma_injector->doInjection() ) return;

    const Geometry& geom = Geom(lev);
    if (!part_box.ok()) part_box = geom.Domain();

    int num_ppc = plasma_injector->num_particles_per_cell;

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    Real scale_fac;

#if BL_SPACEDIM==3
    scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;
#elif BL_SPACEDIM==2
    scale_fac = dx[0]*dx[2]/num_ppc;
#endif

    std::array<Real,PIdx::nattribs> attribs;
    attribs.fill(0.0);
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& tile_box = mfi.tilebox();
        const Box& intersectBox = tile_box & part_box;
        if (!intersectBox.ok()) continue;

        const std::array<Real,3>& tile_corner = WarpX::LowerCorner(intersectBox, lev);

        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();

        const auto& boxlo = intersectBox.smallEnd();
        for (IntVect iv = intersectBox.smallEnd(); iv <= intersectBox.bigEnd(); intersectBox.next(iv))
        {
            for (int i_part=0; i_part<num_ppc;i_part++)
            {
                std::array<Real, 3> r;
                plasma_injector->getPositionUnitBox(r, i_part);
#if ( BL_SPACEDIM == 3 )
                Real x = tile_corner[0] + (iv[0]-boxlo[0] + r[0])*dx[0];
                Real y = tile_corner[1] + (iv[1]-boxlo[1] + r[1])*dx[1];
                Real z = tile_corner[2] + (iv[2]-boxlo[2] + r[2])*dx[2];
#elif ( BL_SPACEDIM == 2 )
                Real x = tile_corner[0] + (iv[0]-boxlo[0] + r[0])*dx[0];
                Real y = 0.;
                Real z = tile_corner[2] + (iv[1]-boxlo[1] + r[2])*dx[2];
#endif
                if (plasma_injector->insideBounds(x, y, z)) {
                    Real weight;
                    std::array<Real, 3> u;
                    plasma_injector->getMomentum(u);
                    weight = plasma_injector->getDensity(x, y, z) * scale_fac;
                    attribs[PIdx::w ] = weight;
                    attribs[PIdx::ux] = u[0];
                    attribs[PIdx::uy] = u[1];
                    attribs[PIdx::uz] = u[2];
                    AddOneParticle(lev, grid_id, tile_id, x, y, z, attribs);
                }
            }
        }
    }
}

void
PhysicalParticleContainer::
FieldGatherES (const amrex::Array<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E)
{

    const int lev = 0;
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);

    BoxArray nba = ba;
    nba.surroundingNodes();

    const Real* dx  = gm.CellSize();
    const Real* plo = gm.ProbLo();
    const int ng = 1;

    BL_ASSERT(OnSameGrids(lev, *E[lev][0]));

    for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti) {
        const Box& box = nba[pti];

        // Particle structs
        const auto& particles = pti.GetArrayOfStructs();
        int nstride = particles.dataShape().first;           
        const long np  = pti.numParticles();

        // Particle attributes
        auto& attribs = pti.GetAttribs();
        auto& Exp = attribs[PIdx::Ex];
        auto& Eyp = attribs[PIdx::Ey];
        auto& Ezp = attribs[PIdx::Ez];
        
        // Data on the grid
        const FArrayBox& exfab = (*E[lev][0])[pti];
        const FArrayBox& eyfab = (*E[lev][1])[pti];
        const FArrayBox& ezfab = (*E[lev][2])[pti];
        
        Exp.assign(np,0.0);
        Eyp.assign(np,0.0);
        Ezp.assign(np,0.0);
        
        //
        // Field Gather
        //
        warpx_interpolate_cic(particles.data(), nstride, np, 
                              Exp.data(), Eyp.data(), Ezp.data(),
                              exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
                              box.loVect(), box.hiVect(), plo, dx);
    }
}

void
PhysicalParticleContainer::FieldGather (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    // WarpX assumes the same number of guard cells for Ex, Ey, Ez, Bx, By, Bz
    long ng = Ex.nGrow();

    BL_ASSERT(OnSameGrids(lev,Ex));

    {
	Array<Real> xp, yp, zp;

	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
	    const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
            auto& Ezp = attribs[PIdx::Ez];
            auto& Bxp = attribs[PIdx::Bx];
            auto& Byp = attribs[PIdx::By];
            auto& Bzp = attribs[PIdx::Bz];

            const long np = pti.numParticles();

	    // Data on the grid
	    const FArrayBox& exfab = Ex[pti];
	    const FArrayBox& eyfab = Ey[pti];
	    const FArrayBox& ezfab = Ez[pti];
	    const FArrayBox& bxfab = Bx[pti];
	    const FArrayBox& byfab = By[pti];
	    const FArrayBox& bzfab = Bz[pti];

	    Exp.assign(np,0.0);
	    Eyp.assign(np,0.0);
	    Ezp.assign(np,0.0);
	    Bxp.assign(np,0.0);
	    Byp.assign(np,0.0);
	    Bzp.assign(np,0.0);

	    //
	    // copy data from particle container to temp arrays
	    //
            pti.GetPosition(xp, yp, zp);

            const std::array<Real,3>& xyzmin = WarpX::LowerCorner(box, lev);

	    //
	    // Field Gather
	    //
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
	    warpx_geteb_energy_conserving(
	       &np, xp.data(), yp.data(), zp.data(),
	       Exp.data(),Eyp.data(),Ezp.data(),
	       Bxp.data(),Byp.data(),Bzp.data(),
	       &xyzmin[0], &xyzmin[1], &xyzmin[2],
	       &dx[0], &dx[1], &dx[2],
	       &WarpX::nox, &WarpX::noy, &WarpX::noz,
	       exfab.dataPtr(), &ng, exfab.length(),
	       eyfab.dataPtr(), &ng, eyfab.length(),
	       ezfab.dataPtr(), &ng, ezfab.length(),
               bxfab.dataPtr(), &ng, bxfab.length(),
	       byfab.dataPtr(), &ng, byfab.length(),
	       bzfab.dataPtr(), &ng, bzfab.length(),
	       &ll4symtry, &l_lower_order_in_v,
	       &lvect_fieldgathe, &WarpX::field_gathering_algo);
        }
    }
}

void
PhysicalParticleContainer::EvolveES (const Array<std::array<std::unique_ptr<MultiFab>, 3> >& E,
                                           Array<std::unique_ptr<MultiFab> >& rho,
                                     Real t, Real dt)
{
    BL_PROFILE("PPC::EvolveES()");

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
	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
	    const Box& box = nba[pti];

            // Particle structs
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;           
            const long np  = pti.numParticles();

            // Particle attribues
            auto& attribs = pti.GetAttribs();
            auto&  wp = attribs[PIdx::w];
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];
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
            warpx_interpolate_cic(particles.data(), nstride, np, 
                                  Exp.data(), Eyp.data(), Ezp.data(),
                                  exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
                                  box.loVect(), box.hiVect(), plo, dx);

	    //
	    // Particle Push
	    //
            warpx_push_leapfrog(particles.data(), nstride, np,
                                uxp.data(), uyp.data(), uzp.data(),
                                Exp.data(), Eyp.data(), Ezp.data(),
                                &this->charge, &this->mass, &dt,
                                prob_domain.lo(), prob_domain.hi());

	    //
	    // Charge Deposition
	    // xxxxx this part needs to be thread safe if we have OpenMP over tiles
	    //
            warpx_deposit_cic(particles.data(), nstride, np,
                              wp.data(), &this->charge,
                              rhofab.dataPtr(), box.loVect(), box.hiVect(), plo, dx);           
	}
    }
}

void
PhysicalParticleContainer::Evolve (int lev,
				   const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
				   const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
				   MultiFab& jx, MultiFab& jy, MultiFab& jz, Real t, Real dt)
{
    BL_PROFILE("PPC::Evolve()");
    BL_PROFILE_VAR_NS("PPC::Evolve::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr_fg);
    BL_PROFILE_VAR_NS("PICSAR::ParticlePush", blp_pxr_pp);
    BL_PROFILE_VAR_NS("PICSAR::CurrentDeposition", blp_pxr_cd);

    const std::array<Real,3>& dx = WarpX::CellSize(lev);

    // WarpX assumes the same number of guard cells for Ex, Ey, Ez, Bx, By, Bz
    long ngE = Ex.nGrow();
    // WarpX assumes the same number of guard cells for Jx, Jy, Jz
    long ngJ = jx.nGrow();

    BL_ASSERT(OnSameGrids(lev,Ex));

    {
	Array<Real> xp, yp, zp, giv;

	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
	    const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto&  wp = attribs[PIdx::w];
            auto& uxp = attribs[PIdx::ux];
            auto& uyp = attribs[PIdx::uy];
            auto& uzp = attribs[PIdx::uz];
            auto& Exp = attribs[PIdx::Ex];
            auto& Eyp = attribs[PIdx::Ey];
            auto& Ezp = attribs[PIdx::Ez];
            auto& Bxp = attribs[PIdx::Bx];
            auto& Byp = attribs[PIdx::By];
            auto& Bzp = attribs[PIdx::Bz];

            const long np = pti.numParticles();

	    // Data on the grid
	    const FArrayBox& exfab = Ex[pti];
	    const FArrayBox& eyfab = Ey[pti];
	    const FArrayBox& ezfab = Ez[pti];
	    const FArrayBox& bxfab = Bx[pti];
	    const FArrayBox& byfab = By[pti];
	    const FArrayBox& bzfab = Bz[pti];
	    FArrayBox&       jxfab = jx[pti];
	    FArrayBox&       jyfab = jy[pti];
	    FArrayBox&       jzfab = jz[pti];

	    Exp.assign(np,0.0);
	    Eyp.assign(np,0.0);
	    Ezp.assign(np,0.0);
	    Bxp.assign(np,0.0);
	    Byp.assign(np,0.0);
	    Bzp.assign(np,0.0);

	    giv.resize(np);

	    //
	    // copy data from particle container to temp arrays
	    //
	    BL_PROFILE_VAR_START(blp_copy);
            pti.GetPosition(xp, yp, zp);
	    BL_PROFILE_VAR_STOP(blp_copy);

            const std::array<Real,3>& xyzmin = WarpX::LowerCorner(box, lev);

	    //
	    // Field Gather
	    //
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
	    BL_PROFILE_VAR_START(blp_pxr_fg);
	    warpx_geteb_energy_conserving(
                &np, xp.data(), yp.data(), zp.data(),
                Exp.data(),Eyp.data(),Ezp.data(),
                Bxp.data(),Byp.data(),Bzp.data(),
                &xyzmin[0], &xyzmin[1], &xyzmin[2],
                &dx[0], &dx[1], &dx[2],
                &WarpX::nox, &WarpX::noy, &WarpX::noz,
                exfab.dataPtr(), &ngE, exfab.length(),
                eyfab.dataPtr(), &ngE, eyfab.length(),
                ezfab.dataPtr(), &ngE, ezfab.length(),
                bxfab.dataPtr(), &ngE, bxfab.length(),
                byfab.dataPtr(), &ngE, byfab.length(),
                bzfab.dataPtr(), &ngE, bzfab.length(),
                &ll4symtry, &l_lower_order_in_v,
                &lvect_fieldgathe, &WarpX::field_gathering_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_fg);

	    //
	    // Particle Push
	    //
	    BL_PROFILE_VAR_START(blp_pxr_pp);
	    warpx_particle_pusher(&np, xp.data(), yp.data(), zp.data(),
				  uxp.data(), uyp.data(), uzp.data(), giv.data(),
				  Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
				  Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
				  &this->charge, &this->mass, &dt,
				  &WarpX::particle_pusher_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_pp);

	    //
	    // Current Deposition
	    // xxxxx this part needs to be thread safe if we have OpenMP over tiles
	    //
	    long lvect = 8;
	    BL_PROFILE_VAR_START(blp_pxr_cd);
	    warpx_current_deposition(
                jxfab.dataPtr(), &ngJ, jxfab.length(),
                jyfab.dataPtr(), &ngJ, jyfab.length(),
                jzfab.dataPtr(), &ngJ, jzfab.length(),
                &np, xp.data(), yp.data(), zp.data(),
                uxp.data(), uyp.data(), uzp.data(),
                giv.data(), wp.data(), &this->charge,
                &xyzmin[0], &xyzmin[1], &xyzmin[2],
                &dt, &dx[0], &dx[1], &dx[2],
                &WarpX::nox,&WarpX::noy,&WarpX::noz,
                &lvect,&WarpX::current_deposition_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_cd);

	    //
	    // copy particle data back
	    //
	    BL_PROFILE_VAR_START(blp_copy);
            pti.SetPosition(xp, yp, zp);
            BL_PROFILE_VAR_STOP(blp_copy);
	}
    }
}
