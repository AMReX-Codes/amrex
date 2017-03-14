#include <limits>
#include <sstream>

#include <ParticleContainer.H>
#include <AMReX_ParmParse.H>
#include <WarpX_f.H>
#include <WarpX.H>

using namespace amrex;

PhysicalParticleContainer::PhysicalParticleContainer (AmrCore* amr_core, int ispecies)
    : WarpXParticleContainer(amr_core, ispecies)
{
    std::stringstream ss;
    ss << "species" << ispecies; 
    ParmParse pp(ss.str());

    std::string plasma_profile_s;
    pp.get("profile", plasma_profile_s);
    std::transform(plasma_profile_s.begin(), plasma_profile_s.end(), plasma_profile_s.begin(), ::tolower);
    if (plasma_profile_s == "constant") {
	plasma_injector.reset(new ConstantPlasmaInjector);
    } 
    
    else if (plasma_profile_s == "double_ramp") {
	plasma_injector.reset(new DoubleRampPlasmaInjector);
    }
    
    else if (plasma_profile_s == "custom") {
	plasma_injector.reset(new CustomPlasmaInjector);
    }
    
    else {
	amrex::Abort("Unknown plasma injector type");
    }
}

void
PhysicalParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because grids have not
    // been built when constructor was called.
    reserveData();
    resizeData();
}

void
PhysicalParticleContainer::FieldGather (int lev,
                                        const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                        const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    const Geometry& gm  = Geom(lev);

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif

#if (BL_SPACEDIM == 3)
    long ngx_eb = Ex.nGrow();
    long ngy_eb = ngx_eb;
    long ngz_eb = ngx_eb;
#elif (BL_SPACEDIM == 2)
    long ngx_eb = Ex.nGrow();
    long ngy_eb = 0;
    long ngz_eb = ngx_eb;
#endif

    BL_ASSERT(OnSameGrids(lev,Ex));

    {
	Array<Real> xp, yp, zp;

	for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
	{
	    const Box& box = pti.validbox();

            auto& attribs = pti.GetAttribs();

            auto&  wp = attribs[PIdx::w];
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
#if (BL_SPACEDIM == 3)
            pti.GetPosition(xp, yp, zp);
#elif (BL_SPACEDIM == 2)
            pti.GetPosition(xp, zp);
            yp.resize(np, std::numeric_limits<Real>::quiet_NaN());
#endif

#if (BL_SPACEDIM == 3)
	    long nx = box.length(0);
	    long ny = box.length(1);
	    long nz = box.length(2);
#elif (BL_SPACEDIM == 2)
	    long nx = box.length(0);
	    long ny = 0;
	    long nz = box.length(1);
#endif
	    RealBox grid_box = RealBox( box, gm.CellSize(), gm.ProbLo() );
#if (BL_SPACEDIM == 3)
	    const Real* xyzmin = grid_box.lo();
#elif (BL_SPACEDIM == 2)
	    Real xyzmin[3] = { grid_box.lo(0), std::numeric_limits<Real>::quiet_NaN(), grid_box.lo(1) };
#endif

	    //
	    // Field Gather
	    //
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
	    warpx_geteb_energy_conserving(&np, xp.data(), yp.data(), zp.data(),
					  Exp.data(),Eyp.data(),Ezp.data(),
					  Bxp.data(),Byp.data(),Bzp.data(),
					  &xyzmin[0], &xyzmin[1], &xyzmin[2],
					  &dx[0], &dx[1], &dx[2],
					  &nx, &ny, &nz, &ngx_eb, &ngy_eb, &ngz_eb,
					  &WarpX::nox, &WarpX::noy, &WarpX::noz,
					  exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
					  bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
					  &ll4symtry, &l_lower_order_in_v,
                                          &lvect_fieldgathe,
		                          &WarpX::field_gathering_algo);
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

    const Geometry& gm  = Geom(lev);

#if (BL_SPACEDIM == 3)
    const Real* dx = gm.CellSize();
#elif (BL_SPACEDIM == 2)
    Real dx[3] = { gm.CellSize(0), std::numeric_limits<Real>::quiet_NaN(), gm.CellSize(1) };
#endif

#if (BL_SPACEDIM == 3)
    long ngx_eb = Ex.nGrow();
    long ngy_eb = ngx_eb;
    long ngz_eb = ngx_eb;
    long ngx_j  = jx.nGrow();
    long ngy_j  = ngx_j;
    long ngz_j  = ngx_j;
#elif (BL_SPACEDIM == 2)
    long ngx_eb = Ex.nGrow();
    long ngy_eb = 0;
    long ngz_eb = ngx_eb;
    long ngx_j  = jx.nGrow();;
    long ngy_j  = 0;
    long ngz_j  = ngx_j;
#endif

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
#if (BL_SPACEDIM == 3)
            pti.GetPosition(xp, yp, zp);
#elif (BL_SPACEDIM == 2)
            pti.GetPosition(xp, zp);
            yp.resize(np, std::numeric_limits<Real>::quiet_NaN());
#endif
	    BL_PROFILE_VAR_STOP(blp_copy);

#if (BL_SPACEDIM == 3)
	    long nx = box.length(0);
	    long ny = box.length(1);
	    long nz = box.length(2);
#elif (BL_SPACEDIM == 2)
	    long nx = box.length(0);
	    long ny = 0;
	    long nz = box.length(1);
#endif
	    RealBox grid_box = RealBox( box, gm.CellSize(), gm.ProbLo() );
#if (BL_SPACEDIM == 3)
	    const Real* xyzmin = grid_box.lo();
#elif (BL_SPACEDIM == 2)
	    Real xyzmin[3] = { grid_box.lo(0), std::numeric_limits<Real>::quiet_NaN(), grid_box.lo(1) };
#endif

	    //
	    // Field Gather
	    //
	    const int ll4symtry          = false;
	    const int l_lower_order_in_v = true;
            long lvect_fieldgathe = 64;
	    BL_PROFILE_VAR_START(blp_pxr_fg);
	    warpx_geteb_energy_conserving(&np, xp.data(), yp.data(), zp.data(),
					  Exp.data(),Eyp.data(),Ezp.data(),
					  Bxp.data(),Byp.data(),Bzp.data(),
					  &xyzmin[0], &xyzmin[1], &xyzmin[2],
					  &dx[0], &dx[1], &dx[2],
					  &nx, &ny, &nz, &ngx_eb, &ngy_eb, &ngz_eb,
					  &WarpX::nox, &WarpX::noy, &WarpX::noz,
					  exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
					  bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
					  &ll4symtry, &l_lower_order_in_v,
                                          &lvect_fieldgathe,
		                          &WarpX::field_gathering_algo);
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
	    warpx_current_deposition(jxfab.dataPtr(), jyfab.dataPtr(), jzfab.dataPtr(),
				     &np, xp.data(), yp.data(), zp.data(),
				     uxp.data(), uyp.data(), uzp.data(),
				     giv.data(), wp.data(), &this->charge,
				     &xyzmin[0], &xyzmin[1], &xyzmin[2],
				     &dt, &dx[0], &dx[1], &dx[2], &nx, &ny, &nz,
				     &ngx_j, &ngy_j, &ngz_j,
                                     &WarpX::nox,&WarpX::noy,&WarpX::noz,
				     &lvect,&WarpX::current_deposition_algo);
	    BL_PROFILE_VAR_STOP(blp_pxr_cd);

	    //
	    // copy particle data back
	    //
	    BL_PROFILE_VAR_START(blp_copy);
#if (BL_SPACEDIM == 3)
            pti.SetPosition(xp, yp, zp);
#elif (BL_SPACEDIM == 2)
            pti.SetPosition(xp, zp);
#endif
            BL_PROFILE_VAR_STOP(blp_copy);
	}
    }
}
