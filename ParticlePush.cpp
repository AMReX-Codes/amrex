
#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::ParticlePush(Real dt)
{
    BL_PROFILE("MyPC::ParticlePush()");
    BL_PROFILE_VAR_NS("MyPC::ParticlePush::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PICSAR::FieldGather", blp_pxr);

    int   lev  = 0; 
    PMap& pmap = m_particles[lev];

    Real half_dt = 0.5 * dt;

    // Mass
    Real mass = 1.0;

    // Charge
    Real q  = 1.0;
    // xxxxx what's weight? Is it the same as charge?
    // xxxxx no negative charges?

    Array<Real>  xp,  yp,  zp;
    Array<Real> uxp, uyp, uzp;

    // Loop over pmap which loops over the grids containing particles
    for (auto& kv : pmap)
    {
	const int gid = kv.first;
        PBox&     pbx = kv.second;

        const long np = pbx.size();

	auto& pdata = partdata[gid];
	const auto& Exp = pdata[PIdx::Ex];
	const auto& Eyp = pdata[PIdx::Ey];
	const auto& Ezp = pdata[PIdx::Ez];
	const auto& Bxp = pdata[PIdx::Bx];
	const auto& Byp = pdata[PIdx::By];
	const auto& Bzp = pdata[PIdx::Bz];
	auto& giv = pdata[PIdx::gaminv];

	BL_ASSERT(Exp->size() == np);

	giv->resize(np);

	BL_PROFILE_VAR_START(blp_copy);

	// 1D Arrays of particle attributes
	 xp.resize(np);
	 yp.resize(np);
	 zp.resize(np);
	uxp.resize(np);
	uyp.resize(np);
	uzp.resize(np);

	// Loop over particles in that box 
        for (auto i = 0; i < np; ++i)
        {
	    const auto& p = pbx[i];
	    BL_ASSERT(p.m_id > 0);
	    xp[i]  = p.m_pos[0];
	    yp[i]  = p.m_pos[1];
	    zp[i]  = p.m_pos[2];
 	    uxp[i] = p.m_data[PIdx::ux]; 
 	    uyp[i] = p.m_data[PIdx::uy]; 
 	    uzp[i] = p.m_data[PIdx::uz]; 
        }

	BL_PROFILE_VAR_STOP(blp_copy);

	BL_PROFILE_VAR_START(blp_pxr);

        pxr_epush_v(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
		    Exp->dataPtr(), Eyp->dataPtr(), Ezp->dataPtr(),
                    &q, &mass, &half_dt);

        pxr_set_gamma(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), giv->dataPtr());

        pxr_bpush_v(&np, uxp.dataPtr(),uyp.dataPtr(),uzp.dataPtr(), giv->dataPtr(),
		    Bxp->dataPtr(), Byp->dataPtr(), Bzp->dataPtr(),
                    &q, &mass, &dt);

        pxr_epush_v(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
		    Exp->dataPtr(), Eyp->dataPtr(), Ezp->dataPtr(),
                    &q, &mass, &half_dt);

        pxr_set_gamma(&np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), giv->dataPtr());

        pxr_pushxyz(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
		    uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), giv->dataPtr(),
		    &dt);

	BL_PROFILE_VAR_STOP(blp_pxr);

	BL_PROFILE_VAR_START(blp_copy);

	// Loop over particles in that box 
        for (auto i = 0; i < np; ++i)
        {
	    auto& p = pbx[i];
	    BL_ASSERT(p.m_id > 0);
	    p.m_pos[0] = xp[i];
	    p.m_pos[1] = yp[i];
	    p.m_pos[2] = zp[i];
 	    p.m_data[PIdx::ux] = uxp[i];
 	    p.m_data[PIdx::uy] = uyp[i];
 	    p.m_data[PIdx::uz] = uzp[i];
        }

	BL_PROFILE_VAR_STOP(blp_copy);
    }
}
