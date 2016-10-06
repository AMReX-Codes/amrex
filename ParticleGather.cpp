
#include <ParticleContainer.H>
#include <PICSAR_f.H>

void
MyParticleContainer::GatherField(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz)
{
    BL_PROFILE("MyPC::GatherField");
    BL_PROFILE_VAR_NS("MyPC::GatherField::Copy", blp_copy);
    BL_PROFILE_VAR_NS("PXR::FieldGather", blp_pxr);

    const int       lev = 0; 
    const Geometry& gm  = m_gdb->Geom(lev);
    const BoxArray& ba  = Ex.boxArray();
    const Real*     dx  = gm.CellSize();

    BL_ASSERT(OnSameGrids(lev,Ex));

    const PMap& pmap = m_particles[lev];

    Array<Real> xp, yp, zp;

    // Loop over the grids containing particles
    for (auto& kv : pmap)
    {
        const  int  gid = kv.first;
        const PBox& pbx = kv.second;

	const long np = pbx.size();

	auto& pdata = partdata[gid];
	auto& Exp = pdata[PIdx::Ex];
	auto& Eyp = pdata[PIdx::Ey];
	auto& Ezp = pdata[PIdx::Ez];
	auto& Bxp = pdata[PIdx::Bx];
	auto& Byp = pdata[PIdx::By];
	auto& Bzp = pdata[PIdx::Bz];

	// 1D Arrays of particle attributes
	xp.resize(np);
	yp.resize(np);
	zp.resize(np);

	Exp->resize(np);
	Eyp->resize(np);
	Ezp->resize(np);
	Bxp->resize(np);
	Byp->resize(np);
	Bzp->resize(np);

	// Data on the grid
        const FArrayBox& exfab = Ex[gid];
        const FArrayBox& eyfab = Ey[gid];
        const FArrayBox& ezfab = Ez[gid];
        const FArrayBox& bxfab = Bx[gid];
        const FArrayBox& byfab = By[gid];
        const FArrayBox& bzfab = Bz[gid];

	BL_PROFILE_VAR_START(blp_copy);

	// Loop over particles in that box 
        for (auto i = 0; i < np; ++i)
        {
	    const auto& p = pbx[i];
	    BL_ASSERT(p.m_id > 0);
	    xp[i] = p.m_pos[0];
	    yp[i] = p.m_pos[1];
	    zp[i] = p.m_pos[2];
        }

	BL_PROFILE_VAR_STOP(blp_copy);

	const Box& box = BoxLib::enclosedCells(ba[gid]);
	long nx = box.length(0);
	long ny = box.length(1);
	long nz = box.length(2); 

	RealBox grid_box = RealBox( box, dx, gm.ProbLo() );
	const Real* xyzmin = grid_box.lo();

	const long order = 1;
	const long field_gathe_algo = 1;
	const long ng = 1;
	const int ll4symtry          = false;
	const int l_lower_order_in_v = true;
	
	BL_ASSERT(ng == Ex.nGrow());
	BL_ASSERT(ng == Ey.nGrow());
	BL_ASSERT(ng == Ez.nGrow());
	BL_ASSERT(ng == Bx.nGrow());
	BL_ASSERT(ng == By.nGrow());
	BL_ASSERT(ng == Bz.nGrow());

	BL_PROFILE_VAR_START(blp_pxr);

        warpx_geteb3d_energy_conserving(&np, xp.dataPtr(), yp.dataPtr(), zp.dataPtr(),
					Exp->dataPtr(),Eyp->dataPtr(),Ezp->dataPtr(),
					Bxp->dataPtr(),Byp->dataPtr(),Bzp->dataPtr(),
					&xyzmin[0], &xyzmin[1], &xyzmin[2],
					&dx[0], &dx[1], &dx[2],
					&nx, &ny, &nz, &ng, &ng, &ng, 
					&order, &order, &order, 
					exfab.dataPtr(), eyfab.dataPtr(), ezfab.dataPtr(),
					bxfab.dataPtr(), byfab.dataPtr(), bzfab.dataPtr(),
					&ll4symtry, &l_lower_order_in_v, &field_gathe_algo);

	BL_PROFILE_VAR_STOP(blp_pxr);
    }
}
