
#include <ParmParse.H>

#include <WarpX.H>

WarpX::WarpX ()
{
    ReadParameters();
    BL_ASSERT(max_level == 0);

    // Geometry
    {
	Geometry::Setup();

	geom_arr.resize(max_level+1);

	Box index_domain(IntVect::TheZeroVector(), IntVect(D_DECL(nx-1,ny-1,nz-1)));

	for (auto& g : geom_arr) {
	    g.define(index_domain);
	    index_domain.refine(ref_ratio);
	}
    }

    // BoxArray
    {
	ba_arr.resize(max_level+1);
	Box index_domain(IntVect::TheZeroVector(), IntVect(D_DECL(nx-1,ny-1,nz-1)));
	// create level 0 BoxArray only
	ba_arr[0].define(index_domain);
	ba_arr[0].maxSize(max_grid_size);
    }

    // DistributionMapping
    {
	dmap_arr.resize(max_level+1);
	// create level 0 DistributionMapping only
	MultiFab foo(ba_arr[0],1,0);
	dmap_arr[0] = foo.DistributionMap();
    }

    // Particle Container
    if (!mypc)
    {
	Array<int> rr(max_level);
	mypc = std::unique_ptr<MyParticleContainer>
	    (new MyParticleContainer(geom_arr, dmap_arr, ba_arr, rr));
	mypc->SetVerbose(0);
    }

    // Assuming there is only a single level

    // PICSAR assumes all fields are nodal plus one ghost cell.
    IntVect nodalflag = IntVect::TheUnitVector();
    int ng = 1;

    for (int i = 0; i < BL_SPACEDIM; ++i) {
	current.push_back(std::unique_ptr<MultiFab>
			  (new MultiFab(ba_arr[0],1,ng,dmap_arr[0],Fab_allocate,nodalflag)));

	Efield.push_back(std::unique_ptr<MultiFab>
			 (new MultiFab(ba_arr[0],1,ng,dmap_arr[0],Fab_allocate,nodalflag)));

	Bfield.push_back(std::unique_ptr<MultiFab>
			 (new MultiFab(ba_arr[0],1,ng,dmap_arr[0],Fab_allocate,nodalflag)));	
    }
}

WarpX::~WarpX ()
{
}

void
WarpX::ReadParameters ()
{
    {
	ParmParse pp;

	pp.query("verbose", verbose);

	pp.query("cfl", cfl);

	pp.get("nx", nx);
	pp.get("ny", ny);
	pp.get("nz", nz);
	pp.get("max_step", max_step);
	
	pp.query("max_grid_size", max_grid_size);
	pp.query("max_level", max_level);
	pp.query("ref_ratio", ref_ratio);

	pp.query("plot_int", plot_int);
    }
}

void 
WarpX::FillBoundary (MultiFab& mf, const Geometry& geom, const IntVect& nodalflag)
{
    const IndexType correct_typ(nodalflag);
    BoxArray ba = mf.boxArray();
    ba.convert(correct_typ);

    MultiFab tmpmf(ba, mf.nComp(), mf.nGrow(), mf.DistributionMap());

    const IndexType& mf_typ = mf.boxArray().ixType();

    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi) {
	FArrayBox& tmpfab = tmpmf[mfi];
	tmpfab.SetBoxType(mf_typ);
	tmpfab.copy(mf[mfi]);
	tmpfab.SetBoxType(correct_typ);
    }

    tmpmf.FillBoundary(geom.periodicity());

    for (MFIter mfi(tmpmf); mfi.isValid(); ++mfi) {
	FArrayBox& tmpfab = tmpmf[mfi];
	tmpfab.SetBoxType(mf_typ);
	mf[mfi].copy(tmpmf[mfi]);
	tmpfab.SetBoxType(correct_typ);
    }
}
