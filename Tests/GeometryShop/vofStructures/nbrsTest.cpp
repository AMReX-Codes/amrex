#include <cmath>
#include <type_traits>

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_FlatPlateGeom.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBISBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_BLProfiler.H"
#include "AMReX_BLFort.H"
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_MeshRefine.H"

#include "AMReX_ArrayLim.H"
#include "AMReX_EBRedist.H"

using namespace amrex;

static void
get_EBLG(EBLevelGrid& eblg)
{
    GridParameters params;
    getGridParameters(params, 0, true);

    ParmParse pp;

    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int which_geom = params.whichGeom;
    bool verbosity=true;
    pp.query("verbosity",verbosity);
    if (which_geom == 0)
    {
        RealVect origin = RealVect::Zero;
        Real radius = 0.5;
        pp.get("sphere_radius", radius);
        std::vector<Real> centervec(SpaceDim);
        pp.getarr("sphere_center", centervec, 0, SpaceDim);
        RealVect center;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            center[idir] = centervec[idir];
        }

        //bool insideRegular = false;
        bool insideRegular = true;
        SphereIF sphere(radius, center, insideRegular);

        GeometryShop gshop(sphere, verbosity);
        ebisPtr->define(params.coarsestDomain, origin, params.coarsestDx, gshop);
    }
    else if (which_geom == 1) 
    {
        RealVect origin = RealVect::Zero;
        std::vector<Real>  platelovec(SpaceDim);
        std::vector<Real>  platehivec(SpaceDim);

        int normalDir;
        Real plateLoc;
    
        pp.getarr("plate_lo", platelovec, 0, SpaceDim);
        pp.getarr("plate_hi", platehivec, 0, SpaceDim);
        pp.get("plate_location", plateLoc);
        pp.get("plate_normal", normalDir);

        RealVect plateLo, plateHi;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            plateLo[idir] = platelovec[idir];
            plateHi[idir] = platehivec[idir];
        }
        FlatPlateGeom flat_plate(normalDir, plateLoc, plateLo, plateHi);
        ebisPtr->define(params.coarsestDomain, origin, params.coarsestDx, flat_plate);
    }
    else
    {
        amrex::Abort("Unknown geom_type");
    }

    std::vector<EBLevelGrid> eblg_tmp;
    getAllIrregEBLG(eblg_tmp, params);
    eblg = eblg_tmp[0];

}

static void
Copy(intDIM& out, const IntVect& in)
{
    
    for (int d=0; d<d; ++d) {
	out[d] = 0;
    }
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

/*
static void
Copy(IntVect& out, const intDIM& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}
*/
extern "C"
{
    void fill_redist_stencil(const int*  lo, const int*  hi,
                             const int* slo, const int* shi,
			     NbrSten* redistSten, const int* Nsten,
			     const BL_FORT_IFAB_ARG_3D(imask),
			     BL_FORT_FAB_ARG_3D(vfrac));
    
    void apply_redist_stencil(const int* lo,  const int*  hi,
                              const int* slo, const int* shi,
                              NbrSten* redistSten, const int* Nsten,
			      const BL_FORT_FAB_ARG_3D(vin),
			      BL_FORT_FAB_ARG_3D(vout));
}

int myTest()
{
    EBLevelGrid eblg;
    get_EBLG(eblg);
    const BoxArray& ba = eblg.getBoxArray();
    const DistributionMapping& dm = eblg.getDM();
    const Box& domain = eblg.getDomain();
    Geometry geom(domain); 

    // Create a mask with the following values 
    // 1=reg, 0=sv, -1=covered, 2=mv, -2=outside domain, -3=crse
    int covered_val = 10; // Will be overwritten with FillBoundary
    int notcovered_val = -3;
    int physbnd_val = -2;
    int interior_val = 10; // Will be overwritten
    FabArray<BaseFab<int> > imask(ba, dm, 1, 1);  // 1 grow for 3-wide stencil
    imask.BuildMask(geom.Domain(),geom.periodicity(),covered_val,notcovered_val,physbnd_val,interior_val);
    
    BaseFab<int> tm;
    for(MFIter  mfi(imask); mfi.isValid(); ++mfi)
    {
	BaseFab<int>& m = imask[mfi];
        const EBISBox& ebis = eblg.getEBISL()[mfi];
	ebis.getEBGraph().fillIntMask(tm);
	m.copy(tm);
	IntVectSet mv = ebis.getMultiCells(mfi.validbox());
	for (IVSIterator it(mv); it.ok(); ++it)
	{
	    m(it(),0) = 2;
	}
    }
    imask.FillBoundary();

    MultiFab vfrac(ba,dm,1,1);
    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis = eblg.getEBISL()[mfi];
        const Box& vbox = mfi.validbox();
	FArrayBox& vfab = vfrac[mfi];
	vfab.setVal(0);
	
        if (ebis.isAllRegular())
	{
	    vfab.setVal(1);
	}
	else
        {
	    for (IntVect iv(vbox.smallEnd()); iv<=vbox.bigEnd(); vbox.next(iv))
	    {
		std::vector<VolIndex> vofs = ebis.getVoFs(iv);
		Real vtot = 0;
		for (int j=0; j<vofs.size(); ++j)
		{
		    vtot += ebis.volFrac(vofs[j]);
		}
		vfab(iv,0) = vtot;
	    }
	}
    }
    vfrac.FillBoundary();

    const Box stenBox(IntVect(D_DECL(-1,-1,-1)),
                      IntVect(D_DECL(+1,+1,+1)));

    std::map<int,std::vector<NbrSten>> redistSten;
    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis = eblg.getEBISL()[mfi];
        const Box& vbox = mfi.validbox();
	int i = mfi.index();
	BaseFab<int>& m = imask[mfi];
	FArrayBox& vfab = vfrac[mfi];
	
        if ( !(ebis.isAllRegular())) 
        {
            IntVectSet isIrreg = ebis.getIrregIVS(vbox);
            for (IVSIterator it(isIrreg); it.ok(); ++it)
            {
                const IntVect& iv=it();
                std::vector<VolIndex> vofs = ebis.getVoFs(iv);
                for (int ivof=0; ivof<vofs.size(); ++ivof)
                {
		    redistSten[i].push_back(NbrSten());
		    NbrSten& nbr = redistSten[i].back();
		    Copy(nbr.iv,iv);
                }
            }

	    int Nsten = redistSten[i].size();
	    fill_redist_stencil(BL_TO_FORTRAN_BOX(vbox),
                                BL_TO_FORTRAN_BOX(stenBox),
                                redistSten[i].data(), &Nsten,
	        		BL_TO_FORTRAN_3D(m),
	        		BL_TO_FORTRAN_3D(vfab));
        }
    }

    MultiFab tin(ba,dm,1,1);
    tin.setVal(2);
    MultiFab tout(ba,dm,1,1);
    tin.setVal(0);
    for(MFIter  mfi(tin); mfi.isValid(); ++mfi)
    {
        const Box& vbox = mfi.validbox();
	int i = mfi.index();
	const FArrayBox& infab = tin[mfi];
	FArrayBox& outfab = tout[mfi];
	int Nsten = redistSten[i].size();
	
        if (Nsten > 0)
        {
	    apply_redist_stencil(BL_TO_FORTRAN_BOX(vbox),
                                 BL_TO_FORTRAN_BOX(stenBox),
                                 redistSten[i].data(), &Nsten,
                                 BL_TO_FORTRAN_3D(infab),
                                 BL_TO_FORTRAN_3D(outfab));
        }
    }

    return 0;
}

int
main(int argc,char **argv)
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE_VAR("main()", pmain);

        int retFlag = myTest();

        if (retFlag != 0)
        {
            amrex::Print() << "non zero return detected = " << retFlag << '\n';
            amrex::Print() << "sphere test failed" << '\n';
        }
        else
        {
            amrex::Print() << "test passed" << '\n';
        }

        BL_PROFILE_VAR_STOP(pmain);
    }
    amrex::Finalize();
    return 0;
}
