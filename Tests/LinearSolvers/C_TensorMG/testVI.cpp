#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>

#include <unistd.h>

#include <AMReX_Utility.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#ifndef NDEBUG
#ifdef BL_USE_ARRAYVIEW
#include <ArrayView.H>
#endif
#endif

#include <TestMCViscBndry.H>
#include <AMReX_DivVis.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MCMultiGrid.H>
#include <AMReX_MCCGSolver.H>
#include <AMReX_ParallelDescriptor.H>

#include <main_F.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    std::cout << std::setprecision(10);

    if (argc < 2)
    {
      std::cerr << "usage:  " << argv[0] << " inputsfile [options]" << '\n';
      exit(-1);
    }

    ParmParse pp;
    
    int n;

    BoxArray bs;
    
#if BL_SPACEDIM == 2
    Box domain(IntVect(0,0),IntVect(11,11));
    std::string boxfile("gr.2_small_a") ;
#elif BL_SPACEDIM == 3
    Box domain(IntVect(0,0,0),IntVect(11,11,11));
    std::string boxfile("grids/gr.3_2x3x4") ;
#endif
    pp.query("boxes", boxfile);

    std::ifstream ifs(boxfile.c_str(), std::ios::in);

    if (!ifs)
    {
        std::string msg = "problem opening grids file: ";
        msg += boxfile.c_str();
        amrex::Abort(msg.c_str());
    }

    ifs >> domain;

    if (ParallelDescriptor::IOProcessor())
	std::cout << "domain: " << domain << std::endl;

    bs.readFrom(ifs);

    if (ParallelDescriptor::IOProcessor())
	std::cout << "grids:\n" << bs << std::endl;

    Geometry geom(domain);
    const Real* H = geom.CellSize();
    int ratio=2; pp.query("ratio", ratio);

    DistributionMapping dm {bs};
    
    // allocate/init soln and rhs
    int Ncomp=BL_SPACEDIM;
    int Nghost=0;
    int Ngrids=bs.size();
    MultiFab soln(bs, dm, Ncomp, Nghost); soln.setVal(0.0);
    MultiFab out (bs, dm, Ncomp, Nghost); 
    MultiFab rhs (bs, dm, Ncomp, Nghost); rhs.setVal(0.0);
    for(MFIter rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi)
    {
	FORT_FILLRHS(rhs[rhsmfi].dataPtr(),
		     ARLIM(rhs[rhsmfi].loVect()),ARLIM(rhs[rhsmfi].hiVect()),
		     H,&Ncomp);
    }
    
    // Create the boundary object
    MCViscBndry vbd(bs,dm,geom);

    BCRec phys_bc;
    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    
    // Create the BCRec's interpreted by ViscBndry objects
#if BL_SPACEDIM==2
    Vector<BCRec> pbcarray(4);
    pbcarray[0] = BCRec(AMREX_D_DECL(REFLECT_ODD,REFLECT_EVEN,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[1] = BCRec(AMREX_D_DECL(REFLECT_EVEN,REFLECT_ODD,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[2] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[3] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
#elif BL_SPACEDIM==3
    Vector<BCRec> pbcarray(12);

#if 1
    pbcarray[0] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
    pbcarray[1] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
    pbcarray[2] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
    pbcarray[3] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[4] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[5] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[6] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[7] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[8] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[9] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[10] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			 AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[11] = BCRec(AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			 AMREX_D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
#else
    for (int i = 0; i < 12; i++)
        pbcarray[i] = phys_bc;
#endif
#endif
    
    Nghost = 1; // need space for bc info
    MultiFab fine(bs,dm,Ncomp,Nghost);
    for(MFIter finemfi(fine); finemfi.isValid(); ++finemfi)
    {
	FORT_FILLFINE(fine[finemfi].dataPtr(),
		      ARLIM(fine[finemfi].loVect()),ARLIM(fine[finemfi].hiVect()),
		      H,&Ncomp);
    }

    // Create "background coarse data"
    Box crse_bx = Box(domain).coarsen(ratio).grow(1);
    BoxArray cba(crse_bx);
    cba.maxSize(32);
    Real h_crse[BL_SPACEDIM];
    for (n=0; n<BL_SPACEDIM; n++) h_crse[n] = H[n]*ratio;

    DistributionMapping cdm{cba};

    MultiFab crse_mf(cba, cdm, Ncomp, 0);
//    FArrayBox crse_fab(crse_bx,Ncomp);

    for (MFIter mfi(crse_mf); mfi.isValid(); ++mfi)
    {
        FORT_FILLCRSE(crse_mf[mfi].dataPtr(),
                      ARLIM(crse_mf[mfi].loVect()),ARLIM(crse_mf[mfi].hiVect()),
                      h_crse,&Ncomp);
    }


    
    // Create coarse boundary register, fill w/data from coarse FAB
    int bndry_InRad=0;
    int bndry_OutRad=1;
    int bndry_Extent=1;
    BoxArray cbs = BoxArray(bs).coarsen(ratio);
    BndryRegister cbr(cbs,dm,bndry_InRad,bndry_OutRad,bndry_Extent,Ncomp);
    for (OrientationIter face; face; ++face)
    {
	Orientation f = face();
	FabSet& bnd_fs(cbr[f]);
	bnd_fs.copyFrom(crse_mf, 0, 0, 0, Ncomp);
    }
  
    // Interpolate crse data to fine boundary, where applicable
    int cbr_Nstart=0;
    int fine_Nstart=0;
    int bndry_Nstart=0;
    vbd.setBndryValues(cbr,cbr_Nstart,fine,fine_Nstart,
		       bndry_Nstart,Ncomp,ratio,pbcarray);
  
    Nghost = 1; // other variables don't need extra space
    
    DivVis lp(vbd,H);
    
    Real a = 0.0;
    Real b[BL_SPACEDIM];
    b[0] = 1.0;
    b[1] = 1.0;
#if BL_SPACEDIM>2
    b[2] = 1.0;
#endif
    MultiFab  acoefs;
    int NcompA = (BL_SPACEDIM == 2  ?  2  :  1);
    acoefs.define(bs, dm, NcompA, Nghost);
    acoefs.setVal(a);
    MultiFab bcoefs[BL_SPACEDIM];
    for (n=0; n<BL_SPACEDIM; ++n)
    {
	BoxArray bsC(bs);
	bcoefs[n].define(bsC.surroundingNodes(n), dm, 1, Nghost);
#if 1
	for(MFIter bmfi(bcoefs[n]); bmfi.isValid(); ++bmfi)
	{
	    FORT_MAKEMU(bcoefs[n][bmfi].dataPtr(),
			ARLIM(bcoefs[n][bmfi].loVect()),ARLIM(bcoefs[n][bmfi].hiVect()),H,n);
	}
#else
	bcoefs[n].setVal(b[n]);
#endif
    } // -->> over dimension
    lp.setCoefficients(acoefs, bcoefs);
#if 1
    lp.maxOrder(4);
#endif
    
    Nghost = 1;
    MultiFab tsoln(bs, dm, Ncomp, Nghost); 
    tsoln.setVal(0.0);
#if 1
    tsoln.copy(fine);
#endif
#if 0
    // testing apply
    lp.apply(out,tsoln);
    Box subbox = out[0].box();
    Real n1 = out[0].norm(subbox,1,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
    ParallelDescriptor::ReduceRealSum(n1);
    if (ParallelDescriptor::IOProcessor())
    {
	cout << "n1 output is "<<n1<<std::endl;
    }
    out.minus(rhs,0,BL_SPACEDIM,0);
    // special to single grid prob
    Real n2 = out[0].norm(subbox,1,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
    ParallelDescriptor::ReduceRealSum(n2);
    if (ParallelDescriptor::IOProcessor())
    {
	cout << "n2 difference is "<<n2<<std::endl;
    }
#if 0
    subbox.grow(-1);
    Real n3 = out[0].norm(subbox,0,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
    ParallelDescriptor::ReduceRealMax(n3);
    if (ParallelDescriptor::IOProcessor())
    {
	cout << "n3 difference is "<<n3<<std::endl;
    }
#endif
    
#endif
    
    const IntVect refRatio(AMREX_D_DECL(2,2,2));
    const Real bgVal = 1.0;
    
#if 1
#ifndef NDEBUG
    // testing flux computation
    BoxArray xfluxbox(bs);
    xfluxbox.surroundingNodes(0);
    MultiFab xflux(xfluxbox,dm,Ncomp,Nghost);
    xflux.setVal(1.e30);
    BoxArray yfluxbox(bs);
    yfluxbox.surroundingNodes(1);
    MultiFab yflux(yfluxbox,dm,Ncomp,Nghost);
    yflux.setVal(1.e30);
#if BL_SPACEDIM>2
    BoxArray zfluxbox(bs);
    zfluxbox.surroundingNodes(2);
    MultiFab zflux(zfluxbox,dm,Ncomp,Nghost);
    zflux.setVal(1.e30);
#endif
    lp.compFlux(xflux,
		yflux,
#if BL_SPACEDIM>2
		zflux,
#endif
		tsoln);
    
    // Write fluxes
    //writeMF(&xflux,"xflux.mfab");
    //writeMF(&yflux,"yflux.mfab");
#if BL_SPACEDIM>2
    //writeMF(&zflux,"zflux.mfab");
#endif
    
#endif
#endif
    
    Real tolerance = 1.0e-10; pp.query("tol", tolerance);
    Real tolerance_abs = 1.0e-10; pp.query("tol_abs", tolerance_abs);

#if 0
    cout << "Bndry Data object:" << std::endl;
    cout << lp.bndryData() << std::endl;
#endif
    
#if 0
    bool use_mg_pre = false;
    MCCGSolver cg(lp,use_mg_pre);
    cg.solve(soln,rhs,tolerance,tolerance_abs);
#else
    MCMultiGrid mg(lp);
    mg.solve(soln,rhs,tolerance,tolerance_abs);
#endif

#if 0
    cout << "MCLinOp object:" << std::endl;
    cout << lp << std::endl;
#endif
    
    VisMF::Write(soln,"soln");
    
#if 0
    // apply operator to soln to see if really satisfies eqn
    tsoln.copy(soln);
    lp.apply(out,tsoln);
    soln.copy(out);
    // Output "apply" results on soln
    VisMF::Write(soln,"apply");

    // Compute truncation
    for (MFIter smfi(soln); smfi.isValid(); ++smfi)
    {
	soln[smfi] -= fine[smfi];
    }
    for( int icomp=0; icomp < BL_SPACEDIM ; icomp++ )
    {
	Real solnMin = soln.min(icomp);
	Real solnMax = soln.max(icomp);
	ParallelDescriptor::ReduceRealMin(solnMin);
	ParallelDescriptor::ReduceRealMax(solnMax);
	if (ParallelDescriptor::IOProcessor())
	{
	    cout << icomp << "  "<<solnMin << " " << solnMax <<std::endl;
	}
    }
    // Output truncation
    VisMF::Write(soln,"trunc");
#endif

    int dumpLp=0; pp.query("dumpLp",dumpLp);
    bool write_lp = (dumpLp == 1 ? true : false);
    if (write_lp)
	std::cout << lp << std::endl;

    // Output trunc
    ParallelDescriptor::EndParallel();
}
