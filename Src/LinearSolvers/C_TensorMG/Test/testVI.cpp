#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>

#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include <Utility.H>
#include <Tracer.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <VisMF.H>
#ifndef NDEBUG
#ifdef BL_USE_ARRAYVIEW
#include <ArrayView.H>
#endif
#include <TV_TempWrite.H>
#endif

#include <TestMCViscBndry.H>
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>
#include <ParallelDescriptor.H>

#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif

#include <main_F.H>
#include <WritePlotFile.H>

BoxList
readBoxList(aString file,
	    BOX&    domain );

int
main (int   argc,
      char* argv[])
{
    //
    // Make sure to catch new failures.
    //
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif

    ParallelDescriptor::StartParallel(&argc,&argv);

    cout << setprecision(10);

    if(argc < 2) {
      cerr << "usage:  " << argv[0] << " inputsfile [options]" << '\n';
      exit(-1);
    }

    ParmParse pp(argc-2,argv+2,NULL,argv[1]); 
    //
    // Initialize random seed after we're running in parallel.
    //
    Utility::InitRandom(ParallelDescriptor::MyProc() + 1);
    
#ifndef WIN32
    int sleeptime = 0; pp.query("sleep", sleeptime);
    sleep(sleeptime);
#endif
    
    TRACER("mcmg");
    
    int n;
    
#if BL_SPACEDIM == 2
    Box container(IntVect(0,0),IntVect(11,11));
    aString boxfile("gr.2_small_a") ;
#elif BL_SPACEDIM == 3
    Box container(IntVect(0,0,0),IntVect(11,11,11));
    aString boxfile("grids/gr.3_small_a") ;
#endif
    pp.query("boxes", boxfile);
    BoxArray bs(readBoxList(boxfile,container));
    Geometry geom(container);
    const Real* H = geom.CellSize();
    int ratio=2; pp.query("ratio", ratio);

    // allocate/init soln and rhs
    int Ncomp=BL_SPACEDIM;
    int Nghost=0;
    int Ngrids=bs.length();
    MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
    MultiFab out(bs, Ncomp, Nghost, Fab_allocate); 
    MultiFab rhs(bs, Ncomp, Nghost, Fab_allocate); rhs.setVal(0.0);
    for(MultiFabIterator rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi)
    {
	FORT_FILLRHS(rhsmfi().dataPtr(),
		     ARLIM(rhsmfi().loVect()),ARLIM(rhsmfi().hiVect()),
		     H,&Ncomp);
    }
    
    // Create the boundary object
    MCViscBndry vbd(bs,geom);
    
    // Create the BCRec's interpreted by ViscBndry objects
#if BL_SPACEDIM==2
    Array<BCRec> pbcarray(4);
    pbcarray[0] = BCRec(D_DECL(REFLECT_ODD,REFLECT_EVEN,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[1] = BCRec(D_DECL(REFLECT_EVEN,REFLECT_ODD,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[2] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[3] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
#elif BL_SPACEDIM==3
    Array<BCRec> pbcarray(12);
#define TEST_EXT 0
#if TEST_EXT
    pbcarray[0] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
    pbcarray[1] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
    pbcarray[2] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
#else
    pbcarray[0] = BCRec(INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR);
    pbcarray[1] = BCRec(INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR);
    pbcarray[2] = BCRec(INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR);
#endif
    // the other bc's don't really matter.  Just init with anything
    pbcarray[3] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[4] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[5] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[6] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[7] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[8] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[9] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[10] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			 D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
    pbcarray[11] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			 D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
#endif
    
    Nghost = 1; // need space for bc info
    MultiFab fine(bs,Ncomp,Nghost,Fab_allocate);
    for(MultiFabIterator finemfi(fine); finemfi.isValid(); ++finemfi)
    {
	FORT_FILLFINE(finemfi().dataPtr(),
		      ARLIM(finemfi().loVect()),ARLIM(finemfi().hiVect()),
		      H,&Ncomp);
    }

    // Create "background coarse data"
    BOX crse_bx = Box(container).coarsen(ratio).grow(1);
    REAL h_crse[BL_SPACEDIM];
    for (n=0; n<BL_SPACEDIM; n++) h_crse[n] = H[n]*ratio;
    FARRAYBOX crse_fab(crse_bx,Ncomp);
    FORT_FILLCRSE(crse_fab.dataPtr(),
		  ARLIM(crse_fab.loVect()),ARLIM(crse_fab.hiVect()),
		  h_crse,&Ncomp);
    
    // Create coarse boundary register, fill w/data from coarse FAB
    int bndry_InRad=0;
    int bndry_OutRad=1;
    int bndry_Extent=1;
    BoxArray cbs = BoxArray(bs).coarsen(ratio);
    BndryRegister cbr(cbs,bndry_InRad,bndry_OutRad,bndry_Extent,Ncomp);
    for (OrientationIter face; face; ++face)
    {
	Orientation f = face();
	FabSet& bnd_fs(cbr[f]);
	bnd_fs.copyFrom(crse_fab);
    }
  
    // Interpolate crse data to fine boundary, where applicable
    int cbr_Nstart=0;
    int fine_Nstart=0;
    int bndry_Nstart=0;
    vbd.setBndryValues(cbr,cbr_Nstart,fine,fine_Nstart,
		       bndry_Nstart,Ncomp,ratio,pbcarray);
  
    Nghost = 1; // other variables don't need extra space
    
    DivVis lp(vbd,H);
    
    REAL a = 0.0;
    REAL b[BL_SPACEDIM];
    b[0] = 1.0;
    b[1] = 1.0;
#if BL_SPACEDIM>2
    b[2] = 1.0;
#endif
    MultiFab  acoefs;
    int NcompA = (BL_SPACEDIM == 2  ?  2  :  1);
    acoefs.define(bs, NcompA, Nghost, Fab_allocate);
    acoefs.setVal(a);
    MultiFab bcoefs[BL_SPACEDIM];
    for (n=0; n<BL_SPACEDIM; ++n)
    {
	BoxArray bsC(bs);
	bcoefs[n].define(bsC.surroundingNodes(n), 1,
			 Nghost, Fab_allocate);
#if 1
	for(MultiFabIterator bmfi(bcoefs[n]); bmfi.isValid(); ++bmfi)
	{
	    FORT_MAKEMU(bmfi().dataPtr(),
			ARLIM(bmfi().loVect()),ARLIM(bmfi().hiVect()),H,n);
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
    MultiFab tsoln(bs, Ncomp, Nghost, Fab_allocate); 
    tsoln.setVal(0.0);
#if 1
    tsoln.copy(fine);
#endif
#if 0
    // testing apply
    lp.apply(out,tsoln);
    Box subbox = out[0].box();
    REAL n1 = out[0].norm(subbox,1,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
    ParallelDescriptor::ReduceRealSum(n1);
    if (ParallelDescriptor::IOProcessor())
    {
	cout << "n1 output is "<<n1<<endl;
    }
    out.minus(rhs,0,BL_SPACEDIM,0);
    // special to single grid prob
    REAL n2 = out[0].norm(subbox,1,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
    ParallelDescriptor::ReduceRealSum(n2);
    if (ParallelDescriptor::IOProcessor())
    {
	cout << "n2 difference is "<<n2<<endl;
    }
#if 0
    subbox.grow(-1);
    REAL n3 = out[0].norm(subbox,0,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
    ParallelDescriptor::ReduceRealMax(n3);
    if (ParallelDescriptor::IOProcessor())
    {
	cout << "n3 difference is "<<n3<<endl;
    }
#endif
    
#endif
    
    const IntVect refRatio(D_DECL(2,2,2));
    const Real bgVal = 1.0;
    
#if 1
#ifndef NDEBUG
    // testing flux computation
    BoxArray xfluxbox(bs);
    xfluxbox.surroundingNodes(0);
    MultiFab xflux(xfluxbox,Ncomp,Nghost,Fab_allocate);
    xflux.setVal(1.e30);
    BoxArray yfluxbox(bs);
    yfluxbox.surroundingNodes(1);
    MultiFab yflux(yfluxbox,Ncomp,Nghost,Fab_allocate);
    yflux.setVal(1.e30);
#if BL_SPACEDIM>2
    BoxArray zfluxbox(bs);
    zfluxbox.surroundingNodes(2);
    MultiFab zflux(zfluxbox,Ncomp,Nghost,Fab_allocate);
    zflux.setVal(1.e30);
#endif
    lp.compFlux(xflux,
		yflux,
#if BL_SPACEDIM>2
		zflux,
#endif
		tsoln);
    
    // Write fluxes
    writeMF(&xflux,"xflux.mfab");
    writeMF(&yflux,"yflux.mfab");
#if BL_SPACEDIM>2
    writeMF(&zflux,"zflux.mfab");
#endif
    
#endif
#endif
    
    REAL tolerance = 1.0e-10; pp.query("tol", tolerance);
    REAL tolerance_abs = 1.0e-10; pp.query("tol_abs", tolerance_abs);

#if 0
    cout << "Bndry Data object:" << endl;
    cout << lp.bndryData() << endl;
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
    cout << "MCLinOp object:" << endl;
    cout << lp << endl;
#endif
    
    // Write solution
    writePlotFile("pltfile",soln,geom,refRatio,bgVal);
    
#if 0
    // apply operator to soln to see if really satisfies eqn
    tsoln.copy(soln);
    lp.apply(out,tsoln);
    soln.copy(out);
    // Output "apply" results on soln
    writePlotFile("plt_apply",soln,geom,refRatio,bgVal);

    // Compute truncation
    for(MultiFabIterator smfi(soln); smfi.isValid(); ++smfi)
    {
	DependentMultiFabIterator fmfi(smfi,fine);
	smfi() -= fmfi();
    }
    for( int icomp=0; icomp < BL_SPACEDIM ; icomp++ )
    {
	Real solnMin = soln.min(icomp);
	Real solnMax = soln.max(icomp);
	ParallelDescriptor::ReduceRealMin(solnMin);
	ParallelDescriptor::ReduceRealMax(solnMax);
	if (ParallelDescriptor::IOProcessor())
	{
	    cout << icomp << "  "<<solnMin << " " << solnMax <<endl;
	}
    }
    // Output truncation
    writePlotFile("plt_trunc",soln,geom,refRatio,bgVal);
#endif

    int dumpLp=0; pp.query("dumpLp",dumpLp);
    bool write_lp = (dumpLp == 1 ? true : false);
    if (write_lp)
	cout << lp << endl;

    // Output trunc
    ParallelDescriptor::EndParallel();
}

BoxList
readBoxList(const aString file, BOX& domain )
{
    BoxList retval;
    ifstream boxspec(file.c_str());
    if( !boxspec )
    {
	BoxLib::Error("readBoxList: unable to open " + *file.c_str());
    }
    boxspec >> domain;
    
    int numbox;
    boxspec >> numbox;

    for(int i=0; i<numbox; i++)
    {
	BOX tmpbox;
	boxspec >> tmpbox;
	if( ! domain.contains(tmpbox))
	{
	    cerr << "readBoxList: bogus box " << tmpbox << endl;
	    exit(1);
        }
	retval.append(tmpbox);
    }
    boxspec.close();
    return retval;
}

