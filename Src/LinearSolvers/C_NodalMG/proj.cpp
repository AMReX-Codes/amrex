//BL_COPYRIGHT_NOTICE

#include "hg_projector.H"

#include <CArena.H>
#include <Utility.H>
#include <ParmParse.H>
#include <RunStats.H>

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <iomanip>
#include <fstream>
#include <new>
#ifndef __GNUC__
#include <sstream>
#else
#include <cstdio>
#endif
using namespace std;

#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdio.h>
#include <new.h>
#endif

#ifndef WIN32
#include <unistd.h>
#endif

#ifdef HG_DEBUG
std::ofstream debug_out;
#endif

// bool HG_is_debugging = false;

void projtest(const Array<BoxArray>& m, Array<IntVect>& ratio, Array<Box>& domain);

void driver(const char* filename);

int pcode = 4;
int nrep  = 1;
Real tol = 2.e-10;
int coordsys = 0;
RegType bcall = refWall;
RegType bc00;
RegType bc10;
RegType bc20;
RegType bc01;
RegType bc11;
RegType bc21;

holy_grail_amr_multigrid::stencil hg_stencil = holy_grail_amr_multigrid::cross;

int
main(int argc, char **argv)
{
#ifndef WIN32
    set_new_handler(Utility::OutOfMemory);
#endif
    ParallelDescriptor::StartParallel(&argc, &argv);
    if ( argc < 2 ) BoxLib::Error("expected inputs file");
    ParmParse pp(argc-2,argv+2, 0, argv[1]);

    HG_is_debugging = true;   
    HG_is_debugging = false;   
#ifdef HG_DEBUG
#ifdef __GNUC__
    char buf[1024];
    sprintf(buf, "guf%d_%d", ParallelDescriptor::NProcs(), ParallelDescriptor::MyProc());
    debug_out.open(buf);
#else
    std::ostringstream fname;
    fname << "gu" << ParallelDescriptor::NProcs()
	  << "_" << ParallelDescriptor::MyProc() << std::ends;
    debug_out.open(fname.str().c_str(), ios::trunc);
#endif
    if ( debug_out.fail() ) BoxLib::Error( "Failed to open debug file" );
    debug_out << std::setprecision(15);
#endif
    {
	int i = ParallelDescriptor::MyProc();
	int j = ParallelDescriptor::NProcs();
	FORT_HGDEBUGINIT(&i, &j);
    }
    HG::MPI_init();
#ifndef WIN32
    int slp = 0;
    pp.query("sleep", slp);
    if ( slp > 0 )
    {
	cout << "going to sleep for " << slp
	     << " seconds, debug pid = " << ::getpid() << endl;
	::sleep(slp);
    }
#endif
    aString ass;
    if ( pp.query("bcall", ass) ) bcall = RegTypeSet(ass);
    bc00 = bc10 = bc20 = bc01 = bc11 = bc21 = bcall;
    if ( pp.query("bc00", ass) ) bc00 = RegTypeSet(ass);
    if ( pp.query("bc10", ass) ) bc10 = RegTypeSet(ass);
    if ( pp.query("bc20", ass) ) bc20 = RegTypeSet(ass);
    if ( pp.query("bc01", ass) ) bc01 = RegTypeSet(ass);
    if ( pp.query("bc11", ass) ) bc11 = RegTypeSet(ass);
    if ( pp.query("bc21", ass) ) bc21 = RegTypeSet(ass);
    pp.query("nrep", nrep);
    pp.query("pcode", pcode);
    pp.query("tol", tol);
    pp.query("coordsys", coordsys);
    aString stencil;
    pp.query("stencil", stencil);
    if ( stencil == "cross" )
    {
	hg_stencil = holy_grail_amr_multigrid::cross;
    }
    else if ( stencil == "terrain" )
    {
	hg_stencil = holy_grail_amr_multigrid::terrain;
    }
    else if ( stencil == "full" )
    {
	hg_stencil = holy_grail_amr_multigrid::full;
    }
    else
    {
      BoxLib::Error("stencil must be cross, terrain, or full");
    }

    if ( ParallelDescriptor::IOProcessor() )
    {
	cout << "nrep = " << nrep << endl;
	cout << "pcode = " << pcode << endl;
	cout << "tol = " << tol << endl;
    }
    cout << setprecision(15);

    int num = pp.countname("file");
    for ( int k = 0; k < num; k++)
    {
	aString filename;
	pp.getkth("file", k, filename, 0);
	if ( ParallelDescriptor::IOProcessor() )
	{
	    cout << "file " << k << " is " << filename << endl;
	}
	driver(filename.c_str());
    }

#ifdef HG_DEBUG
    debug_out.close();
#endif
    RunStats::report(cout);

    if (CArena* arena = dynamic_cast<CArena*>(The_FAB_Arena))
    {
        //
        // A barrier to make sure our output follows that of RunStats.
        //
        ParallelDescriptor::Barrier();
        //
        // We're using a CArena -- output some FAB memory stats.
        //
        // This'll output total # of bytes of heap space in the Arena.
        //
        // It's actually the high water mark of heap space required by FABs.
        //
        char buf[256];

        sprintf(buf,
                "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
                ParallelDescriptor::MyProc(),
                arena->heap_space_used());

        cout << buf << endl;
    }

    HG::MPI_finish();

    ParallelDescriptor::EndParallel();

    return 0;
}

void
driver(const char *filename)
{
    Array<BoxArray> m;
    Array<IntVect> ratio;
    Array<Box> domain;
    
    fstream grid;
    grid.open(filename, ios::in);
    if ( grid.fail() )
    {
	BoxLib::Warning("Failed to open grid file");
	return;
    }
    amr_multigrid::mesh_read(m, ratio, domain, grid);
    grid.close();
    projtest(m, ratio, domain);
}

void
init(PArray<MultiFab> u[], PArray<MultiFab>& p,
     const Array<BoxArray>& m, const Array<IntVect>& ratio)
{
#if (BL_SPACEDIM == 2)
    for (int ilev = 0; ilev < m.length(); ilev++) 
    {
	u[0][ilev].setVal(0.0);
	u[1][ilev].setVal(0.0);
    }
    if (m.length() == 1) 
    {
	for (MultiFabIterator u_mfi(u[0][0]); u_mfi.isValid(); ++u_mfi)
	{
	    (*u_mfi)(m[0][u_mfi.index()].smallEnd() + IntVect(2,2)) = 3.0;
	    //u[0][0][igrid](m[0][igrid].smallEnd() + IntVect(3,3)) = 1.0;
	    //u[1][0][igrid](m[0][igrid].smallEnd() + IntVect(3,3)) = 1.0;
	    //u[0][0][igrid](m[0][igrid].smallEnd() + IntVect(20,2)) = 3.0;
	    //u[1][0][igrid](m[0][igrid].smallEnd() + IntVect(2,2)) = 3.0;
	}
    }
    else if (m.length() == 2) 
    {
	for (MultiFabIterator u_mfi(u[0][1]); u_mfi.isValid(); ++u_mfi)
	{
	    (*u_mfi)(m[1][u_mfi.index()].smallEnd() + IntVect(2,2)) = 3.0;
	}
	//u[0][1][0](IntVect(20,90)) = 1.0;
	//u[0][1][0](IntVect(50,50)) = 1.0;
	//u[0][1][0](IntVect(22,12)) = 1.0;
	if (hg_stencil != holy_grail_amr_multigrid::terrain)
	{
	    if ( is_local(u[0][0], 0) )
		u[0][0][0](IntVect(12,12)) = 3.0;
	}
	else
	{
//	    if ( is_local(u[0][0], 0) )
//		u[0][0][0](IntVect(12,12)) = 3.0 * ratio[0][0];
	    if ( is_local(u[0][0], 0) )
		u[0][0][0](IntVect(12,12)) = 3.0;
	    /*
	      if (m[0].domain().length(0) == 32)
	      u[0][1][0](IntVect(30,30)) = 1.0;
	      else if (m[0].domain().length(0) == 64)
	      u[0][1][0](IntVect(60,60)) = 1.0;
	      else
	      u[0][1][0](IntVect(120,120)) = 1.0;
	    */
	}
    }
    else if (m.length() == 3)
    {
	for (MultiFabIterator u_mfi(u[0][1]); u_mfi.isValid(); ++u_mfi)
	{
	    (*u_mfi)(m[1][u_mfi.index()].smallEnd() + IntVect(2,2)) = 3.0;
	}
	//u[0][1][0](IntVect(20,90)) = 1.0;
	//u[0][1][0](IntVect(50,50)) = 1.0;
	//u[0][1][0](IntVect(22,12)) = 1.0;
    }
    else 
    {
	for (int ilev = 0; ilev < m.length(); ilev++) 
	{
	    for ( MultiFabIterator u_mfi(u[0][ilev]); u_mfi.isValid(); ++u_mfi)
	    {
		DependentMultiFabIterator u_dmfi(u_mfi, u[1][ilev]);
		//u[0][ilev][igrid].setVal(1.e20);
		//u[1][ilev][igrid].setVal(3.e20);
		u_mfi->setVal(0.0, m[ilev][u_mfi.index()], 0);
		u_dmfi->setVal(0.0, m[ilev][u_mfi.index()], 0);
	    }
	}
	if ( is_local(u[0][2], 0) )
	{
	    u[0][2][0](m[2][0].smallEnd() + IntVect(2,2)) = 3.0;
	}
	// for gr2ann
	//u[0][2][0](IntVect(20,20)) = 1.0;
	//u[0][2][0](IntVect(20,20)) = 0.0;
	//u[0][2][2](IntVect(70,80)) = 1.0;
    }
    //u[0][1][0](IntVect(31,31)) = -1.0;
    //u[1][1][0](IntVect(31,30)) = 1.0;
    //u[1][1][0](IntVect(30,31)) = -1.0;
    for (int ilev = 0; ilev < p.length(); ilev++)
    {
	p[ilev].setVal(0.0);
    }
#else
    for (int ilev = 0; ilev < m.length(); ilev++) 
    {
	u[0][ilev].setVal(0.0);
	u[1][ilev].setVal(0.0);
	u[2][ilev].setVal(0.0);
    }
    if (m.length() == 1) 
    {
	//int ioff = m[0].domain().length(0) / 8;
	int ioff = 2;
	for (MultiFabIterator u_mfi(u[0][0]); u_mfi.isValid(); ++u_mfi)
	{
	    // Used for timings3_94:
	    //u[0][0][igrid](m[0][igrid].smallEnd() + IntVect(2,2,2)) = 3.0;
	    (*u_mfi)(m[0][u_mfi.index()].smallEnd() + IntVect(ioff,ioff,ioff)) = 3.0;
	}
    }
    else if (m.length() == 2) 
    {
	// used for convergence-rate tests:
	//int ioff = m[1].domain().length(0) / 8;
	int ioff = 2;
	for (MultiFabIterator u_mfi(u[0][1]); u_mfi.isValid(); ++u_mfi)
	{
	    // Used for timings3_94:
	    //u[0][1][igrid](m[1][igrid].smallEnd() + IntVect(2,2,2)) = 3.0;
	    (*u_mfi)(m[1][u_mfi.index()].smallEnd() + IntVect(ioff,ioff,ioff)) = 3.0;
	}
	if ( is_local(u[0][0], 0) ) u[0][0][0](IntVect(1,1,1)) = 3.0;
    }
    else if (m.length() == 3) 
    {
	int ioff = 2;
	for ( MultiFabIterator u_mfi(u[0][2]); u_mfi.isValid(); ++u_mfi)
	{
	    (*u_mfi)(m[2][u_mfi.index()].smallEnd() + IntVect(ioff,ioff,ioff)) = 3.0;
	}
    }
    for (int ilev = 0; ilev < m.length(); ilev++) 
    {
	p[ilev].setVal(0.0);
    }
#endif
}

#if (BL_SPACEDIM == 2)
void
hb93_test1(PArray<MultiFab> u[], const Array<BoxArray>& m, const Array<Box>& d)
{
    for (int ilev = 0 ; ilev < m.length() ; ilev++) 
    {
	double h = 1.0 / d[ilev].length(0);
	double pi = 3.14159265358979323846;
	for ( MultiFabIterator u_mfi(u[0][ilev]); u_mfi.isValid(); ++u_mfi)
	{
	    DependentMultiFabIterator u_dmfi(u_mfi, u[1][ilev]);
	    int igrid = u_mfi.index();
	    for (int i = m[ilev][igrid].smallEnd(0);
		 i <= m[ilev][igrid].bigEnd(0); i++) 
	    {
		for (int j = m[ilev][igrid].smallEnd(1);
		     j <= m[ilev][igrid].bigEnd(1); j++) 
		{
		    double x = (i + 0.5) * h;
		    double y = (j + 0.5) * h;
		    (*u_mfi)(IntVect(i,j)) = -0.5*(1.0-cos(2*pi*x))*sin(2*pi*y);
		    (*u_dmfi)(IntVect(i,j)) =  0.5*(1.0-cos(2*pi*y))*sin(2*pi*x);
		    //u[0][ilev][igrid](IntVect(i,j)) =  0.2*(x+1)*sin(pi*x)*
		    //  (pi*(y+1)*cos(pi*y)+sin(pi*y));
		    //u[1][ilev][igrid](IntVect(i,j)) = -0.2*(y+1)*sin(pi*y)*
		    //  (pi*(x+1)*cos(pi*x)+sin(pi*x));
		}
	    }
	}
    }
}

void
linear_test(PArray<MultiFab> u[], const Array<BoxArray>& m, const Array<Box>& d)
{
    for (int ilev = 0 ; ilev < m.length() ; ilev++) 
    {
	double h = 1.0 / d[ilev].length(0);
	for ( MultiFabIterator u_mfi(u[0][ilev]); u_mfi.isValid(); ++u_mfi)
	{
	    DependentMultiFabIterator u_dmfi(u_mfi, u[1][ilev]);
	    int igrid = u_mfi.index();
	    for (int i = m[ilev][igrid].smallEnd(0);
		 i <= m[ilev][igrid].bigEnd(0); i++) 
	    {
		for (int j = m[ilev][igrid].smallEnd(1);
		     j <= m[ilev][igrid].bigEnd(1); j++) 
		{
		    double x = (i + 0.5) * h;
		    // double y = (j + 0.5) * h;
		    (*u_mfi)(IntVect(i,j)) = 0.0;
		    (*u_dmfi)(IntVect(i,j)) = x;
		}
	    }
	}
    }
}

void
rz_adj(PArray<MultiFab> u[], PArray<MultiFab>& rhs,
       PArray<MultiFab>& rhoinv, const Array<BoxArray>& m,
       const Array<Box>& d)
{
    for (int ilev = 0 ; ilev < m.length() ; ilev++) 
    {
	double h = 1.0 / d[ilev].length(0);
	// double pi = 3.14159265358979323846;
	for ( MultiFabIterator u_mfi(u[1][ilev]); u_mfi.isValid(); ++u_mfi)
	{
	    DependentMultiFabIterator r_dmfi(u_mfi, rhoinv[ilev]);
	    int igrid = u_mfi.index();
	    for (int i = m[ilev][igrid].smallEnd(0);
		 i <= m[ilev][igrid].bigEnd(0); i++) 
	    {
		for (int j = m[ilev][igrid].smallEnd(1) - 1;
		     j <= m[ilev][igrid].bigEnd(1); j++) 
		{
		    double x = (i + 0.5) * h;
		    // double y = (j + 0.5) * h;
		    // double x0 = i * h, x1 = x0 + h;
		    //u[0][ilev][igrid](IntVect(i,j)) = 0.0;
		    (*u_mfi)(IntVect(i,j)) *= x;
		    //u[1][ilev][igrid](IntVect(i,j)) = x;
		    if (j >= m[ilev][igrid].smallEnd(1))
			(*r_dmfi)(IntVect(i,j)) *= x;
		}
	    }
	}
    }
}
#endif

void
projtest(const Array<BoxArray>& m, Array<IntVect>& ratio, Array<Box>& domain)
{
    // Note:  For terrain and full problems, h is ignored.
    Geometry crse_geom(domain[0]);
  
    Real h[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	h[i] = 1;
    }
    
    RegType bc[BL_SPACEDIM][2];
    
    bc[0][0] = bc00;
    bc[1][0] = bc10;
    bc[0][1] = bc01;
    bc[1][1] = bc11;

#if BL_SPACEDIM==3
    bc[2][0] = bc20;
    bc[2][1] = bc21;
#endif

    // bc[1][0] = refWall;
    // bc[1][1] = refWall;
    // bc[0][0] = refWall;
    // bc[0][1] = refWall;
    // bc[2][0] = periodic;
    // bc[2][1] = periodic;
    // bc[1][0] = inflow;
    // bc[1][1] = outflow;
        
    PArray<MultiFab> u[BL_SPACEDIM];
    PArray<MultiFab> p, rhoinv, rhs;
    
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
	u[i].resize(m.length());
    }
    p.resize(m.length());
    rhoinv.resize(m.length());
    rhs.resize(m.length());
    
    for (int ilev = 0; ilev < m.length(); ilev++) 
    {
	const BoxArray& cmesh = m[ilev];
	BoxArray nmesh = cmesh;
	nmesh.convert(IndexType(IntVect::TheNodeVector()));
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    u[i].set(ilev, new MultiFab(cmesh, 1, 1));
	}
	p.set(ilev, new MultiFab(nmesh, 1, 1));
	if (hg_stencil == holy_grail_amr_multigrid::terrain)
	{
	    rhoinv.set(ilev, new MultiFab(cmesh, 2*BL_SPACEDIM-1, 0));
	    rhoinv[ilev].setVal(1.0);
	    if ( true )
	    {
		rhoinv[ilev].setVal(0.0, BL_SPACEDIM, BL_SPACEDIM-1);
	    }
	    else
	    {
#if (BL_SPACEDIM == 2)
		rhoinv[ilev].setVal(0.2, 2, 1);
#else
		rhoinv[ilev].setVal(0.2, 3, 1);
		rhoinv[ilev].setVal(0.5, 4, 1);
#endif
	    }
	}
	else if ( hg_stencil == holy_grail_amr_multigrid::full )
	{
	    rhoinv.set(ilev, new MultiFab(cmesh, BL_SPACEDIM, 0));
	    rhoinv[ilev].setVal(1.0);
	}
	else
	{
	    rhoinv.set(ilev, new MultiFab(cmesh, 1, 0));
	    rhoinv[ilev].setVal(1.0);
	}
	rhs.set(ilev, new MultiFab(nmesh, 1, 1));
	//rhs.set(ilev, new MultiFab(cmesh, 1, 1));
	rhs[ilev].setVal(0.0);
    }
    
    if (hg_stencil == holy_grail_amr_multigrid::terrain)
    {
	// Adjust sigmas using refinement ratio information.
	// Assume spacing on level 0 already incorporated into values assigned above.
	// (h is ignored.)
	IntVect rat = IntVect::TheUnitVector();
	for (int ilev = 1; ilev < m.length(); ilev++) 
	{
	    rat *= ratio[ilev-1];
#if (BL_SPACEDIM == 2)
	    rhoinv[ilev].mult(Real(rat[0]) / rat[1], 0, 1);
	    rhoinv[ilev].mult(Real(rat[1]) / rat[0], 1, 1);
	    // component 2 remains unchanged
#else
	    rhoinv[ilev].mult(Real(rat[0]) / (rat[1] * rat[2]), 0, 1);
	    rhoinv[ilev].mult(Real(rat[1]) / (rat[0] * rat[2]), 1, 1);
	    rhoinv[ilev].mult(Real(rat[2]) / (rat[0] * rat[1]), 2, 1);
	    rhoinv[ilev].mult(1.0 /rat[1]                     , 3, 1);
	    rhoinv[ilev].mult(1.0 /rat[0]                     , 4, 1);
#endif
	}
    }
    else if (hg_stencil == holy_grail_amr_multigrid::full)
    {
	// Adjust sigmas using refinement ratio information.
	// Assume spacing on level 0 already incorporated into values assigned above.
	// (h is ignored.)
	IntVect rat = IntVect::TheUnitVector();
	for (int ilev = 1; ilev < m.length(); ilev++) 
	{
	    rat *= ratio[ilev-1];
#if (BL_SPACEDIM == 2)
	    rhoinv[ilev].mult(Real(rat[0]) / rat[1], 0, 1);
	    rhoinv[ilev].mult(Real(rat[1]) / rat[0], 1, 1);
	    // component 2 remains unchanged
#else
	    rhoinv[ilev].mult(Real(rat[0]) / (rat[1] * rat[2]), 0, 1);
	    rhoinv[ilev].mult(Real(rat[1]) / (rat[0] * rat[2]), 1, 1);
	    rhoinv[ilev].mult(Real(rat[2]) / (rat[0] * rat[1]), 2, 1);
#endif
	}
    }
    init(u, p, m, ratio);
    for(int ilev = 0; ilev < m.length(); ilev++)
    {
	HG_TEST_NORM( rhoinv[ilev], "proj");
	HG_TEST_NORM(p[ilev], "proj");
	for ( int i = 0; i < BL_SPACEDIM; ++i)
	{
	    HG_TEST_NORM( u[i][ilev], "proj");
	}
    }
	    
    
#if (BL_SPACEDIM == 2)
    //hb93_test1(u, m);
    //linear_test(u, m, domain);
    //rhs[1][0](IntVect(16,51)) = 100.0;
    //rhs[1][0](IntVect(16,50)) = 100.0;
    //rhs[0][0](IntVect(24,32)) = 100.0;
    //rhs[1][0](IntVect(30,40)) = -100.0;
    //rhs[0][0](IntVect(4,13)) = 100.0;
    //rhs[0][0](IntVect(20,90)) = 100.0;
#else
    //rhs[1][0](IntVect(16,20,40)) = 100.0;
    //rhs[1][0](IntVect(16,20,21)) = 100.0;
#endif
    
    /*
    u[0].assign(0.0);
    u[1].assign(-980.0);
    
      //Box bb(0,27,63,36);
      Box bb(27,0,36,63);
      //Box bb(26,0,37,63);
      for (ilev = 0; ilev < m.length(); ilev++) 
      {
      for (int igrid = 0; igrid < m[ilev].length(); igrid++) 
      {
      rhoinv[ilev][igrid].assign(100.0, (bb & m[ilev][igrid]));
      }
      }
    */
    /*
    #if (BL_SPACEDIM == 2)
    //Box bb(6,6,11,11);
    Box bb(0,0,5,5);
    //Box bb(0,0,6,7);
    #else
    Box bb(IntVect(0,0,0),IntVect(8,8,8));
    #endif
    for (ilev = 0; ilev < m.length(); ilev++)
    {
    Box b = refine(bb, m[ilev].sig()/16);
    for (int igrid = 0; igrid < m[ilev].length(); igrid++) 
    {
    rhoinv[ilev][igrid].assign(1000.0, (b & m[ilev][igrid]));
    //rhoinv[ilev][igrid].assign(1.0, (b & m[ilev][igrid]));
    }
    }
    */
    /*
    Box bb(0,1,1,4);
    for (ilev = 0; ilev < m.length(); ilev++) 
    {
    Box b = refine(bb, m[ilev].sig()/8);
    for (int igrid = 0; igrid < m[ilev].length(); igrid++) 
    {
    rhoinv[ilev][igrid].assign(1000.0, (b & m[ilev][igrid]));
    }
    }
    */
    /* Layer
    //Box bb(0,1,1,4);
    Box bb(0,0,63,4);
    for (ilev = 0; ilev < m.length(); ilev++) 
    {
    Box b = refine(bb, m[ilev].sig()/32);
    for (int igrid = 0; igrid < m[ilev].length(); igrid++) 
    {
    rhoinv[ilev][igrid].assign(0.00001, (b & m[ilev][igrid]));
    }
    }
    */
    /* Drop
    Box bbb(16,6,19,9);
    for (ilev = 0; ilev < m.length(); ilev++) 
    {
    Box b = refine(bbb, m[ilev].sig()/32);
    for (int igrid = 0; igrid < m[ilev].length(); igrid++) 
    {
    rhoinv[ilev][igrid].assign(0.00001, (b & m[ilev][igrid]));
    }
    }
    */
    
    int sum = 0;
    if ( ParallelDescriptor::IOProcessor() )
    {
	cout << "Cells by level: ";
	for (int ilev = 0; ilev < m.length(); ilev++) 
	{
	    int lsum = 0;
	    for (int i = 0; i < m[ilev].length(); i++) 
	    {
		lsum += m[ilev][i].numPts();
	    }
	    cout << " " << lsum;
	    sum += lsum;
	}
	cout << "\nTotal cells:  " << sum << endl;
    }
    
    double t0, t1, t2;
    
#ifdef UNICOS
    //int pcode = 1, nrep = 8;
    // int pcode = 1, nrep = 1;
    // Real tol = 1.e-6;
    //Real tol = 2.e-10;
#else
    // int pcode = 4, nrep = 1;
    //Real tol = 1.e-14;
    //int pcode = 1, nrep = 3;
    //Real tol = 1.e-6;
    // for vd tests in May, and most code validation tests:
    // Real tol = 2.e-10;
    //Real tol = 5.e-9;
#endif
    t0 = Utility::second();
    inviscid_fluid_boundary afb(bc);
    holy_grail_amr_projector proj(m, ratio,
				  domain[m.length() - 1], 0,
				  m.length() - 1, m.length() - 1,
				  afb, hg_stencil, pcode);

#if BL_SPACEDIM == 2
    if (   coordsys == amr_multigrid::rz
	&& hg_stencil == holy_grail_amr_multigrid::cross)
    {
	rz_adj(u, rhs, rhoinv, m, domain);
	proj.setCoordSys(amr_multigrid::rz);
    }
#endif
    
    //proj.smoother_mode  = 1;
    //proj.line_solve_dim = BL_SPACEDIM - 1;
    //proj.line_solve_dim = -1;
    
    if (m.length() == 1) 
    {
	t1 = Utility::second();
	proj.project(u, p, null_amr_real, rhoinv,
		     0, 0, crse_geom,
		     h, tol);
	for (int i = 1; i < nrep; i++) 
	{
	    init(u, p, m, ratio);
	    proj.project(u, p, null_amr_real, rhoinv,
			 0, 0, crse_geom,
			 h, tol);
	}
	t2 = Utility::second();
	for(int i = 0; i < BL_SPACEDIM; ++i )
	{
	    HG_TEST_NORM( u[i][0], "proj");
	}
	if ( ParallelDescriptor::IOProcessor() )
	{
	    cout << "Init time was " << t1 - t0 << endl;
	    cout << "Proj time was " << t2 - t1 << endl;
	    cout << "Speed was " << double(t2 - t1) / (nrep * sum) << endl;
	}
	if ( false ) 
	{
	    cout << setprecision(16);
	    for(int i = 0; i < BL_SPACEDIM; ++i)
	    {
		cout << "uvm"[i] << "min = " << u[i][0].min(0)
		     << ", " << "uvw"[i] << "max = " << u[i][0].max(0) << endl;
	    }
	    cout << setprecision(6);
	}
    }
    else if (m.length() == 2) 
    {
	//proj.project(u, p, null_amr_real, rhoinv, h, tol, 0, 0);
	//proj.project(u, p, null_amr_real, rhoinv, h, tol, 1, 1);
	for (int i = 0; i < p.length(); i++)
	{
	    p[i].setVal(0.0);
	}
	t1 = Utility::second();
	proj.project(u, p, null_amr_real, rhoinv,
		     0, 0, crse_geom,
		     h, tol, 0, 1);
	//proj.manual_project(u, p, rhs, null_amr_real, rhoinv, 1, h, tol, 0, 1);
	for (int i = 1; i < nrep; i++) 
	{
	    init(u, p, m, ratio);
	    proj.project(u, p, null_amr_real, rhoinv,
			 0, 0, crse_geom,
			 h, tol, 0, 1);
	}
	t2 = Utility::second();
	for(int i = 0; i < BL_SPACEDIM; ++i )
	{
	    HG_TEST_NORM( u[i][0], "proj");
	}
	if ( ParallelDescriptor::IOProcessor() )
	{
	    cout << "First time is " << t1 - t0 << endl;
	    cout << "Second time is " << t2 - t1 << endl;
	    cout << "Sync speed was " << double(t2 - t1) / (nrep * sum) << endl;
	}
	/*
	for (i = m[1][0].smallEnd(1); i <= m[1][0].bigEnd(1)+1; i++) 
	{
	cout << p[1][0](IntVect(0, i)) << endl;
	}
	proj.project(u, p, null_amr_real, rhoinv, h, tol, 1, 1);
	for (i = m[1][0].smallEnd(1); i <= m[1][0].bigEnd(1)+1; i++) 
	{
	cout << p[1][0](IntVect(0, i)) << endl;
	}
	*/
    }
    else if (m.length() == 3) 
    {
        double t00, t01, t10, t11, t20, t21;
        init(u, p, m, ratio);
	t00 = Utility::second();
	proj.project(u, p, null_amr_real, rhoinv,
		     0, 0, crse_geom,
		     h, tol,
		     2, 2);
	t01 = Utility::second();
	for (int i = 0; i < p.length(); i++)
	{
	    p[i].setVal(0.0);
	}
        init(u, p, m, ratio);
	t10 = Utility::second();
	proj.project(u, p, null_amr_real, rhoinv,
		     0, 0, crse_geom,
		     h, tol,
		     1, 2);
	t11 = Utility::second();
	for (int i = 0; i < p.length(); i++)
	{
	    p[i].setVal(0.0);
	}
        init(u, p, m, ratio);
	t20 = Utility::second();
	proj.project(u, p, null_amr_real, rhoinv,
		     0, 0, crse_geom,
		     h, tol,
		     0, 2);
	t21 = Utility::second();
	for(int i = 0; i < BL_SPACEDIM; ++i )
	{
	    HG_TEST_NORM( u[i][0], "proj");
	}
	if ( ParallelDescriptor::IOProcessor() )
	{
	    cout << "First time is " << t01 - t00 << endl;
	    cout << "Second time is " << t11 - t10 << endl;
	    cout << "Third time is " << t21 - t20 << endl;
	    cout << "Total time (counting inits) was  " << t21 - t00 << endl;
	}
    }
    else 
    {
	proj.make_it_so();
	proj.manual_project(u, p, null_amr_real, rhs, rhoinv,
			    0, 0, crse_geom,
			    true,
			    h, tol, 0, 3);
	t1 = Utility::second();
	for(int i = 0; i < BL_SPACEDIM; ++i )
	{
	    HG_TEST_NORM( u[i][0], "proj");
	}
	cout << "First time is " << t1 - t0 << endl;
    }
    /*
    if (m.length() < 3) 
    {
    for (i = 0; i < p.length(); i++)
    p[i].setVal(0.0);
    holy_grail_amr_projector proj(m, 0, m.length() - 1, afb, pcode);
    proj.project(u, p, null_amr_real, rhoinv, h, 1.e-14);
    t2 = Utility::second();
    cout << "Second time is " << t2 - t1 << endl;
    cout << "Total time was  " << t2 - t0 << endl;
    }
    */
    for (int ilev = 0; ilev < m.length(); ilev++) 
    {
	for (int i = 0; i < BL_SPACEDIM; i++)
	{
	    delete u[i].remove(ilev);
	}
	delete rhoinv.remove(ilev);
	delete p.remove(ilev);
	delete rhs.remove(ilev);
    }
}
