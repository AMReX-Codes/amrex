//
// $Id: main.cpp,v 1.13 2000-08-02 16:11:53 car Exp $
//

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
#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include <BndryData.H>
#include <MultiGrid.H>
#include <CGSolver.H>
#include <Laplacian.H>
#include <ABecLaplacian.H>
#include <ParallelDescriptor.H>
#include <VisMF.H>

#ifdef BL3_PTHREADS
#include <BoxLib3/WorkQueue.H>
BoxLib3::WorkQueue wrkq;
#endif
#ifdef BL3_PROFILING
#include <BoxLib3/Parallel.H>
#include <BoxLib3/Timer.H>
#endif


#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

#include <WritePlotFile.H>
#ifdef MG_USE_HYPRE
#include <HypreABec.H>
#endif

static
Real
mfnorm_0_valid (const MultiFab& mf)
{
    Real r = 0;
    for ( ConstMultiFabIterator cmfi(mf); cmfi.isValid(); ++cmfi )
    {
	Real s = cmfi->norm(cmfi.validbox(), 0, 0, cmfi->nComp());
	r = (r > s) ? r : s;
    }
    ParallelDescriptor::ReduceRealMax(r);
    return r;
}

static
Real
mfnorm_2_valid (const MultiFab& mf)
{
    Real r = 0;
    for ( ConstMultiFabIterator cmfi(mf); cmfi.isValid(); ++cmfi )
    {
	Real s = cmfi->norm(cmfi.validbox(), 2, 0, cmfi->nComp());
	r += s*s;
    }
    ParallelDescriptor::ReduceRealSum(r);
    return ::sqrt(r);
}

BoxList readBoxList (aString file, BOX& domain);

int
main (int   argc, char* argv[])
{
  //
  // Make sure to catch new failures.
  //
#ifndef WIN32
  set_new_handler(Utility::OutOfMemory);
#endif
#ifdef BL3_PROFILING
  BoxLib3::Profiler::Initialize(argc, argv);
  BL3_PROFILE_TIMER(mg_main, "main()");
  BL3_PROFILE_START(mg_main);
#endif
#ifdef BL3_USE_MPI
  BoxLib3::Parallel::Initialize(argc, argv);
#endif
  ParallelDescriptor::StartParallel(&argc, &argv);

  cout << setprecision(10);

  if ( argc < 2 )
    {
      cerr << "usage:  " << argv[0] << " inputsfile [options]" << '\n';
      exit(-1);
    }

  ParmParse pp(argc-2,argv+2,NULL,argv[1]); 
  //
  // Initialize random seed after we're running in parallel.
  //
  Utility::InitRandom(ParallelDescriptor::MyProc() + 1);
  //
  // Instantiate after we're running in Parallel.
  //
#ifndef WIN32
  int sleeptime = 0; pp.query("sleep", sleeptime);
  sleep(sleeptime);
#endif
#ifdef BL3_PTHREADS
  int maxthreads = 1; pp.query("maxthreads", maxthreads);
  wrkq.max_threads(maxthreads);
#endif
    
  TRACER("mg");
    
  // Obtain prob domain and box-list, set H per phys domain [0:1]Xn
  BOX container;
#if (BL_SPACEDIM == 2)
  aString boxfile("grids/gr.2_small_a") ; pp.query("boxes", boxfile);
#elif (BL_SPACEDIM == 3)
  aString boxfile("grids/gr.3_small_a") ; pp.query("boxes", boxfile);
#endif
  BoxArray bs(readBoxList(boxfile,container));
  Geometry geom( container );
  Real H[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; n++ )
    {
      H[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/container.length(n);
    } // -->> over dimension
    
  // Allocate/initialize solution and right-hand-side, reset
  //  rhs=1 at each box center
  int Ncomp=1;
  int Nghost=0;
  MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
  MultiFab  rhs(bs, Ncomp, Nghost, Fab_allocate);  rhs.setVal(0.0);
  for ( MultiFabIterator rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi )
    {
      INTVECT ivmid = (rhsmfi().smallEnd() + rhsmfi().bigEnd())/2;
      //ivmid -= IntVect(0, (rhsmfi().bigEnd()[1]-rhsmfi().smallEnd()[1])/2);
      //ivmid = rhsmfi().smallEnd();
      rhsmfi().operator()(ivmid,0) = 1;
      ivmid += IntVect::TheUnitVector();
      rhsmfi().operator()(ivmid,0) = -1;
    } // -->> over boxes in domain

  // Initialize boundary data, set boundary condition flags and locations:
  // (phys boundaries set to dirichlet on cell walls)
  BndryData bd(bs, 1, geom);
  int comp = 0;
  for ( int n=0; n<BL_SPACEDIM; ++n )
    {
      for ( MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi )
	{
	  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
	  bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
	  bd.setBoundLoc(Orientation(n, Orientation::high),i,0.0 );
	  bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	  bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
	  bd.setValue(Orientation(n, Orientation::low) ,i,1.0);
	  bd.setValue(Orientation(n, Orientation::high),i,1.0);
	} // -->> over boxes in domain
    } // -->> over dimension

  // Choose operator (Laplacian or ABecLaplacian), get tolerance, numiter
  bool ABec=false; pp.query("ABec",ABec);
  bool Hypre=false; pp.query("Hypre", Hypre);
  Real tolerance = 1.0e-10; pp.query("tol", tolerance);
  Real tolerance_abs = 1.0e-10; pp.query("tol_abs", tolerance_abs);
  int numiter = 41; pp.query("numiter", numiter);
  int maxiter = 40; pp.query("maxiter", maxiter);
  bool mg = true; pp.query("mg", mg);
  bool cg = false; pp.query("cg", cg);
  bool use_mg_pre=false; pp.query("mg_pre",use_mg_pre);
  bool new_bc=false; pp.query("new_bc",new_bc);
  bool dump_norm=true; pp.query("dump_norm", dump_norm);
  bool dump_Lp=false; pp.query("dump_Lp",dump_Lp);
  bool dump_Mf=false; pp.query("dump_Mf", dump_Mf);
  bool dump_VisMF=false; pp.query("dump_VisMF", dump_VisMF);
  bool dump_ascii=false; pp.query("dump_ascii", dump_ascii);
  if ( !ABec && !Hypre )
    {
      // Build Laplacian operator, solver, then solve 
      Laplacian lp(bd, H[0]);
      if ( mg )
	{
	  MultiGrid mg(lp);
	  mg.setNumIter(numiter);
	  mg.setMaxIter(maxiter);
	  mg.solve(soln, rhs, tolerance, tolerance_abs);
	  if ( new_bc )
	    {
	      for ( MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi )
		{
		  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
		  for (int n=0; n<BL_SPACEDIM; ++n)
		    {
		      bd.setValue(Orientation(n, Orientation::low) ,i,2.0);
		      bd.setValue(Orientation(n, Orientation::high),i,2.0);
                    } // -->> over dimensions
                } // -->> over boxes in domain
	      lp.bndryData(bd);
	      mg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }
      if ( cg )
	{
	  CGSolver cg(lp,use_mg_pre);
	  cg.setMaxIter(maxiter);
	  cg.solve(soln, rhs, tolerance, tolerance_abs);
	  if ( new_bc )
	    {
	      for ( MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi )
		{
		  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
		  for ( int n=0; n<BL_SPACEDIM; ++n )
		    {
		      bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
		      bd.setValue(Orientation(n, Orientation::high),i,4.0);
                    } // -->> over dimensions
                } // -->> over boxes in domain
	      lp.bndryData(bd);
	      cg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }

      // Look at operator
      if ( dump_Lp )
	cout << lp << endl;
	
        
    }
  else
    {
      // Allocate space for ABecLapacian coeffs, fill with values
      Real alpha = 1.0; pp.query("alpha",alpha);
      Real beta = -1.0; pp.query("beta",beta);
      Real a=0.0; pp.query("a",  a);
      Tuple<Real, BL_SPACEDIM> b;
      b[0]=1.0; pp.query("b0", b[0]);
      b[1]=1.0; pp.query("b1", b[1]);
#if (BL_SPACEDIM > 2)
      b[2]=1.0; pp.query("b2", b[2]);
#endif
        
      MultiFab  acoefs;
      acoefs.define(bs, Ncomp, Nghost, Fab_allocate);
      acoefs.setVal(a);
        
      MultiFab bcoefs[BL_SPACEDIM];
      for ( int n=0; n<BL_SPACEDIM; ++n )
	{
	  BoxArray bsC(bs);
	  bcoefs[n].define(bsC.surroundingNodes(n), Ncomp,
			   Nghost, Fab_allocate);
	  bcoefs[n].setVal(b[n]);
	} // -->> over dimension
        
      // Build operator, set coeffs, build solver, solve
      if ( Hypre )
	{
#ifdef MG_USE_HYPRE
	  ParmParse pp("hy");
	  int solver_flag = 0; pp.query("solver", solver_flag);
	  bool inhom = true; pp.query("inhom", inhom);
	  HypreABec hp(bs, bd, H, solver_flag, false);
	  hp.setScalars(alpha, beta);
	  hp.setCoefficients(acoefs, bcoefs);
	  hp.setup_solver(tolerance, tolerance_abs, maxiter);
	  hp.solve(soln, rhs, inhom);
	  hp.clear_solver();
#else
	  BoxLib::Error("No Hypre in this build");
#endif
	}
      else
	{
	  ABecLaplacian lp(bd, H);
	  lp.setScalars(alpha, beta);
	  lp.setCoefficients(acoefs, bcoefs);

	  if ( mg )
	    {
	      MultiGrid mg(lp);
	      mg.setNumIter(numiter);
	      mg.setMaxIter(maxiter);
	      mg.solve(soln, rhs, tolerance, tolerance_abs);
	      if ( new_bc )
		{
		  for ( int i=0; i < bs.length(); ++i )
		    {
		      for ( int n=0; n<BL_SPACEDIM; ++n )
			{
			  bd.setValue(Orientation(n, Orientation::low) ,i,2.0);
			  bd.setValue(Orientation(n, Orientation::high),i,2.0);
			} // -->> over dimensions
		    } // -->> over boxes in domain
		  lp.bndryData(bd);
		  mg.solve(soln, rhs, tolerance, tolerance_abs);
		}
	    }
	  if ( cg )
	    {
	      CGSolver cg(lp,use_mg_pre);
	      cg.setMaxIter(maxiter);
	      cg.solve(soln, rhs, tolerance, tolerance_abs);
	      if ( new_bc )
		{
		  for ( int i=0; i < bs.length(); ++i )
		    {
		      for ( int n=0; n<BL_SPACEDIM; ++n )
			{
			  bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
			  bd.setValue(Orientation(n, Orientation::high),i,4.0);
			} // -->> over dimensions
		    } // -->> over boxes in domain
		  lp.bndryData(bd);
		  cg.solve(soln, rhs, tolerance, tolerance_abs);
		}
	    }

	  // Look at operator
	  if ( dump_Lp )
	    cout << lp << endl;
	}
    } // -->> solve D^2(soln)=rhs   or   (alpha*a - beta*D.(b.G))soln=rhs

  // Write solution, and rhs
  if ( dump_norm )
  {
    double d1 = mfnorm_2_valid(soln);
    double d2 = mfnorm_0_valid(soln);
    if ( ParallelDescriptor::IOProcessor() )
      {
	cout << "solution norm = " << d1 << "/" << d2 << endl;
      }
  }
  if ( dump_Mf || dump_VisMF )
    {
      MultiFab temp(bs, 2, Nghost, Fab_allocate);
      temp.setVal(0.0);
      temp.copy(soln, 0, 0, 1);
      temp.copy(rhs,  0, 1, 1);
      if ( dump_Mf )
	{
	  writePlotFile("soln", temp, geom, IntVect(D_DECL(2,2,2)), temp.min(0));
	}
      if ( dump_VisMF )
	{
	  VisMF::Write(temp, "soln_vismf", VisMF::OneFilePerCPU);
	}
    }
  
  if ( dump_ascii )
    {
      for ( MultiFabIterator mfi(soln); mfi.isValid(); ++mfi )
	{
	  cout << *mfi << endl;
	} // -->> over boxes in domain
    }

#ifdef BL3_PROFILING
  BL3_PROFILE_STOP(mg_main);
  BoxLib3::Profiler::Finalize();
#endif
  ParallelDescriptor::EndParallel();
} // -->> main fnc

BoxList
readBoxList(const aString file, BOX& domain)
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

  for ( int i=0; i<numbox; i++ )
    {
      BOX tmpbox;
      boxspec >> tmpbox;
      if( ! domain.contains(tmpbox))
	{
	  cerr << "readBoxList: bogus box " << tmpbox << '\n';
	  exit(1);
        }
      retval.append(tmpbox);
    }
  boxspec.close();
  return retval;
}

