//
// $Id: main.cpp,v 1.21 2001-08-01 21:51:04 lijewski Exp $
//

#include <fstream>
#include <iomanip>

#include <Utility.H>
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

#include <WritePlotFile.H>
#ifdef MG_USE_HYPRE
#include <HypreABec.H>
#endif

static
Real
mfnorm_0_valid (const MultiFab& mf)
{
  Real r = 0;
  for ( MFIter cmfi(mf); cmfi.isValid(); ++cmfi )
    {
      Real s = mf[cmfi].norm(cmfi.validbox(), 0, 0, mf[cmfi].nComp());
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
  for ( MFIter cmfi(mf); cmfi.isValid(); ++cmfi )
    {
      Real s = mf[cmfi].norm(cmfi.validbox(), 2, 0, mf[cmfi].nComp());
      r += s*s;
    }
  ParallelDescriptor::ReduceRealSum(r);
  return ::sqrt(r);
}

BoxList readBoxList (std::string file, Box& domain);

int
main (int   argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  std::cout << std::setprecision(15);

  ParmParse pp;

#ifdef BL3_PTHREADS
  int maxthreads = 1; pp.query("maxthreads", maxthreads);
  wrkq.max_threads(maxthreads);
  std::cout << "maxthreads = " << wrkq.max_threads() << std::endl;
#endif

  // Obtain prob domain and box-list, set H per phys domain [0:1]Xn
  Box container;
#if (BL_SPACEDIM == 2)
  std::string boxfile("grids/gr.2_small_a") ; pp.query("boxes", boxfile);
#elif (BL_SPACEDIM == 3)
  std::string boxfile("grids/gr.3_small_a") ; pp.query("boxes", boxfile);
#endif
  BoxArray bs(readBoxList(boxfile,container));
  Geometry geom( container );
  Real H[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; n++ )
    {
      H[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/container.length(n);
    }
    
  // Allocate/initialize solution and right-hand-side, reset
  //  rhs=1 at each box center
  int Ncomp=1;
  int Nghost=0;
  MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
  MultiFab  rhs(bs, Ncomp, Nghost, Fab_allocate);  rhs.setVal(0.0);
  for ( MFIter rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi )
    {
      IntVect ivmid = (rhs[rhsmfi].smallEnd() + rhs[rhsmfi].bigEnd())/2;
      //ivmid -= IntVect(0, (rhsmfi().bigEnd()[1]-rhsmfi().smallEnd()[1])/2);
      //ivmid = rhsmfi().smallEnd();
      rhs[rhsmfi].operator()(ivmid,0) = 1;
      ivmid += IntVect::TheUnitVector();
      rhs[rhsmfi].operator()(ivmid,0) = -1;
      // rhsmfi->setVal(1.0);
    }

  // Initialize boundary data, set boundary condition flags and locations:
  // (phys boundaries set to dirichlet on cell walls)
  BndryData bd(bs, 1, geom);
  int comp = 0;
  for ( int n=0; n<BL_SPACEDIM; ++n )
    {
      for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
	{
	  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
	  bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
	  bd.setBoundLoc(Orientation(n, Orientation::high),i,0.0 );
	  bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	  bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
	  bd.setValue(Orientation(n, Orientation::low) ,i, 0.0);
	  bd.setValue(Orientation(n, Orientation::high),i, 0.0);
	}
    }

  // Choose operator (Laplacian or ABecLaplacian), get tolerance, numiter
  bool ABec=false           ; pp.query("ABec",ABec);
  bool Hypre=false          ; pp.query("Hypre", Hypre);
  Real tolerance = 1.0e-10  ; pp.query("tol", tolerance);
  Real tolerance_abs = -1.0 ; pp.query("tol_abs", tolerance_abs);
  int numiter = 41          ; pp.query("numiter", numiter);
  int maxiter = 40          ; pp.query("maxiter", maxiter);
  bool mg = true            ; pp.query("mg", mg);
  bool cg = false           ; pp.query("cg", cg);
  bool bicg = false         ; pp.query("bicg", bicg);
  bool acg = false          ; pp.query("acg", acg);
  bool use_mg_pre=false     ; pp.query("mg_pre",use_mg_pre);
  bool new_bc=false         ; pp.query("new_bc",new_bc);
  bool dump_norm=true       ; pp.query("dump_norm", dump_norm);
  bool dump_Lp=false        ; pp.query("dump_Lp",dump_Lp);
  bool dump_Mf=false        ; pp.query("dump_Mf", dump_Mf);
  bool dump_VisMF=false     ; pp.query("dump_VisMF", dump_VisMF);
  bool dump_ascii=false     ; pp.query("dump_ascii", dump_ascii);
  int res;
  if ( !ABec && !Hypre )
    {
      // Build Laplacian operator, solver, then solve 
      Laplacian lp(bd, H[0]);
      {
	double d = lp.norm();
	if ( ParallelDescriptor::IOProcessor() )
	  {
	    std::cout << "Norm = " << d << std::endl;
	  }
      }
      if ( mg )
	{
	  MultiGrid mg(lp);
	  mg.setNumIter(numiter);
	  mg.setMaxIter(maxiter);
	  mg.solve(soln, rhs, tolerance, tolerance_abs);
	  if ( new_bc )
	    {
	      for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
		{
		  int i = mfi.index(); //   ^^^ using rhs to get mfi.index() yes, this is a hack
		  for (int n=0; n<BL_SPACEDIM; ++n)
		    {
		      bd.setValue(Orientation(n, Orientation::low) ,i,2.0);
		      bd.setValue(Orientation(n, Orientation::high),i,2.0);
                    }
                }
	      lp.bndryData(bd);
	      mg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }
      if ( cg )
	{
	  CGSolver cg(lp,use_mg_pre); cg.setCGSolver(CGSolver::CG);
	  cg.setMaxIter(maxiter);
	  res = cg.solve(soln, rhs, tolerance, tolerance_abs);
	  std::cout << "CG Result = " << res << std::endl;
	  if ( new_bc )
	    {
	      for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
		{
		  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
		  for ( int n=0; n<BL_SPACEDIM; ++n )
		    {
		      bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
		      bd.setValue(Orientation(n, Orientation::high),i,4.0);
                    }
                }
	      lp.bndryData(bd);
	      res  = cg.solve(soln, rhs, tolerance, tolerance_abs);
	      std::cout << "CG (new_bc) Result = " << res << std::endl;
            }
        }
      if ( bicg )
	{
	  CGSolver cg(lp,use_mg_pre); cg.setCGSolver(CGSolver::BiCGStab);
	  cg.setMaxIter(maxiter);
	  res = cg.solve(soln, rhs, tolerance, tolerance_abs);
	  std::cout << "BiCGStab Result = " << res << std::endl;
	  if ( new_bc )
	    {
	      for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
		{
		  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
		  for ( int n=0; n<BL_SPACEDIM; ++n )
		    {
		      bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
		      bd.setValue(Orientation(n, Orientation::high),i,4.0);
                    }
                }
	      lp.bndryData(bd);
	      res = cg.solve(soln, rhs, tolerance, tolerance_abs);
	      std::cout << "BiCGStab (new_bc) Result = " << res << std::endl;
            }
        }
      if ( acg )
	{
	  CGSolver cg(lp,use_mg_pre); cg.setCGSolver(CGSolver::CG_Alt);
	  cg.setMaxIter(maxiter);
	  res = cg.solve(soln, rhs, tolerance, tolerance_abs);
	  std::cout << "aCG Result = " << res << std::endl;
	  if ( new_bc )
	    {
	      for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
		{
		  int i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
		  for ( int n=0; n<BL_SPACEDIM; ++n )
		    {
		      bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
		      bd.setValue(Orientation(n, Orientation::high),i,4.0);
                    } 
                } 
	      lp.bndryData(bd);
	      res = cg.solve(soln, rhs, tolerance, tolerance_abs);
	      std::cout << "aCG (new_bc) Result = " << res << std::endl;
            }
        }

      // Look at operator
      if ( dump_Lp )
	std::cout << lp << std::endl;
	
        
    }
  else
    {
      // Allocate space for ABecLapacian coeffs, fill with values
      Real alpha = 1.0; pp.query("alpha",alpha);
      Real beta =  1.0; pp.query("beta",beta);
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
	} 

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
          {
	    double d = lp.norm();
	    if ( ParallelDescriptor::IOProcessor() )
	      {
		std::cout << "Norm = " << d << std::endl;
	      }
          }

	  if ( mg )
	    {
	      MultiGrid mg(lp);
	      mg.setNumIter(numiter);
	      mg.setMaxIter(maxiter);
	      mg.solve(soln, rhs, tolerance, tolerance_abs);
	      if ( new_bc )
		{
		  for ( int i=0; i < bs.size(); ++i )
		    {
		      for ( int n=0; n<BL_SPACEDIM; ++n )
			{
			  bd.setValue(Orientation(n, Orientation::low) ,i,2.0);
			  bd.setValue(Orientation(n, Orientation::high),i,2.0);
			} 
		    }
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
		  for ( int i=0; i < bs.size(); ++i )
		    {
		      for ( int n=0; n<BL_SPACEDIM; ++n )
			{
			  bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
			  bd.setValue(Orientation(n, Orientation::high),i,4.0);
			}
		    }
		  lp.bndryData(bd);
		  cg.solve(soln, rhs, tolerance, tolerance_abs);
		}
	    }
	  if ( bicg )
	    {
	      CGSolver cg(lp,use_mg_pre); cg.setCGSolver(CGSolver::BiCGStab);
	      cg.setMaxIter(maxiter);
	      cg.solve(soln, rhs, tolerance, tolerance_abs);
	      if ( new_bc )
		{
		  for ( int i=0; i < bs.size(); ++i )
		    {
		      for ( int n=0; n<BL_SPACEDIM; ++n )
			{
			  bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
			  bd.setValue(Orientation(n, Orientation::high),i,4.0);
			}
		    }
		  lp.bndryData(bd);
		  cg.solve(soln, rhs, tolerance, tolerance_abs);
		}
	    }
	  if ( acg )
	    {
	      CGSolver cg(lp,use_mg_pre); cg.setCGSolver(CGSolver::CG_Alt);
	      cg.setMaxIter(maxiter);
	      cg.solve(soln, rhs, tolerance, tolerance_abs);
	      if ( new_bc )
		{
		  for ( int i=0; i < bs.size(); ++i )
		    {
		      for ( int n=0; n<BL_SPACEDIM; ++n )
			{
			  bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
			  bd.setValue(Orientation(n, Orientation::high),i,4.0);
			}
		    }
		  lp.bndryData(bd);
		  cg.solve(soln, rhs, tolerance, tolerance_abs);
		}
	    }

	  // Look at operator
	  if ( dump_Lp )
	    std::cout << lp << std::endl;
	}
    } // -->> solve D^2(soln)=rhs   or   (alpha*a - beta*D.(b.G))soln=rhs

  // Write solution, and rhs
  if ( dump_norm )
    {
      double d1 = mfnorm_2_valid(soln);
      double d2 = mfnorm_0_valid(soln);
      if ( ParallelDescriptor::IOProcessor() )
	{
	  std::cout << "solution norm = " << d1 << "/" << d2 << std::endl;
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
      for ( MFIter mfi(soln); mfi.isValid(); ++mfi )
	{
	  std::cout << soln[mfi] << std::endl;
	}
    }

  BoxLib::Finalize();

#ifdef BL3_PROFILING
  BL3_PROFILE_STOP(mg_main);
  BoxLib3::Profiler::Finalize();
#endif

}

BoxList
readBoxList(const std::string file, Box& domain)
{
  BoxList retval;
  std::ifstream boxspec(file.c_str());
  if( !boxspec )
    {
      BoxLib::Error("readBoxList: unable to open " + *file.c_str());
    }
  boxspec >> domain;
    
  int numbox = 0;
  boxspec >> numbox;

  for ( int i=0; i<numbox; i++ )
    {
      Box tmpbox;
      boxspec >> tmpbox;
      if( ! domain.contains(tmpbox))
	{
	  std::cerr << "readBoxList: bogus box " << tmpbox << '\n';
	  exit(1);
        }
      retval.append(tmpbox);
    }
  boxspec.close();
  return retval;
}

