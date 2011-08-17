
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
#include <COEF_F.H>

#ifdef USE_F90_MG
#include <MacBndry.H>
#include <MGT_Solver.H>
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

static
BoxList
readBoxList(const std::string file, Box& domain)
{
    BoxList retval;

    std::ifstream boxspec;

    boxspec.open(file.c_str(), std::ios::in);

    if( !boxspec )
    {
        std::string msg = "readBoxList: unable to open ";
        msg += file;
        BoxLib::Error(msg.c_str());
    }
    boxspec >> domain;
    
    int numbox = 0;
    boxspec >> numbox;

    for ( int i=0; i<numbox; i++ )
    {
        Box tmpbox;
        boxspec >> tmpbox;
        if( !domain.contains(tmpbox) )
	{
            std::cerr << "readBoxList: bogus box " << tmpbox << '\n';
            exit(1);
        }
        retval.push_back(tmpbox);
    }

    return retval;
}

static
void
writePlotFile (const std::string& dir,
               const MultiFab&    mf,
               const Geometry&    geom)
{
    BL_ASSERT(mf.nComp() == 2);
    //
    // Only let 64 CPUs be writing at any one time.
    //
    VisMF::SetNOutFiles(64);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(dir, 0755))
            BoxLib::CreateDirectoryFailed(dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    std::string HeaderFileName = dir + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc|std::ios::binary);
        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);
        HeaderFile << "NavierStokes-V1.1\n";
        HeaderFile << 2 << '\n';
        HeaderFile << "soln\nrhs\n";
        HeaderFile << BL_SPACEDIM << '\n';
        HeaderFile << 0 << '\n';
        HeaderFile << 0 << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom.ProbLo(i) << ' ';
        HeaderFile << '\n';
        for (int i = 0; i < BL_SPACEDIM; i++)
            HeaderFile << geom.ProbHi(i) << ' ';
        HeaderFile << '\n';
        HeaderFile << '\n';
        HeaderFile << geom.Domain() << ' ';
        HeaderFile << '\n';
        HeaderFile << 0 << ' ';
        HeaderFile << '\n';
        for (int k = 0; k < BL_SPACEDIM; k++)
            HeaderFile << geom.CellSize()[k] << ' ';
        HeaderFile << '\n';
        HeaderFile << geom.Coord() << '\n';
        HeaderFile << "0\n";
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", 0);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        HeaderFile << 0 << ' ' << mf.boxArray().size() << ' ' << 0 << '\n';
        HeaderFile << 0 << '\n';

        for (int i = 0; i < mf.boxArray().size(); ++i)
        {
            RealBox loc = RealBox(mf.boxArray()[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < BL_SPACEDIM; n++)
                HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
        }

        std::string PathNameInHeader = Level;
        PathNameInHeader += BaseName;
        HeaderFile << PathNameInHeader << '\n';
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(mf,TheFullPath);
}

int
main (int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  std::cout << std::setprecision(15);

  ParmParse pp;
  //
  // Obtain prob domain and box-list, set H per phys domain [0:1]Xn
  //

  int n_cell;
  int max_grid_size;

  pp.get("n_cell",n_cell);
  pp.get("max_grid_size",max_grid_size);

  // Define a single box covering the domain
  IntVect dom_lo(0,0,0);
  IntVect dom_hi(n_cell-1,n_cell-1,n_cell-1);
  Box container(dom_lo,dom_hi);

  // Initialize the boxarray "bs" from the single box "bx"
  BoxArray bs(container);

  // Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
  bs.maxSize(max_grid_size);

  Geometry geom( container );
  Real dx[BL_SPACEDIM];
  for ( int n=0; n<BL_SPACEDIM; n++ )
      dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/container.length(n);

  //
  // Allocate/initialize solution and right-hand-side, reset
  // rhs=1 at each box center.
  //

  int Ncomp=1;

#ifdef USE_F90_MG
  int Nghost=1;
#else
  int Nghost=0;
#endif

  MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
  MultiFab  rhs(bs, Ncomp,      0, Fab_allocate);  rhs.setVal(0.0);

  // Initialize solution to zero
  soln.setVal(0.0);

  // Define rhs in Fortran routine.
  rhs.setVal(0.0);
  for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
  {
    int i = mfi.index();
    const int* rlo = rhs[mfi].loVect();
    const int* rhi = rhs[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_SET_RHS(rhs[mfi].dataPtr(),ARLIM(rlo),ARLIM(rhi),
                 bx.loVect(),bx.hiVect(),dx,geom.ProbLo(),geom.ProbHi());
  }

  //
  // Initialize boundary data, set boundary condition flags and locations:
  // (phys boundaries set to dirichlet on cell walls).
  //

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
  //
  // Choose operator (Laplacian or ABecLaplacian), get tolerance, numiter.
  //
  bool ABec=false           ; pp.query("ABec",ABec);
  Real tolerance = 1.0e-12  ; pp.query("tol", tolerance);
  Real tolerance_abs = -1.0 ; pp.query("tol_abs", tolerance_abs);
  int numiter = 41          ; pp.query("numiter", numiter);
  int maxiter = 40          ; pp.query("maxiter", maxiter);
  bool mg = true            ; pp.query("mg", mg);
  bool cg = false           ; pp.query("cg", cg);
  bool bicg = false         ; pp.query("bicg", bicg);
  bool use_mg_pre=false     ; pp.query("mg_pre",use_mg_pre);
  bool new_bc=false         ; pp.query("new_bc",new_bc);
  bool dump_norm=true       ; pp.query("dump_norm", dump_norm);
  bool dump_Lp=false        ; pp.query("dump_Lp",dump_Lp);
  bool dump_MF=false        ; pp.query("dump_MF", dump_MF);
  bool dump_VisMF=false     ; pp.query("dump_VisMF", dump_VisMF);
  bool use_fboxlib=false    ; pp.query("use_fboxlib", use_fboxlib);

  bool use_variable_coef=false; pp.query("use_variable_coef", use_variable_coef);

  int res;


  {
      //
      // Allocate space for ABecLapacian coeffs, fill with values.
      //
      Real alpha = 1.0; pp.query("alpha",alpha);
      Real beta =  1.0; pp.query("beta",beta);

      Real a=0.0; pp.query("a",  a);

      Tuple<Real, BL_SPACEDIM> b;

      b[0]=1.0; pp.query("b0", b[0]);
      b[1]=1.0; pp.query("b1", b[1]);
      b[2]=1.0; pp.query("b2", b[2]);
        
      MultiFab  acoefs;
      acoefs.define(bs, Ncomp, Nghost, Fab_allocate);
      acoefs.setVal(a);
        
      MultiFab bcoefs[BL_SPACEDIM];

      if (use_variable_coef) {
        MultiFab cc_coef(bs,1,1);
        for ( MFIter mfi(cc_coef); mfi.isValid(); ++mfi )
        {
          int i = mfi.index();
          const int* clo = cc_coef[mfi].loVect();
          const int* chi = cc_coef[mfi].hiVect();
          const Box& bx = mfi.validbox();
  
          FORT_SET_CC_COEF(cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),bx.loVect(),bx.hiVect(),dx,geom.ProbLo(),geom.ProbHi());
        }

        VisMF::Write(cc_coef,"COEF");

        for ( int n=0; n<BL_SPACEDIM; ++n )
        {
  	  BoxArray bsC(bs);
	  bcoefs[n].define(bsC.surroundingNodes(n), Ncomp, Nghost, Fab_allocate);
          for ( MFIter mfi(bcoefs[n]); mfi.isValid(); ++mfi )
          {
            int i = mfi.index();
            Box bx(bs[i]);
            const int* clo = cc_coef[mfi].loVect();
            const int* chi = cc_coef[mfi].hiVect();
            const int* edgelo = bcoefs[n][mfi].loVect();
            const int* edgehi = bcoefs[n][mfi].hiVect();

            FORT_COEF_TO_EDGES(&n,bcoefs[n][mfi].dataPtr(),ARLIM(edgelo),ARLIM(edgehi),
                               cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
                               bx.loVect(),bx.hiVect());
          }
        }
      } else {
        for ( int n=0; n<BL_SPACEDIM; ++n )
        {
  	  BoxArray bsC(bs);
	  bcoefs[n].define(bsC.surroundingNodes(n), Ncomp, Nghost, Fab_allocate);
  	  bcoefs[n].setVal(b[n]);
        }
      }

      //
      // Build operator, set coeffs, build solver, solve
      //
      if ( use_fboxlib )
      {
#ifdef USE_F90_MG
          // Translate into F90 solver
          std::vector<BoxArray> bav(1);
          bav[0] = rhs.boxArray();
          std::vector<DistributionMapping> dmv(1);
          dmv[0] = rhs.DistributionMap();
          std::vector<Geometry> fgeom(1);
          fgeom[0] = geom;

          // Hard-wired to periodic
          int mg_bc[2*BL_SPACEDIM];
          for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
          {
             mg_bc[2*dir + 0] = 0;
             mg_bc[2*dir + 1] = 0;
          }

          MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false);

          MultiFab* soln_p[1]; soln_p[0] = &soln;
          MultiFab*  rhs_p[1];  rhs_p[0] = &rhs;

          Array< PArray<MultiFab> > coeffs(1);
          coeffs[0].resize(BL_SPACEDIM,PArrayManage);

          Array< Array<Real> > xa(1);
          Array< Array<Real> > xb(1);

          xa[0].resize(BL_SPACEDIM);
          xb[0].resize(BL_SPACEDIM);

          for (int i = 0; i < BL_SPACEDIM ; i++) {

              BoxArray edge_boxes(rhs.boxArray());
              edge_boxes.surroundingNodes(i);
              coeffs[0].set(i, new MultiFab(edge_boxes,1,0,Fab_allocate));
              coeffs[0][i].setVal(-1.0);

              xa[0][i] = 0.;
              xb[0][i] = 0.;
          }

          mgt_solver.set_gravity_coefficients(coeffs,xa,xb,1);

          BCRec* phys_bc = new BCRec;
          for (int i = 0; i < BL_SPACEDIM; i++)
          {
              phys_bc->setLo(i,0);
              phys_bc->setHi(i,0);
          }

          MacBndry bndry(bs,1,geom);
          bndry.setBndryValues(soln,0,0,1,*phys_bc);

          const Real run_strt = ParallelDescriptor::second();
         
          Real final_resnorm;
          mgt_solver.solve(soln_p, rhs_p, tolerance, tolerance_abs, bndry, final_resnorm);

          const int IOProc   = ParallelDescriptor::IOProcessorNumber();
	  Real      run_stop = ParallelDescriptor::second() - run_strt;

	  ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

	  if (ParallelDescriptor::IOProcessor())
              std::cout << "Run time = " << run_stop << std::endl;
#endif
      }
      else 
      {
	  ABecLaplacian lp(bd, dx);
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
	      const Real run_strt = ParallelDescriptor::second();

	      MultiGrid mg(lp);
	      mg.solve(soln, rhs, tolerance, tolerance_abs);

	      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
	      Real      run_stop = ParallelDescriptor::second() - run_strt;

	      ParallelDescriptor::ReduceRealMax(run_stop,IOProc);

	      if (ParallelDescriptor::IOProcessor())
                  std::cout << "Run time = " << run_stop << std::endl;
          }

	  if ( dump_Lp )
              std::cout << lp << std::endl;
      }
  } // -->> solve D^2(soln)=rhs   or   (alpha*a - beta*D.(b.G))soln=rhs

  //
  // Write solution, and rhs.
  //
  if ( dump_norm )
  {
      double d1 = mfnorm_2_valid(soln);
      double d2 = mfnorm_0_valid(soln);
      if ( ParallelDescriptor::IOProcessor() )
      {
	  std::cout << "solution norm = " << d1 << "/" << d2 << std::endl;
      }

      if (false)
      {
          double mean = 0;
          for (MFIter mfi(soln); mfi.isValid(); ++mfi)
              mean += soln[mfi].sum(0);

          ParallelDescriptor::ReduceRealSum(mean);

          mean /= soln.boxArray().numPts();

          for (MFIter mfi(soln); mfi.isValid(); ++mfi)
              soln[mfi].plus(-mean);

          double d1 = mfnorm_2_valid(soln);
          double d2 = mfnorm_0_valid(soln);
          if ( ParallelDescriptor::IOProcessor() )
          {
              std::cout << "solution norm (w/mean subtracted off) = " << d1 << "/" << d2 << std::endl;
          }
      }
  }

  if ( dump_MF || dump_VisMF )
  {
      MultiFab temp(bs, 2, 0);
      temp.setVal(0.0);
      temp.copy(soln, 0, 0, 1);
      temp.copy(rhs,  0, 1, 1);
      if ( dump_MF )
      {
	  writePlotFile("soln_pf", temp, geom);
      }
      if ( dump_VisMF )
      {
	  VisMF::Write(temp, "soln_vismf", VisMF::OneFilePerCPU);
      }
  }
  
  BoxLib::Finalize();

}

