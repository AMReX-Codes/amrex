//
// $Id: main.cpp,v 1.2 1998-04-01 00:32:27 car Exp $
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
#include <WriteMultiFab.H>
#include <ParallelDescriptor.H>

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

const int NPROCS = 1;

BoxList
readBoxList (aString file,
	     BOX&    domain);

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

    TRACER("mg");
    
    if(argc < 2) {
      cerr << "usage:  " << argv[0] << " inputsfile [options]" << '\n';
      exit(-1);
    }

    // Parse command line
    ParmParse pp(argc-2,argv+2,NULL,argv[1]); 

    int nprocs = NPROCS; pp.query("nprocs", nprocs);
#ifndef BL_USE_BSP
    if (nprocs > 1)
    {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
    }
#endif
    ParallelDescriptor::StartParallel(nprocs);

    int i, n;
    
      // Obtain prob domain and box-list, set H per phys domain [0:1]Xn
    BOX container;
#if (BL_SPACEDIM == 2)
    aString boxfile("grids/gr.2_small_a") ; pp.query("boxes", boxfile);
#endif
#if (BL_SPACEDIM == 3)
    aString boxfile("grids/gr.3_small_a") ; pp.query("boxes", boxfile);
#endif
    BoxArray bs(readBoxList(boxfile,container));
    Geometry geom( container );
    REAL H[BL_SPACEDIM];
    for (n=0; n<BL_SPACEDIM; n++) {
        H[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/container.length(n);
    } // -->> over dimension
    
      // Allocate/initialize solution and right-hand-side, reset
      //  rhs=1 at each box center
    int Ncomp=1;
    int Nghost=0;
    MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
    MultiFab  rhs(bs, Ncomp, Nghost, Fab_allocate);  rhs.setVal(0.0);
    //for (i=0; i < bs.length(); i++) {
    for(MultiFabIterator rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi) {
        INTVECT ivmid = (rhsmfi().smallEnd() + rhsmfi().bigEnd())/2;
        rhsmfi().operator()(ivmid,0) = 1;
        ivmid += IntVect::TheUnitVector();
        rhsmfi().operator()(ivmid,0) = -1;
    } // -->> over boxes in domain
    
      // Initialize boundary data, set boundary condition flags and locations:
      // (phys boundaries set to dirichlet on cell walls)
    BndryData bd(bs, 1, geom);
    for(n=0; n<BL_SPACEDIM; ++n) {
      //for(i=0; i < bs.length(); ++i) {
      for(MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi) {
        i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
            bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
            bd.setBoundLoc(Orientation(n, Orientation::high),i,0.0 );
            bd.setBoundCond(Orientation(n, Orientation::low) ,i,LO_DIRICHLET);
            bd.setBoundCond(Orientation(n, Orientation::high),i,LO_DIRICHLET);
            bd.setValue(Orientation(n, Orientation::low) ,i,1.0);
            bd.setValue(Orientation(n, Orientation::high),i,1.0);
      } // -->> over boxes in domain
    } // -->> over dimension

      // Choose operator (Laplacian or ABecLaplacian), get tolerance, numiter
    int ABec=0; pp.query("ABec",ABec);
    REAL tolerance = 1.0e-10; pp.query("tol", tolerance);
    REAL tolerance_abs = 1.0e-10; pp.query("tol_abs", tolerance_abs);
    int numiter = 41; pp.query("numiter", numiter);
    int maxiter = 40; pp.query("maxiter", maxiter);
    int mg = 1; pp.query("mg", mg);
    int cg = 0; pp.query("cg", cg);
    int mg_pre=0; pp.query("mg_pre",mg_pre);
    bool use_mg_pre = (mg_pre == 1 ? true : false );
    int new_bc=0; pp.query("new_bc",new_bc);
    if( !ABec ) {
          // Build Laplacian operator, solver, then solve 
        Laplacian lp(bs, bd, H[0]);
        if (mg) {
            MultiGrid mg(lp);
            mg.setNumIter(numiter);
            mg.setMaxIter(maxiter);
            mg.solve(soln, rhs, tolerance, tolerance_abs);
            if (new_bc) {
                //for (i=0; i < bs.length(); ++i) {
      for(MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi) {
        i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
                    for (n=0; n<BL_SPACEDIM; ++n) {
                        bd.setValue(Orientation(n, Orientation::low) ,i,2.0);
                        bd.setValue(Orientation(n, Orientation::high),i,2.0);
                    } // -->> over dimensions
                } // -->> over boxes in domain
                lp.bndryData(bd);
                mg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }
        if (cg) {
            CGSolver cg(lp,use_mg_pre);
            cg.setMaxIter(maxiter);
            cg.solve(soln, rhs, tolerance, tolerance_abs);
            if (new_bc) {
                //for (i=0; i < bs.length(); ++i) {
      for(MultiFabIterator mfi(rhs); mfi.isValid(); ++mfi) {
        i = mfi.index();  //   ^^^ using rhs to get mfi.index() yes, this is a hack
                    for (n=0; n<BL_SPACEDIM; ++n) {
                        bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
                        bd.setValue(Orientation(n, Orientation::high),i,4.0);
                    } // -->> over dimensions
                } // -->> over boxes in domain
                lp.bndryData(bd);
                cg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }
        
    } else {
          // Allocate space for ABecLapacian coeffs, fill with values
        REAL alpha=1.0; pp.query("alpha",alpha);
        REAL beta=-1.0; pp.query("beta",beta);
        REAL a=0.0; pp.query("a",  a);
        Tuple<REAL, BL_SPACEDIM> b;
        b[0]=1.0; pp.query("b0", b[0]);
        b[1]=1.0; pp.query("b1", b[1]);
#if (BL_SPACEDIM > 2)
        b[2]=1.0; pp.query("b2", b[2]);
#endif
        int new_b=0; pp.query("new_b",new_b);
        
        MultiFab  acoefs;
        acoefs.define(bs, Ncomp, Nghost, Fab_allocate);
        acoefs.setVal(a);
        
        MultiFab bcoefs[BL_SPACEDIM];
        for (n=0; n<BL_SPACEDIM; ++n) {
            BoxArray bsC(bs);
            bcoefs[n].define(bsC.surroundingNodes(n), Ncomp,
                             Nghost, Fab_allocate);
            bcoefs[n].setVal(b[n]);
        } // -->> over dimension
        
          // Build operator, set coeffs, build solver, solve
        ABecLaplacian lp(bs, bd, H);
        lp.setScalars(alpha, beta);
        lp.setCoefficients(acoefs, bcoefs);

        if (mg) {
            MultiGrid mg(lp);
            mg.setNumIter(numiter);
            mg.setMaxIter(maxiter);
            mg.solve(soln, rhs, tolerance, tolerance_abs);
            if (new_bc) {
                for (i=0; i < bs.length(); ++i) {
                    for (n=0; n<BL_SPACEDIM; ++n) {
                        bd.setValue(Orientation(n, Orientation::low) ,i,2.0);
                        bd.setValue(Orientation(n, Orientation::high),i,2.0);
                    } // -->> over dimensions
                } // -->> over boxes in domain
                lp.bndryData(bd);
                mg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }
        if (cg) {
            CGSolver cg(lp,use_mg_pre);
            cg.setMaxIter(maxiter);
            cg.solve(soln, rhs, tolerance, tolerance_abs);
            if (new_bc) {
                for (i=0; i < bs.length(); ++i) {
                    for (n=0; n<BL_SPACEDIM; ++n) {
                        bd.setValue(Orientation(n, Orientation::low) ,i,4.0);
                        bd.setValue(Orientation(n, Orientation::high),i,4.0);
                    } // -->> over dimensions
                } // -->> over boxes in domain
                lp.bndryData(bd);
                cg.solve(soln, rhs, tolerance, tolerance_abs);
            }
        }
    } // -->> solve D^2(soln)=rhs   or   (alpha*a - beta*D.(b.G))soln=rhs

    //ofstream os("pltfile");
    int ratio=2;
    REAL bg_val=1.0;
    FArrayBox::setFormat(FABio::FAB_NATIVE);
    BoxList bl(bs);
    if ( !ABec ) {
        WriteMultiFab("pltfile",soln, H[0], bl, container, ratio, bg_val);
    } else {
        WriteMultiFab("pltfile",soln, H, bl, container, ratio, bg_val);
    }

    EndParallel();
    
} // -->> main fnc


BoxList
readBoxList(const aString file, BOX& domain )
{
    BoxList retval;
    ifstream boxspec(file.c_str());
    if( !boxspec ) {
        BoxLib::Error("readBoxList: unable to open " + *file.c_str());
    }
    boxspec >> domain;
    
    int numbox;
    boxspec >> numbox;

    for(int i=0; i<numbox; i++) {
        BOX tmpbox;
        boxspec >> tmpbox;
        if( ! domain.contains(tmpbox)) {
            cerr << "readBoxList: bogus box " << tmpbox << '\n';
            exit(1);
        }
        retval.append(tmpbox);
    }
    boxspec.close();
    return retval;
}

