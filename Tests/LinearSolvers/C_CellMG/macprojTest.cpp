//
// $Id: macprojTest.cpp,v 1.2 2001-08-21 22:15:37 car Exp $
//
// In order to use this utility to mimic a mac projection solve from the
// full AMR code, the lines of code below need to be added to the top
// of mac_level_driver in MacOperator.cpp.  This will dump the required 
// information.  The data sent to cout needs to be set in the input and 
// grid files for this test utility.  The boundary data and multifabs will 
// be read directly.
//
//cout << grids << endl;
//cout << "use_cg_solve = " << use_cg_solve << endl;
//cout << "level = " << level << endl;
//cout << "Density = " << Density << endl;
//cout << "dx = " << dx[0] << " " << dx[1] << " " << dx[2] << endl;
//cout << "dt = " << dt << endl;
//cout << "mac_tol = " << mac_tol << endl;
//cout << "mac_abs_tol = " << mac_abs_tol << endl;
//cout << "rhs_scale = " << rhs_scale << endl;
//ofstream macOS;
//macOS.open("mac_bndry_OS",ios::out|ios::binary);
//mac_bndry.writeOn(macOS);
//macOS.close();
//
//
//

#include <Utility.H>
#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include <MacBndry.H>
#include <MultiGrid.H>
#include <CGSolver.H>
#include <Laplacian.H>
#include <MacOperator.H>
#include <ParallelDescriptor.H>
#include <VisMF.H>
#include <TV_TempWrite.H>
#include <WorkQueue.H>

#include <WritePlotFile.H>

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



BoxList readBoxList(aString file, 
                    Box& domain);

void mac_driver (const MacBndry& mac_bndry,
                 const BoxArray& grids,
                 int             use_cg_solve,
                 int             Density,
                 const Real*     dx,
                 Real            dt,
                 Real            mac_tol,
                 Real            mac_abs_tol,
                 Real            rhs_scale,
                 MultiFab*       area,
                 MultiFab&       volume,
                 MultiFab&       S,
                 MultiFab&       Rhs,
                 MultiFab*       u_mac,
                 MultiFab*       mac_phi);

int
main (int   argc, char* argv[])
{
    BoxLib::Initialize(argc, argv);

    //
    // Instantiate after we're running in Parallel.
    //
    
    //
    // Obtain prob domain and box-list, set dx
    //
    Box container;
#if (BL_SPACEDIM == 2)
    aString boxfile("grids/gr.2_mac_tst");
#elif (BL_SPACEDIM == 3)
    aString boxfile("grids/gr.3_mac_tst");
#endif
    ParmParse pp;
    pp.query("boxes", boxfile);

    BoxArray bs(readBoxList(boxfile, container));
    Geometry geom( container );

    Real dx[BL_SPACEDIM];
    for ( int n=0; n<BL_SPACEDIM; n++ )
    {
        dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/container.length(n);
    }
    
    //
    // Read the problem information
    //
    
    int use_cg_solve;
    if (!pp.query("use_cg_solve", use_cg_solve))
        BoxLib::Abort("Must specify use_cg_solve");

    int Density;
    if (!pp.query("Density", Density))
        BoxLib::Abort("Must specify Density");

    Real dt;
    if (!pp.query("dt", dt))
        BoxLib::Abort("Must specify dt");

    Real mac_tol;
    if (!pp.query("mac_tol", mac_tol))
        BoxLib::Abort("Must specify mac_tol");

    Real mac_abs_tol;
    if (!pp.query("mac_abs_tol", mac_abs_tol))
        BoxLib::Abort("Must specify mac_abs_tol");

    Real rhs_scale;
    if (!pp.query("rhs_scale", rhs_scale))
        BoxLib::Abort("Must specify rhs_scale");

    bool dump_norm = false;
    pp.query("dump_norm", dump_norm);

    //
    // Read the MultiFabs Defining the Problem
    //
    MultiFab volume, S, Rhs, mac_phi;
    MultiFab area[3], u_mac[3];
    MacBndry mac_bndry;

    readMF(area[0], "area0_MF");
    readMF(area[1], "area1_MF");
    readMF(area[2], "area2_MF");
    readMF(volume, "volume_MF");
    readMF(S, "S_MF");
    readMF(Rhs, "Rhs_MF");
    readMF(u_mac[0], "u_mac0_MF");
    readMF(u_mac[1], "u_mac1_MF");
    readMF(u_mac[2], "u_mac2_MF");
    readMF(mac_phi, "mac_phi_MF");

    ifstream macOS;
    macOS.open("mac_bndry_OS",ios::in|ios::binary);
    mac_bndry.readFrom(macOS);
    macOS.close();

    //
    // Solve System
    //
    mac_driver (mac_bndry, bs, use_cg_solve, Density, dx, dt,
                mac_tol, mac_abs_tol, rhs_scale, area, volume, S,
                Rhs, u_mac, &mac_phi);

    //
    // Write solution, and rhs
    //
    if ( dump_norm )
    {
        double d1 = mfnorm_2_valid(mac_phi);
        double d2 = mfnorm_0_valid(mac_phi);

        if ( ParallelDescriptor::IOProcessor() )
        {
            cout << "solution norm = " << d1 << "/" << d2 << endl;
        }
    }
  
    BoxLib::Finalize();
}

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


//
//  This routine is intended to mimic the mac_level_driver in MacOperator.cpp.
//  So, we can read in the appropriate MultiFabs and then pass them to this 
//  routine
//
void
mac_driver (const MacBndry& mac_bndry,
            const BoxArray& grids,
            int             use_cg_solve,
            int             Density,
            const Real*     dx,
            Real            dt,
            Real            mac_tol,
            Real            mac_abs_tol,
            Real            rhs_scale,
            MultiFab*       area,
            MultiFab&       volume,
            MultiFab&       S,
            MultiFab&       Rhs,
            MultiFab*       u_mac,
            MultiFab*       mac_phi)
{
    MacOperator mac_op(mac_bndry,dx);
    mac_op.setCoefficients(area,S,Density,dx);
    mac_op.defRHS(area,volume,Rhs,u_mac,rhs_scale);
    mac_op.maxOrder(2);

    if (use_cg_solve && mac_op.maxOrder() != 2)
    {
        BoxLib::Error("Can't use CGSolver with maxorder > 2");
    }
    //
    // Construct MultiGrid or CGSolver object and solve system.
    //
    if (use_cg_solve)
    {
        bool use_mg_precond = true;
        CGSolver mac_cg(mac_op,use_mg_precond);
        mac_cg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
    else
    {
        MultiGrid mac_mg(mac_op);
        mac_mg.solve(*mac_phi,Rhs,mac_tol,mac_abs_tol);
    }
}


