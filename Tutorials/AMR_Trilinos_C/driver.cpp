//stand alone solver for opal space charges calculation

#include "ml_include.h"

//the following code cannot be compiled without these ML Trilinos packages
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

#include "Solver.H"

#include <Utility.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <MultiFab.H>
#include "Epetra_MpiComm.h"
#include <VisMF.H>
#include <writePlotFile.H>

#include <vector>
#include <string>
#include <iostream>

using namespace Teuchos;

void   DumpMF       (std::string filename, MultiFab& soln, Box& domain, std::vector<double> hr);
void   compute_rhs  (                      MultiFab&  rhs, Box& domain, std::vector<double> hr);

double x_center, y_center, z_center;

int main(int argc, char *argv[]) {

    BoxLib::Initialize(argc,argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    // These don't need default values because they must be read from the inputs file.
    int    nx, ny;
    double hx, hy; 
#if (BL_SPACEDIM == 3)
    int    nz;
    double hz;
#endif
    bool verbose;

    // Default values
    int maxiters = 1000;
    int numBlocks = 1;
    int recycleBlocks = 0;
    int maxOldLHS = 1;
    double tol = 1e-8;

    ParmParse pp;

    pp.get("nx", nx);
    pp.get("ny", ny);
#if (BL_SPACEDIM == 3)
    pp.get("nz", nz);
#endif

    pp.get("tol", tol); 
    pp.get("maxiters", maxiters); 
    pp.get("numBlocks", numBlocks);
    pp.get("recycleBlocks", recycleBlocks);
    pp.get("maxOldLHS", maxOldLHS);

    pp.get("verbose", verbose);

    // Define the problem size to be a [0,1] x [0,1] x [0,1] box.
    hx = 1.0 / double(nx);
    hy = 1.0 / double(ny);

    x_center = 0.5;
    y_center = 0.5;

#if (BL_SPACEDIM == 3)
    hz = 1.0 / double(nz);
    z_center = 0.5;
#endif

#if (BL_SPACEDIM == 2)
    std::vector<int> nr(2);
    std::vector<double> hr(2);
    nr[0] = nx; nr[1] = ny; 
    hr[0] = hx; hr[1] = hy; 
#elif (BL_SPACEDIM == 3)
    std::vector<int> nr(3);
    std::vector<double> hr(3);
    nr[0] = nx; nr[1] = ny; nr[2] = nz;
    hr[0] = hx; hr[1] = hy; hr[2] = hz;
#endif
 
    int max_grid_size;
    pp.get("max_grid_size",max_grid_size);
 
    // Define a single box covering the domain
    IntVect dom_lo(D_DECL(0,0,0));
    IntVect dom_hi(D_DECL(nr[0]-1,nr[1]-1,nr[2]-1));
    Box domain(dom_lo,dom_hi);
 
    // Initialize the boxarray "bs" from the single box "bx"
    BoxArray bs(domain);
 
    // Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
    bs.maxSize(max_grid_size);
  
    int ncomps = 1;
    int nghost = 0;
    MultiFab rhs(bs,ncomps,nghost);
    MultiFab soln(bs,ncomps,nghost);
 
    // Set the initial guess to 1
    soln.setVal(0.);
 
    compute_rhs(rhs,domain,hr);

    // This sets the boundary conditions to be periodic or not
    int is_per[BL_SPACEDIM];
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "BL_SPACEDIM : " << BL_SPACEDIM << std::endl; 
#if (BL_SPACEDIM == 2)
        std::cout << "nr = " << nr[0] << " " << nr[1] << std::endl;
        std::cout << "hr = " << hr[0] << " " << hr[1] << std::endl;
#elif (BL_SPACEDIM == 3)
        std::cout << "nr = " << nr[0] << " " << nr[1] << " " << nr[2] << std::endl;
        std::cout << "hr = " << hr[0] << " " << hr[1] << " " << hr[2] << std::endl;
#endif
        std::cout << "verbose: " << verbose << std::endl;
	std::cout << "Max_grid_size   : " << max_grid_size << std::endl;
        std::cout << "Number of grids : " << bs.size() << std::endl; 
    }

    Solver* s = new Solver(domain, hr, Comm, verbose, 
                           tol, maxiters, numBlocks, recycleBlocks, maxOldLHS, 
                           rhs, soln);

    s->Compute();
    s->CopySolution(domain,soln);
    DumpMF("SOLN",soln,domain,hr);
    DumpMF("RHS", rhs ,domain,hr);

    delete s;

    BoxLib::Finalize();
    exit(EXIT_SUCCESS);
}

void compute_rhs(MultiFab& rhs, Box& domain, std::vector<double> hr)
{
  // Define the RHS
    double rsq;

    const int* dom_lo = domain.loVect();
    const int* dom_hi = domain.hiVect();

    // Initialize the RHS to zero everywhere.
    rhs.setVal(0.0);

    Real radius_sq = 0.25*0.25;

    // Note that the values are NODE-based
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.validbox();
       const int* bx_lo = bx.loVect();
       const int* bx_hi = bx.hiVect();
       FArrayBox& fab = rhs[mfi.index()];
       for(int i = bx_lo[0]; i<= bx_hi[0]; i++)
          for(int j = bx_lo[1]; j<= bx_hi[1]; j++)
             if (i > dom_lo[0] && i < dom_hi[0] && j > dom_lo[1] && j < dom_hi[1])
             {
#if (BL_SPACEDIM == 3)
             for(int k = bx_lo[2]; k<= bx_hi[2]; k++) 
                if (k > dom_lo[2] && k < dom_hi[2])
#endif
                {
                    double xx = (i+0.5)*hr[0] - x_center;
                    double yy = (j+0.5)*hr[1] - y_center;
                    rsq = xx*xx+yy*yy;
#if (BL_SPACEDIM == 3)
                    double zz = (k+0.5)*hr[2] - z_center;
                    rsq = rsq + zz*zz;
#endif
                    IntVect cells(D_DECL(i,j,k));
                    if(rsq < radius_sq)
                      fab(cells,0) = 1.0;

                }
             }
    }
}

void 
DumpMF(std::string filename, MultiFab& mf, Box& domain, std::vector<double> hr)
{
    int output_type = 2;

    if (output_type == 1)
    {
       VisMF::Write(mf,filename);
    }
    else  if (output_type == 2)
    {
       // This is all in order to define a Geometry object which is needed for writing the plotfiles
       int coord = 0;
       int is_per[BL_SPACEDIM];
       for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 0;
        RealBox real_box;
        for (int n = 0; n < BL_SPACEDIM; n++) {
          real_box.setLo(n, 0.0);
          real_box.setHi(n, 1.0);
        }

       Geometry geom(domain,&real_box,coord,is_per);

       writePlotFile(filename, mf, geom);
    }
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO) */
