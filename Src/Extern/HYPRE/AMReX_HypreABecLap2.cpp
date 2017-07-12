#include <AMReX_HypreABecLap2.H>
#include <AMReX_HypreABec_F.H>

namespace {
    static int ispow2(int i)
    {
        return (i == 1) ? 1 : (((i <= 0) || (i & 1)) ? 0 : ispow2(i / 2));
    }
}

namespace amrex
{

HypreABecLap2::HypreABecLap2 (const BoxArray& grids,
                              const DistributionMapping& dmap,
                              const Geometry& geom_,
                              MPI_Comm comm_)
    : comm(comm_),
      geom(geom_),
      verbose(1),
      hgrid(NULL), stencil(NULL), graph(NULL),
      A(NULL), A0(NULL), b(NULL), x(NULL)
{
    const int ncomp = 1;
    const int ngrow = 0;
    acoefs.define(grids, dmap, ncomp, ngrow);
    acoefs.setVal(0.0);

    for (int i = 0; i < BL_SPACEDIM; ++i) {
        BoxArray edge_boxes(grids);
        edge_boxes.surroundingNodes(i);
        bcoefs[i].define(edge_boxes, dmap, ncomp, ngrow);        
        bcoefs[i].setVal(0.0);
    }

    bd.define(grids, dmap, ncomp, geom);

    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs );
    MPI_Comm_rank(comm, &myid );
 
    HYPRE_SStructGridCreate(comm, BL_SPACEDIM, 1, &hgrid);

    if (geom.isAnyPeriodic()) {
        int is_periodic[BL_SPACEDIM];
        for (int i = 0; i < BL_SPACEDIM; i++) {
            is_periodic[i] = 0;
            if (geom.isPeriodic(i)) {
                is_periodic[i] = geom.period(i);
                AMREX_ASSERT(ispow2(is_periodic[i]));
                AMREX_ASSERT(geom.Domain().smallEnd(i) == 0);
            }
        }
        HYPRE_SStructGridSetPeriodic(hgrid, 0, is_periodic);
    }

    for (int i = 0; i < grids.size(); i++) {
        if (dmap[i] == myid) {
            const Box& bx = grids[i];
            HYPRE_SStructGridSetExtents(hgrid, 0,
                                        const_cast<int*>(bx.loVect()),
                                        const_cast<int*>(bx.hiVect()));
        }
    }

    // All variables are cell-centered
    HYPRE_SStructVariable vars[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
    HYPRE_SStructGridSetVariables(hgrid, 0, 1, vars);

    HYPRE_SStructGridAssemble(hgrid);

    // Setup stencils
    static_assert(BL_SPACEDIM > 1, "HypreABecLap2: 1D not supported");
#if (BL_SPACEDIM == 2)
    int offsets[5][2] = {{ 0,  0},
                         {-1,  0},
                         { 1,  0},
                         { 0, -1},
                         { 0,  1}};
#elif (BL_SPACEDIM == 3)
    int offsets[7][3] = {{ 0,  0,  0},
                         {-1,  0,  0},
                         { 1,  0,  0},
                         { 0, -1,  0},
                         { 0,  1,  0},
                         { 0,  0, -1},
                         { 0,  0,  1}};
#endif

    HYPRE_SStructStencilCreate(BL_SPACEDIM, 2 * BL_SPACEDIM + 1, &stencil);

    for (int i = 0; i < 2 * BL_SPACEDIM + 1; i++) {
        HYPRE_SStructStencilSetEntry(stencil, i, offsets[i], 0);
    }

    HYPRE_SStructGraphCreate(comm, hgrid, &graph);
    HYPRE_SStructGraphSetObjectType(graph, HYPRE_PARCSR);

    HYPRE_SStructGraphSetStencil(graph, 0, 0, stencil);

    HYPRE_SStructGraphAssemble(graph);

    HYPRE_SStructMatrixCreate(comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_SStructMatrixInitialize(A);

    HYPRE_SStructMatrixCreate(comm, graph, &A0);
    HYPRE_SStructMatrixSetObjectType(A0, HYPRE_PARCSR);
    HYPRE_SStructMatrixInitialize(A0);

    HYPRE_SStructVectorCreate(comm, hgrid, &b);
    HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
    
    HYPRE_SStructVectorCreate(comm, hgrid, &x);
    HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);

    HYPRE_SStructVectorInitialize(b);
    HYPRE_SStructVectorInitialize(x);

    HYPRE_SStructVectorAssemble(b);
    HYPRE_SStructVectorAssemble(x);
}

HypreABecLap2::~HypreABecLap2 ()
{
    HYPRE_SStructVectorDestroy(b);
    HYPRE_SStructVectorDestroy(x);
    
    HYPRE_SStructMatrixDestroy(A);
    HYPRE_SStructMatrixDestroy(A0);
    
    HYPRE_SStructGraphDestroy(graph);
    HYPRE_SStructStencilDestroy(stencil);
    HYPRE_SStructGridDestroy(hgrid);
}

void
HypreABecLap2::setScalars (Real sa, Real sb)
{
    scalar_a = sa;
    scalar_b = sb;
}

void
HypreABecLap2::setACoeffs (const MultiFab& alpha)
{
    MultiFab::Copy(acoefs, alpha, 0, 0, 1, 0);
}

void
HypreABecLap2::setBCoeffs (const std::array<const MultiFab*,BL_SPACEDIM>& beta)
{
    for (int idim=0; idim<BL_SPACEDIM; idim++) {
        MultiFab::Copy(bcoefs[idim], *beta[idim], 0, 0, 1, 0);
    }
}

void
HypreABecLap2::setVerbose (int _verbose)
{
    verbose = _verbose;
}

void
HypreABecLap2::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                      int max_iter, LinOpBCType bc_type, Real bc_value)
{
    loadBndryData(bc_type, bc_value);

    loadMatrix();

    // finalize maxtrix

    amrex::Abort("solver to do");
}

void
HypreABecLap2::loadBndryData (LinOpBCType bc_type, Real bc_value)
{
    const int comp = 0;
    const Real* dx = geom.CellSize();
    for (int n=0; n<BL_SPACEDIM; ++n) {
        for (MFIter mfi(acoefs); mfi.isValid(); ++mfi ) {
            int i = mfi.index(); 
            
            const Box& bx = mfi.validbox();
            
            // Our default will be that the face of this grid is either touching another grid
            //  across an interior boundary or a periodic boundary.  We will test for the other
            //  cases below.
            {
                // Define the type of boundary conditions to be Dirichlet (even for periodic)
                bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
                bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
                
                // Set the boundary conditions to the cell centers outside the domain
                bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
                bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
            }
            
            // Now test to see if we should override the above with Dirichlet or Neumann physical bc's
            if (bc_type != LinOpBCType::interior)
            { 
                int ibnd = static_cast<int>(bc_type); // either LO_DIRICHLET or LO_NEUMANN
                
                // We are on the low side of the domain in coordinate direction n
                if (bx.smallEnd(n) == geom.Domain().smallEnd(n)) {
                    // Set the boundary conditions to live exactly on the faces of the domain
                    bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
                    
                    // Set the Dirichlet/Neumann boundary values 
                    bd.setValue(Orientation(n, Orientation::low) ,i, bc_value);
                    
                    // Define the type of boundary conditions 
                    bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,ibnd);
                }
                
                // We are on the high side of the domain in coordinate direction n
                if (bx.bigEnd(n) == geom.Domain().bigEnd(n)) {
                    // Set the boundary conditions to live exactly on the faces of the domain
                    bd.setBoundLoc(Orientation(n, Orientation::high) ,i,0.0 );
                    
                    // Set the Dirichlet/Neumann boundary values
                    bd.setValue(Orientation(n, Orientation::high) ,i, bc_value);
                    
                    // Define the type of boundary conditions 
                    bd.setBoundCond(Orientation(n, Orientation::high) ,i,comp,ibnd);
                }
            }
        }
    }
}

void
HypreABecLap2::loadMatrix ()
{
    static_assert(BL_SPACEDIM > 1, "HypreABecLap2: 1D not supported");
#if (BL_SPACEDIM == 2)
    int offsets[5][2] = {{ 0,  0},
                         {-1,  0},
                         { 1,  0},
                         { 0, -1},
                         { 0,  1}};
#elif (BL_SPACEDIM == 3)
    int offsets[7][3] = {{ 0,  0,  0},
                         {-1,  0,  0},
                         { 1,  0,  0},
                         { 0, -1,  0},
                         { 0,  1,  0},
                         { 0,  0, -1},
                         { 0,  0,  1}};
#endif

    const int size = 2 * BL_SPACEDIM + 1;

    int stencil_indices[size];    
    for (int i = 0; i < size; i++) {
        stencil_indices[i] = i;
    }
    
    const int part = 0;
    const int bho = 0;
    const Real* dx = geom.CellSize();

    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        int i = mfi.index();
        const Box &reg = mfi.validbox();

        int volume = reg.numPts();
        Real* mat = hypre_CTAlloc(double, size*volume);

        // build matrix interior

        hmac(mat,
             BL_TO_FORTRAN(acoefs[mfi]),
             ARLIM(reg.loVect()), ARLIM(reg.hiVect()), scalar_a);
        
        for (int idim = 0; idim < BL_SPACEDIM; idim++) {
            hmbc(mat, 
                 BL_TO_FORTRAN(bcoefs[idim][mfi]),
                 ARLIM(reg.loVect()), ARLIM(reg.hiVect()), scalar_b,
                 geom.CellSize(), idim);
        }

        // add b.c.'s to matrix diagonal, and
        // zero out offdiag values at domain boundaries

        const Array< Array<BoundCond> > & bcs_i = bd.bndryConds(i);
        const BndryData::RealTuple      & bcl_i = bd.bndryLocs(i);
        
        const Box& domain = geom.Domain();
        for (OrientationIter oitr; oitr; oitr++) {
            int cdir(oitr());
            int idim = oitr().coordDir();
            const int bctype = bcs_i[cdir][0];
            const Real &bcl  = bcl_i[cdir];
            const Mask &msk  = bd.bndryMasks(oitr())[mfi];

            // Treat an exposed grid edge here as a boundary condition
            // for the linear solver:
            
            if (reg[oitr()] == domain[oitr()]) {
                const FArrayBox &fs = bd.bndryValues(oitr())[mfi];
                hmmat3(mat, ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                       cdir, bctype, bho, bcl,
                       ARLIM(fs.loVect()), ARLIM(fs.hiVect()),
                       BL_TO_FORTRAN(msk),
                       BL_TO_FORTRAN(bcoefs[idim][mfi]),
                       scalar_b, dx);
            }
            else {
                hmmat(mat, ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
                      cdir, bctype, bho, bcl,
                      BL_TO_FORTRAN(msk),
                      BL_TO_FORTRAN(bcoefs[idim][mfi]),
                      scalar_b, dx);
            }
            
        }

        // initialize matrix
        HYPRE_SStructMatrixSetBoxValues(A, part,
                                        const_cast<int*>(reg.loVect()), 
                                        const_cast<int*>(reg.hiVect()),
                                        0, size, stencil_indices, mat);

        hypre_TFree(mat);
    }
}

}
