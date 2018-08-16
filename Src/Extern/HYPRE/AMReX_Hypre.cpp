
#include <AMReX_Hypre.H>
#include <AMReX_HypreABecLap.H>
#include <AMReX_HypreABecLap2.H>
#include <AMReX_HypreABecLap3.H>
#include <HYPRE_config.h>
#include <type_traits>

namespace amrex {

constexpr HYPRE_Int Hypre::regular_stencil_size;
constexpr HYPRE_Int Hypre::eb_stencil_size;

std::unique_ptr<Hypre>
makeHypre (const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, MPI_Comm comm_, Hypre::Interface interface)
{
    if (interface == Hypre::Interface::structed) {
        return std::unique_ptr<Hypre>(new HypreABecLap(grids, dmap, geom, comm_));
    } else if (interface == Hypre::Interface::semi_structed) {
        return std::unique_ptr<Hypre>(new HypreABecLap2(grids, dmap, geom, comm_));
    } else {
        return std::unique_ptr<Hypre>(new HypreABecLap3(grids, dmap, geom, comm_));
    }    
}

Hypre::Hypre (const BoxArray& grids, const DistributionMapping& dmap,
              const Geometry& geom_, MPI_Comm comm_)
    : comm(comm_),
      geom(geom_)
{
    static_assert(AMREX_SPACEDIM > 1, "Hypre: 1D not supported");

    static_assert(std::is_same<Real, HYPRE_Real>::value, "amrex::Real != HYPRE_Real");
#ifdef HYPRE_BIGINT
    static_assert(std::is_same<long long int, HYPRE_Int>::value, "long long int != HYPRE_Int");
#else
    static_assert(std::is_same<int, HYPRE_Int>::value, "int != HYPRE_Int");
#endif

    const int ncomp = 1;
    int ngrow = 0;
    acoefs.define(grids, dmap, ncomp, ngrow);
    acoefs.setVal(0.0);
    
#ifdef AMREX_USE_EB
    ngrow = 1;
#endif
    
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        BoxArray edge_boxes(grids);
        edge_boxes.surroundingNodes(i);
        bcoefs[i].define(edge_boxes, dmap, ncomp, ngrow);
        bcoefs[i].setVal(0.0);
    }

    diaginv.define(grids,dmap,ncomp,0);
}

Hypre::~Hypre ()
{
    m_factory = nullptr;
    m_bndry = nullptr;
    m_maxorder = -1;
}

void
Hypre::setScalars (Real sa, Real sb)
{
    scalar_a = sa;
    scalar_b = sb;
}

void
Hypre::setACoeffs (const MultiFab& alpha)
{
    MultiFab::Copy(acoefs, alpha, 0, 0, 1, 0);
}

void
Hypre::setBCoeffs (const Array<const MultiFab*, BL_SPACEDIM>& beta)
{
    for (int idim=0; idim < AMREX_SPACEDIM; idim++) {
        const int ng = std::min(bcoefs[idim].nGrow(), beta[idim]->nGrow());
        MultiFab::Copy(bcoefs[idim], *beta[idim], 0, 0, 1, ng);
    }
}

void
Hypre::setVerbose (int _verbose)
{
    verbose = _verbose;
}

}  // namespace amrex
