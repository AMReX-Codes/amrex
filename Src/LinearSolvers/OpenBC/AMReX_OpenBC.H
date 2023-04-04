#ifndef AMREX_OPENBC_H_
#define AMREX_OPENBC_H_
#include <AMReX_Config.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

namespace amrex
{

namespace openbc {

    static constexpr int M = 7; // highest order of moments
    static constexpr int P = 3;

    struct Moments
    {
        typedef GpuArray<Real,(M+2)*(M+1)/2> array_type;
        array_type mom;
        Real x, y, z;
        Orientation face;
    };

    struct MomTag
    {
        Array4<Real const> gp;
        Box b2d;
        Orientation face;
        int offset;
    };

    std::ostream& operator<< (std::ostream& os, Moments const& mom);
}

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
template<>
struct Gpu::SharedMemory<openbc::Moments::array_type>
{
    AMREX_GPU_DEVICE openbc::Moments::array_type* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(openbc::Moments::array_type,amrex_openbc_momarray);,
                          extern __shared__  openbc::Moments::array_type amrex_openbc_momarray[];)
            return amrex_openbc_momarray;
    }
};
#endif

/**
 * \brief Open Boundary Poisson Solver
 *
 * References:
 *    (1) The Solution of Poisson's Equation for Isolated Source
 *        Distributions, R. A. James, 1977, JCP 25, 71
 *    (2) A Local Corrections Algorithm for Solving Poisson's Equation in Three
 *        Dimensions, P. McCorquodale, P. Colella, G. T. Balls, & S. B. Baden,
 *        2007, Communications in Applied Mathematics and Computational Science,
 *        2, 1, 57-81
 */
class OpenBCSolver
{
public:
    OpenBCSolver () = default;

    OpenBCSolver (const Vector<Geometry>& a_geom,
                  const Vector<BoxArray>& a_grids,
                  const Vector<DistributionMapping>& a_dmap,
                  const LPInfo& a_info = LPInfo());

    ~OpenBCSolver () = default;

    OpenBCSolver (const OpenBCSolver&) = delete;
    OpenBCSolver (OpenBCSolver&&) = delete;
    OpenBCSolver& operator= (const OpenBCSolver&) = delete;
    OpenBCSolver& operator= (OpenBCSolver&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo());

    void setVerbose (int v) noexcept;
    void setBottomVerbose (int v) noexcept;

    void useHypre (bool use_hypre) noexcept;

    Real solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
                Real a_tol_rel, Real a_tol_abs);

// public for cuda

    void compute_moments (Gpu::DeviceVector<openbc::Moments>& moments);
    void compute_potential (Gpu::DeviceVector<openbc::Moments> const& moments);
    void interpolate_potential (MultiFab& solg);

private:

#ifdef AMREX_USE_MPI
    void bcast_moments (Gpu::DeviceVector<openbc::Moments>& moments);
#endif

    int m_verbose = 0;
    int m_bottom_verbose = 0;
    Vector<Geometry> m_geom;
    Vector<BoxArray> m_grids;
    Vector<DistributionMapping> m_dmap;
    LPInfo m_info;
    std::unique_ptr<MLPoisson> m_poisson_1;
    std::unique_ptr<MLPoisson> m_poisson_2;
    std::unique_ptr<MLMG> m_mlmg_1;
    std::unique_ptr<MLMG> m_mlmg_2;
    BottomSolver m_bottom_solver_type = BottomSolver::bicgstab;

    int m_coarsen_ratio = 0;
    Array<MultiFab,AMREX_SPACEDIM> m_dpdn;
    Gpu::PinnedVector<openbc::MomTag> m_momtags_h;
#ifdef AMREX_USE_GPU
    Gpu::DeviceVector<openbc::MomTag> m_momtags_d;
    Gpu::PinnedVector<int> m_ngpublocks_h;
    Gpu::DeviceVector<int> m_ngpublocks_d;
    int m_nthreads_momtag;
#endif

    int m_nblocks_local = 0;
    int m_nblocks = 0;
#ifdef AMREX_USE_MPI
    Vector<int> m_countvec;
    Vector<int> m_offset;
#endif

    IntVect m_ngrowdomain;
    MultiFab m_crse_grown_faces_phi;
    MultiFab m_phind;
    BoxArray m_bag;

    Vector<IntVect> m_box_offset;
    Vector<BoxArray> m_ba_all;
    Vector<DistributionMapping> m_dm_all;
    Vector<Geometry> m_geom_all;
};

}

#endif
