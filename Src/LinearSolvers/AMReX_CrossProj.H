#ifndef AMREX_CROSS_PROJ_H_
#define AMREX_CROSS_PROJ_H_
#include <AMReX_Config.H>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>

#if (AMREX_SPACEDIM == 2)

namespace amrex::GMG
{

/**
 * \file AMReX_CrossProj.H
 *
 * Cross-stencil projection for periodic domains, solved with geometric
 * multigrid.
 */

/**
 * \brief Projection of a cell-centered velocity field onto the space of
 * fields that are divergence free with respect to the node-centered
 * cross-stencil divergence.
 *
 * Same projection and interface as FFT::CrossProj, in 2D only. The
 * divergence at each node is computed from the 4 cells surrounding it.
 * Currently only all-periodic boundaries are supported, with an even number
 * (at least 4) of cells in each direction starting at index 0, and even box
 * lengths. The first two components of the MFs are the velocity.
 */
class CrossProj
{
public:

    //! \param geom Periodic geometry (all directions must be periodic).
    explicit CrossProj (Geometry const& geom);

    ~CrossProj ();
    CrossProj (CrossProj const&) = delete;
    CrossProj (CrossProj&&) = delete;
    CrossProj& operator= (CrossProj const&) = delete;
    CrossProj& operator= (CrossProj&&) = delete;

    /**
     * \brief Project the velocity in \p a_rhs and store the result in \p a_soln.
     *
     * The two MFs must have the same BoxArray and DistributionMapping, and
     * they could be the same MF. If \p a_soln has ghost cells, one layer of
     * ghost cells will be filled.
     *
     * \param a_soln Projected velocity.
     * \param a_rhs  Velocity to be projected.
     */
    void solve (MultiFab& a_soln, MultiFab const& a_rhs);

    //! Relative tolerance (default 1e-10, 1e-5 in single precision).
    void setRelTol (Real rtol) { m_reltol = rtol; }
    //! Absolute tolerance (default 0).
    void setAbsTol (Real atol) { m_abstol = atol; }
    void setMaxIter (int n) { m_maxiter = n; }
    void setVerbose (int v) { m_verbose = v; }
    //! Number of iterations of the last solve.
    [[nodiscard]] int getNumIters () const { return m_niters; }

    struct Level;

private:

    void setup (BoxArray const& ba, DistributionMapping const& dm);
    void precond ();
    void vcycle (int h, int lev);

    Geometry m_geom;
    Real m_reltol = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
    Real m_abstol = Real(0.0);
    int m_maxiter = 200;
    int m_verbose = 0;
    int m_niters = 0;

    BoxArray m_grids; // cell-centered
    DistributionMapping m_dmap;
    bool m_reduced = false; // dx == dy: PCG on the red/black reduced problem
    MultiFab m_vel, m_x, m_p, m_q, m_rold;
    std::unique_ptr<Level> m_top;
    // Coarse hierarchies: 0 for smooth modes, 1 for the checkerboard.
    Array<Vector<std::unique_ptr<Level>>,2> m_hier;
};

}

#endif
#endif
