#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLNodeLap_K.H>

#include <limits>

namespace amrex {

#if defined(AMREX_USE_HYPRE) && (AMREX_SPACEDIM > 1)

void
MLNodeLaplacian::fillIJMatrix (MFIter const& mfi,
                               Array4<HypreNodeLap::AtomicInt const> const& gid,
                               Array4<int const> const& lid,
                               HypreNodeLap::Int* const ncols,
                               HypreNodeLap::Int* const cols,
                               Real* const mat) const
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        fillIJMatrix_gpu(mfi,gid,lid,ncols,cols,mat);
    } else
#endif
    {
        fillIJMatrix_cpu(mfi,gid,lid,ncols,cols,mat);
    }
}

#ifdef AMREX_USE_GPU

void
MLNodeLaplacian::fillIJMatrix_gpu (MFIter const& mfi,
                                   Array4<HypreNodeLap::AtomicInt const> const& gid,
                                   Array4<int const> const& lid,
                                   HypreNodeLap::Int* const ncols,
                                   HypreNodeLap::Int* const cols,
                                   Real* const mat) const
{
    const int amrlev = 0;
    const int mglev  = m_num_mg_levels[amrlev]-1;

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const Box& ndbx = mfi.validbox();
    const auto ndlo = amrex::lbound(ndbx);
    const auto ndlen = amrex::length(ndbx);
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE
        (static_cast<Long>(ndbx.numPts())*AMREX_D_TERM(3,*3,*3) <
         static_cast<Long>(std::numeric_limits<int>::max()),
         "The Box is too big.  We could use Long here, but it would much slower.");
    const int nmax = ndbx.numPts() * AMREX_D_TERM(3,*3,*3);

    int nelems;

    if (m_coarsening_strategy == CoarseningStrategy::RAP)
    {
        const auto& sten = stencil->const_array(mfi);
        nelems = amrex::Scan::PrefixSum<int>
            (nmax,
             [=] AMREX_GPU_DEVICE (int offset) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 Dim3 node2 = GetNode2()(offset, node);
                 return (lid(node.x,node.y,node.z) >= 0 &&
                         gid(node2.x,node2.y,node2.z)
                         < std::numeric_limits<HypreNodeLap::AtomicInt>::max());
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 mlndlap_fillijmat_sten_gpu(ps, node.x, node.y, node.z, offset, gid, lid,
                                            ncols, cols, mat, sten);
             },
             amrex::Scan::Type::exclusive);
    }
    else if (sigma[0] == nullptr) // const sigma
    {
        Real const_sigma = m_const_sigma;
        nelems = amrex::Scan::PrefixSum<int>
            (nmax,
             [=] AMREX_GPU_DEVICE (int offset) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 Dim3 node2 = GetNode2()(offset, node);
                 return (lid(node.x,node.y,node.z) >= 0 &&
                         gid(node2.x,node2.y,node2.z)
                         < std::numeric_limits<HypreNodeLap::AtomicInt>::max());
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 mlndlap_fillijmat_cs_gpu(ps, node.x, node.y, node.z, offset,
                                          ndbx, gid, lid, ncols, cols, mat,
                                          const_sigma, dxinvarr, domain
#if (AMREX_SPACEDIM == 2)
                                          , is_rz
#endif
                     );
             },
             amrex::Scan::Type::exclusive);
    }
    else if ( (m_use_harmonic_average && mglev > 0) ||
              (m_use_mapped) )
    {
        AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                     Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                     Array4<Real const> const& szarr = sigma[2]->const_array(mfi));
        nelems = amrex::Scan::PrefixSum<int>
            (nmax,
             [=] AMREX_GPU_DEVICE (int offset) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 Dim3 node2 = GetNode2()(offset, node);
                 return (lid(node.x,node.y,node.z) >= 0 &&
                         gid(node2.x,node2.y,node2.z)
                         < std::numeric_limits<HypreNodeLap::AtomicInt>::max());
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 mlndlap_fillijmat_ha_gpu(ps, node.x, node.y, node.z, offset,
                                          ndbx, gid, lid, ncols, cols, mat,
                                          AMREX_D_DECL(sxarr, syarr, szarr),
                                          dxinvarr, domain
#if (AMREX_SPACEDIM == 2)
                                          , is_rz
#endif
                     );
             },
             amrex::Scan::Type::exclusive);
    }
    else
    {
        Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
        nelems = amrex::Scan::PrefixSum<int>
            (nmax,
             [=] AMREX_GPU_DEVICE (int offset) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 Dim3 node2 = GetNode2()(offset, node);
                 return (lid(node.x,node.y,node.z) >= 0 &&
                         gid(node2.x,node2.y,node2.z)
                         < std::numeric_limits<HypreNodeLap::AtomicInt>::max());
             },
             [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
             {
                 Dim3 node = GetNode()(ndlo, ndlen, offset);
                 mlndlap_fillijmat_aa_gpu(ps, node.x, node.y, node.z, offset,
                                          ndbx, gid, lid, ncols, cols, mat,
                                          sarr, dxinvarr, domain
#if (AMREX_SPACEDIM == 2)
                                          , is_rz
#endif
                     );
             },
             amrex::Scan::Type::exclusive);
    }

    amrex::ignore_unused(nelems);
}

#endif

void
MLNodeLaplacian::fillIJMatrix_cpu (MFIter const& mfi,
                                   Array4<HypreNodeLap::AtomicInt const> const& gid,
                                   Array4<int const> const& lid,
                                   HypreNodeLap::Int* const ncols,
                                   HypreNodeLap::Int* const cols,
                                   Real* const mat) const
{
    const int amrlev = 0;
    const int mglev  = m_num_mg_levels[amrlev]-1;

    const auto& sigma = m_sigma[amrlev][mglev];
    const auto& stencil = m_stencil[amrlev][mglev];
    const auto dxinvarr = m_geom[amrlev][mglev].InvCellSizeArray();
#if (AMREX_SPACEDIM == 2)
    bool is_rz = m_is_rz;
#endif

    const Box& ndbx = mfi.validbox();
    const Box& domain = m_geom[amrlev][mglev].growPeriodicDomain(1);

    if (m_coarsening_strategy == CoarseningStrategy::RAP)
    {
        const auto& sten = stencil->const_array(mfi);
        mlndlap_fillijmat_sten_cpu(ndbx, gid, lid, ncols, cols, mat, sten);
    }
    else if (sigma[0] == nullptr) // const sigma
    {
        Real const_sigma = m_const_sigma;
        mlndlap_fillijmat_cs_cpu(ndbx, gid, lid, ncols, cols, mat,
                                 const_sigma, dxinvarr, domain
#if (AMREX_SPACEDIM == 2)
                                 , is_rz
#endif
            );
    }
    else if ( (m_use_harmonic_average && mglev > 0) ||
              (m_use_mapped) )
    {
        AMREX_D_TERM(Array4<Real const> const& sxarr = sigma[0]->const_array(mfi);,
                     Array4<Real const> const& syarr = sigma[1]->const_array(mfi);,
                     Array4<Real const> const& szarr = sigma[2]->const_array(mfi));
        mlndlap_fillijmat_ha_cpu(ndbx, gid, lid, ncols, cols, mat,
                                 AMREX_D_DECL(sxarr,syarr,szarr), dxinvarr, domain
#if (AMREX_SPACEDIM == 2)
                                 , is_rz
#endif
            );
    }
    else
    {
        Array4<Real const> const& sarr = sigma[0]->const_array(mfi);
        mlndlap_fillijmat_aa_cpu(ndbx, gid, lid, ncols, cols, mat,
                                 sarr, dxinvarr, domain
#if (AMREX_SPACEDIM == 2)
                                 , is_rz
#endif
            );
    }
}

void
MLNodeLaplacian::fillRHS (MFIter const& mfi, Array4<int const> const& lid,
                          Real* const rhs, Array4<Real const> const& bfab) const
{
    const int amrlev = 0;
    const int mglev  = m_num_mg_levels[amrlev]-1;
    const Box& nddom = amrex::surroundingNodes(Geom(amrlev,mglev).Domain());
    const Box& bx = mfi.validbox();
    const auto lobc = LoBC();
    const auto hibc = HiBC();
    GpuArray<int,AMREX_SPACEDIM> neumann_lo{AMREX_D_DECL(std::numeric_limits<int>::lowest(),
                                                         std::numeric_limits<int>::lowest(),
                                                         std::numeric_limits<int>::lowest())};
    GpuArray<int,AMREX_SPACEDIM> neumann_hi{AMREX_D_DECL(std::numeric_limits<int>::lowest(),
                                                         std::numeric_limits<int>::lowest(),
                                                         std::numeric_limits<int>::lowest())};
    if (m_coarsening_strategy == CoarseningStrategy::Sigma) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (lobc[idim] == LinOpBCType::Neumann || lobc[idim] == LinOpBCType::inflow) {
                if (bx.smallEnd(idim) == nddom.smallEnd(idim)) {
                    neumann_lo[idim] = bx.smallEnd(idim);
                }
            }
            if (hibc[idim] == LinOpBCType::Neumann || hibc[idim] == LinOpBCType::inflow) {
                if (bx.bigEnd(idim) == nddom.bigEnd(idim)) {
                    neumann_hi[idim] = bx.bigEnd(idim);
                }
            }
        }
    }

    AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
    {
        if (lid(i,j,k) >= 0) {
            Real fac = Real(1.0);
            AMREX_D_TERM(if ((neumann_lo[0] == i) || neumann_hi[0] == i) { fac *= 0.5; },
                         if ((neumann_lo[1] == j) || neumann_hi[1] == j) { fac *= 0.5; },
                         if ((neumann_lo[2] == k) || neumann_hi[2] == k) { fac *= 0.5; })
            rhs[lid(i,j,k)] = fac * bfab(i,j,k);
        }
    });
}

#endif

}
