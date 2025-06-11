#ifndef AMREX_ML_NODE_LINOP_K_H_
#define AMREX_ML_NODE_LINOP_K_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>

namespace amrex::nodelap_detail {

#ifdef AMREX_USE_HYPRE

    struct GetNode {
        AMREX_GPU_DEVICE Dim3 operator() (Dim3 const& lo, Dim3 const& len, int& offset)
        {
            Dim3 node;
            constexpr int nsten = AMREX_D_TERM(3,*3,*3);
            int icell = offset / nsten;
            node.z =  icell /        (len.x*len.y);
            node.y = (icell - node.z*(len.x*len.y)) /        len.x;
            node.x = (icell - node.z*(len.x*len.y)) - node.y*len.x;
            node.x += lo.x;
            node.y += lo.y;
            node.z += lo.z;
            offset -= icell*nsten;
            return node;
        }
    };

    struct GetNode2 {
        AMREX_GPU_DEVICE Dim3 operator() (int offset, Dim3 const& node)
        {
            // In 2D the offsets are
            //   6 7 8
            //   4 0 5
            //   1 2 3
            constexpr int nstenhalf = AMREX_SPACEDIM == 2 ? 4 : 13;
            if (offset == 0) {
                return node;
            } else {
                if (offset <= nstenhalf) { --offset; }
                Dim3 node2;
                node2.z = offset / 9;
                node2.y = (offset - node2.z*9) / 3;
                node2.x = (offset - node2.z*9) - node2.y*3;
                AMREX_D_TERM(node2.x += node.x-1;,
                             node2.y += node.y-1;,
                             node2.z += node.z-1);
                return node2;
            }
        }
    };

#endif /* AMREX_USE_HYPRE */

    constexpr int crse_cell = 0; // Do NOT change the values
    constexpr int fine_cell = 1;
    constexpr int crse_node = 0;
    constexpr int crse_fine_node = 1;
    constexpr int fine_node = 2;
}

namespace amrex {

template <typename T>
void mlndlap_fillbc_cc (Box const& vbx, Array4<T> const& sigma, Box const& domain,
                        GpuArray<LinOpBCType, AMREX_SPACEDIM> bclo,
                        GpuArray<LinOpBCType, AMREX_SPACEDIM> bchi) noexcept
{
    GpuArray<bool,AMREX_SPACEDIM> bflo{{AMREX_D_DECL(bclo[0] != LinOpBCType::Periodic,
                                                     bclo[1] != LinOpBCType::Periodic,
                                                     bclo[2] != LinOpBCType::Periodic)}};
    GpuArray<bool,AMREX_SPACEDIM> bfhi{{AMREX_D_DECL(bchi[0] != LinOpBCType::Periodic,
                                                     bchi[1] != LinOpBCType::Periodic,
                                                     bchi[2] != LinOpBCType::Periodic)}};
    mlndlap_bc_doit(vbx, sigma, domain, bflo, bfhi);
}

template <typename T>
void mlndlap_applybc (Box const& vbx, Array4<T> const& phi, Box const& domain,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> bclo,
                      GpuArray<LinOpBCType, AMREX_SPACEDIM> bchi) noexcept
{
    GpuArray<bool,AMREX_SPACEDIM> bflo{{AMREX_D_DECL(bclo[0] == LinOpBCType::Neumann ||
                                                     bclo[0] == LinOpBCType::inflow,
                                                     bclo[1] == LinOpBCType::Neumann ||
                                                     bclo[1] == LinOpBCType::inflow,
                                                     bclo[2] == LinOpBCType::Neumann ||
                                                     bclo[2] == LinOpBCType::inflow)}};
    GpuArray<bool,AMREX_SPACEDIM> bfhi{{AMREX_D_DECL(bchi[0] == LinOpBCType::Neumann ||
                                                     bchi[0] == LinOpBCType::inflow,
                                                     bchi[1] == LinOpBCType::Neumann ||
                                                     bchi[1] == LinOpBCType::inflow,
                                                     bchi[2] == LinOpBCType::Neumann ||
                                                     bchi[2] == LinOpBCType::inflow)}};
    mlndlap_bc_doit(vbx, phi, domain, bflo, bfhi);
}

}

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLNodeLinOp_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLNodeLinOp_2D_K.H>
#else
#include <AMReX_MLNodeLinOp_3D_K.H>
#endif

#endif
