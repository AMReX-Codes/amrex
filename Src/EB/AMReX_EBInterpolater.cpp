
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_Geometry.H>

namespace amrex {

EBCellConservativeLinear  eb_lincc_interp;
EBCellConservativeLinear  eb_cell_cons_interp(0);

EBCellConservativeLinear::EBCellConservativeLinear (bool do_linear_limiting_)
    : CellConservativeLinear(do_linear_limiting_)
{
}

EBCellConservativeLinear::~EBCellConservativeLinear ()
{
}

void
EBCellConservativeLinear::interp (const FArrayBox& crse,
                                  int              crse_comp,
                                  FArrayBox&       fine,
                                  int              fine_comp,
                                  int              ncomp,
                                  const Box&       fine_region,
                                  const IntVect&   ratio,
                                  const Geometry&  crse_geom,
                                  const Geometry&  fine_geom,
                                  Vector<BCRec> const&  bcr,
                                  int              actual_comp,
                                  int              actual_state,
                                  RunOn            runon)
{
    CellConservativeLinear::interp(crse, crse_comp, fine, fine_comp, ncomp, fine_region, ratio,
                                   crse_geom, fine_geom, bcr, actual_comp, actual_state, runon);

    const Box& target_fine_region = fine_region & fine.box();

    if (crse.getType() == FabType::regular)
    {
        BL_ASSERT(amrex::getEBCellFlagFab(fine).getType(target_fine_region) == FabType::regular);
    }
    else
    {
        const EBFArrayBox& crse_eb = static_cast<EBFArrayBox const&>(crse);
        EBFArrayBox&       fine_eb = static_cast<EBFArrayBox      &>(fine);
        
        const EBCellFlagFab& crse_flag = crse_eb.getEBCellFlagFab();
        const EBCellFlagFab& fine_flag = fine_eb.getEBCellFlagFab();
        
        const Box& crse_bx = CoarseBox(target_fine_region,ratio);
    
        const FabType ftype = fine_flag.getType(target_fine_region);
        const FabType ctype = crse_flag.getType(crse_bx);

        if (ftype == FabType::multivalued || ctype == FabType::multivalued)
        {
            amrex::Abort("EBCellConservativeLinear::interp: multivalued not implemented");
        }
        else if (ftype == FabType::covered)
        {
            ; // don't need to do anything special
        }
        else
        {
            const auto& cflag = crse_flag.const_array();
            auto const& fa = fine.array();
            auto const& ca = crse.const_array();
            bool run_on_gpu = (Gpu::inLaunchRegion() && runon == RunOn::Gpu);
            AMREX_HOST_DEVICE_FOR_4D_FLAG(run_on_gpu, target_fine_region, ncomp, i, j, k, n,
            {
                Dim3 cxyz = amrex::coarsen(Dim3{i,j,k}, ratio);
                if (cflag(cxyz.x,cxyz.y,cxyz.z).numNeighbors() < AMREX_D_TERM(3,*3,*3)) {
                    fa(i,j,k,n+fine_comp) = ca(cxyz.x,cxyz.y,cxyz.z,n+crse_comp);
                }
            });
        }
    }        
}

}
