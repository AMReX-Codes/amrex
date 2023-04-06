
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Interpolater.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>

#include <climits>

namespace amrex {

/*
 * PCInterp, NodeBilinear, FaceLinear, CellConservativeLinear, and
 * CellBilinear are supported for all dimensions on cpu and gpu.
 *
 * CellConservativeProtected only works in 2D and 3D on cpu and gpu.
 *
 * CellQuadratic only works in 2D and 3D on cpu and gpu.
 *
 * CellQuartic works in 1D, 2D and 3D on cpu and gpu with ref ratio of 2
 *
 * CellConservativeQuartic only works with ref ratio of 2 on cpu and gpu.
 *
 * FaceDivFree works in 2D and 3D on cpu and gpu.
 * The algorithm is restricted to ref ratio of 2.
 */

//
// CONSTRUCT A GLOBAL OBJECT OF EACH VERSION.
//
PCInterp                  pc_interp;
NodeBilinear              node_bilinear_interp;
FaceLinear                face_linear_interp;
FaceDivFree               face_divfree_interp;
CellConservativeLinear    lincc_interp;
CellConservativeLinear    cell_cons_interp(false);
CellConservativeProtected protected_interp;
CellConservativeQuartic   quartic_interp;
CellBilinear              cell_bilinear_interp;
CellQuadratic             quadratic_interp;
CellQuartic               cell_quartic_interp;

Box
NodeBilinear::CoarseBox (const Box& fine,
                         int        ratio)
{
    Box b = amrex::coarsen(fine,ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

Box
NodeBilinear::CoarseBox (const Box&     fine,
                         const IntVect& ratio)
{
    Box b = amrex::coarsen(fine,ratio);

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (b.length(i) < 2)
        {
            //
            // Don't want degenerate boxes.
            //
            b.growHi(i,1);
        }
    }

    return b;
}

void
NodeBilinear::interp (const FArrayBox&  crse,
                      int               crse_comp,
                      FArrayBox&        fine,
                      int               fine_comp,
                      int               ncomp,
                      const Box&        fine_region,
                      const IntVect&    ratio,
                      const Geometry& /*crse_geom */,
                      const Geometry& /*fine_geom */,
                      Vector<BCRec> const& /*bcr*/,
                      int               /*actual_comp*/,
                      int               /*actual_state*/,
                      RunOn             runon)
{
    BL_PROFILE("NodeBilinear::interp()");

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
    {
        mf_nodebilin_interp(i,j,k,n, finearr, fine_comp, crsearr, crse_comp, ratio);
    });
}

Box
FaceLinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine, IntVect(ratio));
}

Box
FaceLinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box b = amrex::coarsen(fine,ratio);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (b.type(i) == IndexType::NODE && b.length(i) < 2) {
            // Don't want degenerate boxes in nodal direction.
            b.growHi(i,1);
        }
    }
    return b;
}

void
FaceLinear::interp (const FArrayBox&  crse,
                    int               crse_comp,
                    FArrayBox&        fine,
                    int               fine_comp,
                    int               ncomp,
                    const Box&        fine_region,
                    const IntVect&    ratio,
                    const Geometry&   crse_geom ,
                    const Geometry&   fine_geom ,
                    Vector<BCRec> const& bcr,
                    int               actual_comp,
                    int               /*actual_state*/,
                    RunOn             runon)
{
    //
    // This version is called from InterpFromCoarseLevel
    //
    BL_PROFILE("FaceLinear::interp()");

    // pass unallocated IArrayBox for solve_mask, so all fine values get filled.
    interp_face(crse, crse_comp, fine, fine_comp, ncomp, fine_region,
                ratio, IArrayBox(), crse_geom, fine_geom, bcr, actual_comp, runon);
}

void
FaceLinear::interp_face (const FArrayBox&  crse,
                         const int         crse_comp,
                         FArrayBox&        fine,
                         const int         fine_comp,
                         const int         ncomp,
                         const Box&        fine_region,
                         const IntVect&    ratio,
                         const IArrayBox&  solve_mask,
                         const Geometry& /*crse_geom */,
                         const Geometry& /*fine_geom */,
                         Vector<BCRec> const& /*bcr*/,
                         const int         /*bccomp*/,
                         RunOn             runon)
{
    BL_PROFILE("FaceLinear::interp_face()");

    AMREX_ASSERT(AMREX_D_TERM(fine_region.type(0),+fine_region.type(1),+fine_region.type(2)) == 1);

    Array4<Real> const& fine_arr = fine.array(fine_comp);
    Array4<Real const> const& crse_arr = crse.const_array(crse_comp);
    Array4<const int> mask_arr;
    if (solve_mask.isAllocated()) {
        mask_arr = solve_mask.const_array();
    }

    //
    // Fill fine ghost faces with piecewise-constant interpolation of coarse data.
    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
    // coarse face, leave the in-between faces alone.
    // The mask ensures we do not overwrite valid fine cells.
    //
    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_face_interp_x(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_face_interp_y(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_face_interp_z(i,j,k,n,fine_arr,crse_arr,mask_arr,ratio);
        });
    }
#endif
#endif

    //
    // Interpolate unfilled grow cells using best data from
    // surrounding faces of valid region, and pc-interpd data
    // on fine faces overlaying coarse edges.
    //
    if (fine_region.type(0) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_x(i,j,k,n,fine_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM >= 2)
    else if (fine_region.type(1) == IndexType::NODE)
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_y(i,j,k,n,fine_arr,ratio);
        });
    }
#if (AMREX_SPACEDIM == 3)
    else
    {
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
        {
            face_linear_interp_z(i,j,k,n,fine_arr,ratio);
        });
    }
#endif
#endif
}

void FaceLinear::interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                             const int         crse_comp,
                             Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                             const int         fine_comp,
                             const int         ncomp,
                             const Box&        fine_region,
                             const IntVect&    ratio,
                             Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                             const Geometry&   /*crse_geom*/,
                             const Geometry&   /*fine_geom*/,
                             Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                             const int         /*actual_comp*/,
                             const int         /*actual_state*/,
                             const RunOn       runon)
{
    BL_PROFILE("FaceLinear::interp_arr()");

    Array<IndexType, AMREX_SPACEDIM> types;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        { types[d].set(d); }

    GpuArray<Array4<const Real>, AMREX_SPACEDIM> crse_arr;
    GpuArray<Array4<Real>, AMREX_SPACEDIM> fine_arr;
    GpuArray<Array4<const int>, AMREX_SPACEDIM> mask_arr;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        crse_arr[d] = crse[d]->const_array(crse_comp);
        fine_arr[d] = fine[d]->array(fine_comp);
        if (solve_mask[d] != nullptr)
            { mask_arr[d] = solve_mask[d]->const_array(0); }
    }

    //
    // Fill fine ghost faces with piecewise-constant interpolation of coarse data.
    // Operate only on faces that overlap--ie, only fill the fine faces that make up each
    // coarse face, leave the in-between faces alone.
    // The mask ensures we do not overwrite valid fine cells.
    //
    // Fuse the launches, 1 for each dimension, into a single launch.
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_face_interp_x(i,j,k,n,fine_arr[0],crse_arr[0],mask_arr[0],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_face_interp_y(i,j,k,n,fine_arr[1],crse_arr[1],mask_arr[1],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_face_interp_z(i,j,k,n,fine_arr[2],crse_arr[2],mask_arr[2],ratio);
                      }
                  });
              });

    //
    // Interpolate unfilled grow cells using best data from
    // surrounding faces of valid region, and pc-interpd data
    // on fine faces overlaying coarse edges.
    //
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_x(i,j,k,n,fine_arr[0],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_y(i,j,k,n,fine_arr[1],ratio);
                      }
                  });
              },
              amrex::convert(fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          face_linear_interp_z(i,j,k,n,fine_arr[2],ratio);
                      }
                  });
              });
}

Box
CellBilinear::CoarseBox (const Box& fine, int ratio)
{
    return CoarseBox(fine, IntVect(ratio));
}

Box
CellBilinear::CoarseBox (const Box& fine, const IntVect& ratio)
{
    const int* lo = fine.loVect();
    const int* hi = fine.hiVect();

    Box crse(amrex::coarsen(fine,ratio));
    const int* clo = crse.loVect();
    const int* chi = crse.hiVect();

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if ((lo[i]-clo[i]*ratio[i])*2 < ratio[i]) {
            crse.growLo(i,1);
        }
        if ((hi[i]-chi[i]*ratio[i])*2 >= ratio[i]) {
            crse.growHi(i,1);
        }
    }
    return crse;
}

void
CellBilinear::interp (const FArrayBox&  crsefab,
                      int               crse_comp,
                      FArrayBox&        finefab,
                      int               fine_comp,
                      int               ncomp,
                      const Box&        fine_region,
                      const IntVect &   ratio,
                      const Geometry& /*crse_geom*/,
                      const Geometry& /*fine_geom*/,
                      Vector<BCRec> const& /*bcr*/,
                      int               /*actual_comp*/,
                      int               /*actual_state*/,
                      RunOn             runon)
{
    BL_PROFILE("CellBilinear::interp()");

    auto const& crse = crsefab.const_array();
    auto const& fine = finefab.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,fine_region,ncomp,i,j,k,n,
    {
        mf_cell_bilin_interp(i,j,k,n, fine, fine_comp, crse, crse_comp, ratio);
    });
}


CellConservativeLinear::CellConservativeLinear (bool do_linear_limiting_)
    : do_linear_limiting(do_linear_limiting_)
{}

Box
CellConservativeLinear::CoarseBox (const Box&     fine,
                                   const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellConservativeLinear::CoarseBox (const Box& fine,
                                   int        ratio)
{
    Box crse(amrex::coarsen(fine,ratio));
    crse.grow(1);
    return crse;
}

void
CellConservativeLinear::interp (const FArrayBox& crse,
                                int              crse_comp,
                                FArrayBox&       fine,
                                int              fine_comp,
                                int              ncomp,
                                const Box&       fine_region,
                                const IntVect&   ratio,
                                const Geometry&  crse_geom,
                                const Geometry&  fine_geom,
                                Vector<BCRec> const& bcr,
                                int              /*actual_comp*/,
                                int              /*actual_state*/,
                                RunOn             runon)
{
    BL_PROFILE("CellConservativeLinear::interp()");
    BL_ASSERT(bcr.size() >= ncomp);

    AMREX_ASSERT(fine.box().contains(fine_region));

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());

    Box const& cdomain = crse_geom.Domain();
    amrex::ignore_unused(fine_geom);

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();

    const Box& crse_region = CoarseBox(fine_region,ratio);
    const Box& cslope_bx = amrex::grow(crse_region,-1);

    AsyncArray<BCRec> async_bcr(bcr.data(), (run_on_gpu) ? ncomp : 0);
    BCRec const* bcrp = (run_on_gpu) ? async_bcr.data() : bcr.data();

    FArrayBox ccfab(cslope_bx, ncomp*AMREX_SPACEDIM);
    Elixir cceli;
    if (run_on_gpu) cceli = ccfab.elixir();
    Array4<Real> const& tmp = ccfab.array();
    Array4<Real const> const& ctmp = ccfab.const_array();

#if (AMREX_SPACEDIM == 1)
    if (crse_geom.IsSPHERICAL()) {
        Real drf = fine_geom.CellSize(0);
        Real rlo = fine_geom.Offset(0);
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                amrex::ignore_unused(j,k);
                mf_cell_cons_lin_interp_mcslope_sph(i, n, tmp, crsearr, crse_comp, ncomp,
                                                    cdomain, ratio, bcrp, drf, rlo);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            amrex::ignore_unused(j,k);
            mf_cell_cons_lin_interp_sph(i, n, finearr, fine_comp, ctmp,
                                        crsearr, crse_comp, ncomp, ratio, drf, rlo);
        });
    } else
#elif (AMREX_SPACEDIM == 2)
    if (crse_geom.IsRZ()) {
        Real drf = fine_geom.CellSize(0);
        Real rlo = fine_geom.Offset(0);
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                amrex::ignore_unused(k);
                mf_cell_cons_lin_interp_mcslope_rz(i, j, n, tmp, crsearr, crse_comp, ncomp,
                                                   cdomain, ratio, bcrp, drf, rlo);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            amrex::ignore_unused(k);
            mf_cell_cons_lin_interp_rz(i, j, n, finearr, fine_comp, ctmp,
                                       crsearr, crse_comp, ncomp, ratio, drf, rlo);
        });
    } else
#endif
    {
        if (do_linear_limiting) {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cslope_bx, i, j, k,
            {
                mf_cell_cons_lin_interp_llslope(i,j,k, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, bcrp);
            });
        } else {
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
            {
                mf_cell_cons_lin_interp_mcslope(i,j,k,n, tmp, crsearr, crse_comp, ncomp,
                                                cdomain, ratio, bcrp);
            });
        }

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i, j, k, n,
        {
            mf_cell_cons_lin_interp(i,j,k,n, finearr, fine_comp, ctmp,
                                    crsearr, crse_comp, ncomp, ratio);
        });
    }
}

Box
CellQuadratic::CoarseBox (const Box&     fine,
                          const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

Box
CellQuadratic::CoarseBox (const Box& fine,
                          int        ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(1);
    return crse;
}

void
CellQuadratic::interp (const FArrayBox& crse,
                       int              crse_comp,
                       FArrayBox&       fine,
                       int              fine_comp,
                       int              ncomp,
                       const Box&       fine_region,
                       const IntVect&   ratio,
                       const Geometry&  crse_geom,
                       const Geometry&  fine_geom,
                       Vector<BCRec> const&  bcr,
                       int              /* actual_comp */,
                       int              /* actual_state */,
                       RunOn            runon)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(crse,crse_comp,fine,fine_comp,ncomp,fine_region,
                         ratio,crse_geom,fine_geom,bcr,runon);
    amrex::Abort("1D CellQuadratic::interp not supported");
#else

    BL_PROFILE("CellQuadratic::interp()");
    BL_ASSERT(bcr.size() >= ncomp);

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    // Make Box for slopes.
    Box cslope_bx = amrex::coarsen(target_fine_region,ratio);
    BL_ASSERT(crse.box().contains(cslope_bx));

    // Are we running on GPU?
    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());

    // Set up domain for coarse geometry
    Box const& cdomain = crse_geom.Domain();

    // Set up AsyncArray for boundary conditions
    AsyncArray<BCRec> async_bcr(bcr.data(), (run_on_gpu) ? ncomp : 0);
    BCRec const* bcrp = (run_on_gpu) ? async_bcr.data() : bcr.data();

    // Set up temporary fab (with elixir, as needed) for coarse grid slopes
#if (AMREX_SPACEDIM == 2)
    int nslp = 5; // x, y, x^2, y^2, xy, in that order.
#else  /* AMREX_SPACEDIM == 3 */
    int nslp = 9; // x, y, z, x^2, y^2, z^2, xy, xz, yz, in that order.
#endif /* AMREX_SPACEDIM == 2 */
    FArrayBox sfab(cslope_bx, nslp*ncomp);
    Elixir seli;
    if (run_on_gpu) seli = sfab.elixir();

    // Extract pointers to fab data
    Array4<Real>       const&   finearr = fine.array();
    Array4<Real const> const&   crsearr = crse.const_array();
    Array4<Real>       const&  slopearr = sfab.array();
    Array4<Real const> const& cslopearr = sfab.const_array();

    // Compute slopes.
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, cslope_bx, ncomp, i, j, k, n,
    {
        mf_cell_quadratic_calcslope(i, j, k, n,
                                    crsearr, crse_comp,
                                    slopearr,
                                    cdomain, bcrp);
    });

#if (AMREX_SPACEDIM == 2)
    if (crse_geom.IsRZ()) {

        // Get coarse and fine geometry data.
        GeometryData const& cs_geomdata = crse_geom.data();
        GeometryData const& fn_geomdata = fine_geom.data();

        // Compute fine correction.
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                               i, j, k, n,
        {
            mf_cell_quadratic_interp_rz(i, j, k, n,
                                        finearr, fine_comp,
                                        crsearr, crse_comp,
                                        cslopearr,
                                        ratio,
                                        cs_geomdata, fn_geomdata);
        });

    } else { /* crse_geom.IsCartesian() */
#endif /* AMREX_SPACEDIM == 2 */

        // No need for fine geometry data if using Cartesian coordinates.
        amrex::ignore_unused(fine_geom);

        // Compute fine correction.
        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                               i, j, k, n,
        {
            mf_cell_quadratic_interp(i, j, k, n,
                                     finearr, fine_comp,
                                     crsearr, crse_comp,
                                     cslopearr,
                                     ratio);
        });

#if (AMREX_SPACEDIM == 2)
    } // geom
#endif /* AMREX_SPACEDIM == 2 */

#endif /*(AMREX_SPACEDIM == 1)*/
}

Box
PCInterp::CoarseBox (const Box& fine,
                     int        ratio)
{
    return amrex::coarsen(fine,ratio);
}

Box
PCInterp::CoarseBox (const Box&     fine,
                     const IntVect& ratio)
{
    return amrex::coarsen(fine,ratio);
}

void
PCInterp::interp (const FArrayBox& crse,
                  int              crse_comp,
                  FArrayBox&       fine,
                  int              fine_comp,
                  int              ncomp,
                  const Box&       fine_region,
                  const IntVect&   ratio,
                  const Geometry& /*crse_geom*/,
                  const Geometry& /*fine_geom*/,
                  Vector<BCRec> const& /*bcr*/,
                  int               /*actual_comp*/,
                  int               /*actual_state*/,
                  RunOn             runon)
{
    BL_PROFILE("PCInterp::interp()");

    Array4<Real const> const& crsearr = crse.const_array();
    Array4<Real> const& finearr = fine.array();;

    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG ( runon, fine_region, tbx,
    {
        amrex::pcinterp_interp(tbx,finearr,fine_comp,ncomp,crsearr,crse_comp,ratio);
    });
}

CellConservativeProtected::CellConservativeProtected ()
    : CellConservativeLinear(true) {}

void
CellConservativeProtected::protect (const FArrayBox& /*crse*/,
                                    int              /*crse_comp*/,
                                    FArrayBox&       fine,
                                    int              /*fine_comp*/,
                                    FArrayBox&       fine_state,
                                    int              /*state_comp*/,
                                    int              ncomp,
                                    const Box&       fine_region,
                                    const IntVect&   ratio,
                                    const Geometry&  crse_geom,
                                    const Geometry&  fine_geom,
                                    Vector<BCRec>&   /*bcr*/,
                                    RunOn            runon)
{
#if (AMREX_SPACEDIM == 1)
    amrex::ignore_unused(fine,fine_state,
                         ncomp,fine_region,ratio,
                         crse_geom,fine_geom,runon);
    amrex::Abort("1D CellConservativeProtected::protect not supported");
#else
    BL_PROFILE("CellConservativeProtected::protect()");

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    //
    // crse_bx is coarsening of target_fine_region, grown by 1.
    //
    Box crse_bx = CoarseBox(target_fine_region,ratio);

    //
    // cs_bx is coarsening of target_fine_region.
    //
    Box cs_bx(crse_bx);
    cs_bx.grow(-1);

#if (AMREX_SPACEDIM == 2)
    /*
     * Get coarse and fine geometry data.
     */
    GeometryData cs_geomdata = crse_geom.data();
    GeometryData fn_geomdata = fine_geom.data();
#else
    amrex::ignore_unused(crse_geom, fine_geom);
#endif

    // Extract box from fine fab
    const Box& fnbx = fine.box();

    // Extract pointers to fab data
    Array4<Real>       const&   fnarr = fine.array();
    Array4<Real const> const& fnstarr = fine_state.const_array();

    /*
     * Loop over coarse indices.
     */
#if (AMREX_SPACEDIM == 2)
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cs_bx, ic, jc, kc,
    {
        ccprotect_2d(ic, jc, kc, ncomp,
                     fnbx, ratio,
                     cs_geomdata, fn_geomdata,
                     fnarr, fnstarr);
    }); // cs_bx
#else
    AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(runon, cs_bx, ic, jc, kc,
    {
        ccprotect_3d(ic, jc, kc, ncomp,
                     fnbx, ratio,
                     fnarr, fnstarr);
    }); // cs_bx
#endif

#endif /*(AMREX_SPACEDIM == 1)*/

}

Box
CellConservativeQuartic::CoarseBox (const Box& fine,
                                    int        ratio)
{
    Box crse(amrex::coarsen(fine,ratio));
    crse.grow(2);
    return crse;
}

Box
CellConservativeQuartic::CoarseBox (const Box&     fine,
                                    const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(2);
    return crse;
}

void
CellConservativeQuartic::interp (const FArrayBox&  crse,
                                 int               crse_comp,
                                 FArrayBox&        fine,
                                 int               fine_comp,
                                 int               ncomp,
                                 const Box&        fine_region,
                                 const IntVect&    ratio,
                                 const Geometry&   /* crse_geom */,
                                 const Geometry&   /* fine_geom */,
                                 Vector<BCRec> const& /*bcr*/,
                                 int               /* actual_comp */,
                                 int               /* actual_state */,
                                 RunOn             runon)
{
    BL_PROFILE("CellConservativeQuartic::interp()");
    BL_ASSERT(ratio[0] == 2);
#if (AMREX_SPACEDIM >= 2)
    BL_ASSERT(ratio[0] == ratio[1]);
#endif
#if (AMREX_SPACEDIM == 3)
    BL_ASSERT(ratio[1] == ratio[2]);
#endif
    amrex::ignore_unused(ratio);

    //
    // Make box which is intersection of fine_region and domain of fine.
    //
    Box target_fine_region = fine_region & fine.box();

    // Extract pointers to fab data
    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp, i, j, k, n,
    {
        ccquartic_interp(i, j, k, n,
                         crsearr, finearr);
    });
}

Box
FaceDivFree::CoarseBox (const Box& fine,
                        int        ratio)
{
    Box b = amrex::coarsen(fine,ratio).grow(1);
    return b;
}

Box
FaceDivFree::CoarseBox (const Box&     fine,
                        const IntVect& ratio)
{
    Box b = amrex::coarsen(fine,ratio).grow(1);
    return b;
}

void
FaceDivFree::interp (const FArrayBox&  /*crse*/,
                     int               /*crse_comp*/,
                     FArrayBox&        /*fine*/,
                     int               /*fine_comp*/,
                     int               /*ncomp*/,
                     const Box&        /*fine_region*/,
                     const IntVect&    /*ratio*/,
                     const Geometry&   /*crse_geom*/,
                     const Geometry&   /*fine_geom*/,
                     Vector<BCRec> const& /*bcr*/,
                     int               /*actual_comp*/,
                     int               /*actual_state*/,
                     RunOn             /*runon*/)
{
    amrex::Abort("FaceDivFree does not work on a single MultiFab. Call 'interp_arr' instead.");
}

void
FaceDivFree::interp_arr (Array<FArrayBox*, AMREX_SPACEDIM> const& crse,
                         const int         crse_comp,
                         Array<FArrayBox*, AMREX_SPACEDIM> const& fine,
                         const int         fine_comp,
                         const int         ncomp,
                         const Box&        fine_region,
                         const IntVect&    ratio,
                         Array<IArrayBox*, AMREX_SPACEDIM> const& solve_mask,
                         const Geometry&   /*crse_geom */,
                         const Geometry&   fine_geom,
                         Vector<Array<BCRec, AMREX_SPACEDIM> > const& /*bcr*/,
                         const int         /*actual_comp*/,
                         const int         /*actual_state*/,
                         const RunOn       runon)
{
    BL_PROFILE("FaceDivFree::interp()");

    Array<IndexType, AMREX_SPACEDIM> types;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        { types[d].set(d); }

    // This is currently only designed for octree, where ratio = 2.
    AMREX_ALWAYS_ASSERT(ratio == 2);

    const Box c_fine_region = amrex::coarsen(fine_region, ratio);
    GpuArray<Real, AMREX_SPACEDIM> cell_size = fine_geom.CellSizeArray();

    GpuArray<Array4<const Real>, AMREX_SPACEDIM> crsearr;
    GpuArray<Array4<Real>, AMREX_SPACEDIM> finearr;
    GpuArray<Array4<const int>, AMREX_SPACEDIM> maskarr;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
    {
        crsearr[d] = crse[d]->const_array(crse_comp);
        finearr[d] = fine[d]->array(fine_comp);
        if (solve_mask[d] != nullptr)
            { maskarr[d] = solve_mask[d]->const_array(0); }
    }

    // Fuse the launches, 1 for each dimension, into a single launch.
    AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(runon,
              amrex::convert(c_fine_region,types[0]), bx0,
              {
                  AMREX_LOOP_3D(bx0, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          amrex::facediv_face_interp<Real> (i,j,k,crse_comp+n,fine_comp+n, 0,
                                                            crsearr[0], finearr[0], maskarr[0], ratio);
                      }
                  });
              },
              amrex::convert(c_fine_region,types[1]), bx1,
              {
                  AMREX_LOOP_3D(bx1, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          amrex::facediv_face_interp<Real> (i,j,k,crse_comp+n,fine_comp+n, 1,
                                                            crsearr[1], finearr[1], maskarr[1], ratio);
                      }
                  });
              },
              amrex::convert(c_fine_region,types[2]), bx2,
              {
                  AMREX_LOOP_3D(bx2, i, j, k,
                  {
                      for (int n=0; n<ncomp; ++n)
                      {
                          amrex::facediv_face_interp<Real> (i,j,k,crse_comp+n,fine_comp+n, 2,
                                                            crsearr[2], finearr[2], maskarr[2], ratio);
                      }
                  });
              });

    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon,c_fine_region,ncomp,i,j,k,n,
    {
        amrex::facediv_int<Real>(i, j, k, fine_comp+n, finearr, ratio, cell_size);
    });
}

Box
CellQuartic::CoarseBox (const Box& fine, const IntVect& ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(2);
    return crse;
}

Box
CellQuartic::CoarseBox (const Box& fine, int ratio)
{
    Box crse = amrex::coarsen(fine,ratio);
    crse.grow(2);
    return crse;
}

void
CellQuartic::interp (const FArrayBox& crse,
                     int              crse_comp,
                     FArrayBox&       fine,
                     int              fine_comp,
                     int              ncomp,
                     const Box&       fine_region,
                     const IntVect&   ratio,
                     const Geometry&  /*crse_geom*/,
                     const Geometry&  /*fine_geom*/,
                     Vector<BCRec> const&  /*bcr*/,
                     int              /* actual_comp */,
                     int              /* actual_state */,
                     RunOn            runon)
{
    BL_PROFILE("CellQuartic::interp()");
    amrex::ignore_unused(ratio);
    AMREX_ASSERT(ratio == 2);

    Box target_fine_region = fine_region & fine.box();

    bool run_on_gpu = (runon == RunOn::Gpu && Gpu::inLaunchRegion());
    amrex::ignore_unused(run_on_gpu);

    Array4<Real const> const& crsearr = crse.const_array(crse_comp);
    Array4<Real>       const& finearr = fine.array(fine_comp);

#if (AMREX_SPACEDIM == 3)
    Box bz = amrex::coarsen(target_fine_region, IntVect(2,2,1));
    bz.grow(IntVect(2,2,0));
    FArrayBox tmpz(bz, ncomp);
    Elixir tmpz_eli;
    if (run_on_gpu) tmpz_eli = tmpz.elixir();
    Array4<Real> const& tmpzarr = tmpz.array();
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, bz, ncomp, i, j, k, n,
    {
        cell_quartic_interp_z(i,j,k,n,tmpzarr,crsearr);
    });
#endif

#if (AMREX_SPACEDIM >= 2)
    Box by = amrex::coarsen(target_fine_region, IntVect(AMREX_D_DECL(2,1,1)));
    by.grow(IntVect(AMREX_D_DECL(2,0,0)));
    FArrayBox tmpy(by, ncomp);
    Elixir tmpy_eli;
    if (run_on_gpu) tmpy_eli = tmpy.elixir();
    Array4<Real> const& tmpyarr = tmpy.array();
#if (AMREX_SPACEDIM == 2)
    Array4<Real const> srcarr = crsearr;
#else
    Array4<Real const> srcarr = tmpz.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, by, ncomp, i, j, k, n,
    {
        cell_quartic_interp_y(i,j,k,n,tmpyarr,srcarr);
    });
#endif

#if (AMREX_SPACEDIM == 1)
    Array4<Real const> srcarr = crsearr;
#else
    srcarr = tmpy.const_array();
#endif
    AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, target_fine_region, ncomp,
                                           i, j, k, n,
    {
        cell_quartic_interp_x(i,j,k,n,finearr,srcarr);
    });
}

}
