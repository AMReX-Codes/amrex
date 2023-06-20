#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;
void main_main ();

// ================================================

void setupMF(MultiFab& mf, const int type = 0, const BoxArray& exclude = BoxArray())
{
    if ((type < 0) || (type > 2))
        { amrex::Abort("Invalid setup type."); }

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        Box fbx = mfi.fabbox();
        Array4<Real> arr = mf.array(mfi);

        BoxArray ba;
        ba = amrex::complementIn(fbx, exclude);

        for (int bid=0; bid<ba.size(); ++bid)
        {
            Box bx = ba[bid];

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (type == 0)
                    { arr(i,j,k) = amrex::Random()*10; }
                else if (type == 1)
                    { arr(i,j,k) = double(i)+double(j)+double(k); }
                else if (type == 2)
                    { arr(i,j,k) = double(i)*double(i)+double(j)*double(j)+double(k)*double(k); }
            });
        }
    }
}

void calcDiv(Array<MultiFab, AMREX_SPACEDIM>& faces,
             MultiFab& mf_div,
             const GpuArray<Real, AMREX_SPACEDIM> cellsize)
{
    for (MFIter mfi(mf_div); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.fabbox();
        Array4<Real> div = mf_div.array(mfi);

        AMREX_D_TERM( const Array4<Real>& face_x = faces[0].array(mfi);,
                      const Array4<Real>& face_y = faces[1].array(mfi);,
                      const Array4<Real>& face_z = faces[2].array(mfi);  );

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            div(i,j,k) = AMREX_D_TERM(   (face_x(i+1,j  ,k  ) - face_x(i,j,k))/cellsize[0],
                                       + (face_y(i  ,j+1,k  ) - face_y(i,j,k))/cellsize[1],
                                       + (face_z(i  ,j  ,k+1) - face_z(i,j,k))/cellsize[2]);
        });
    }
}

void CoarsenToFine(MultiFab& div_refined_coarse,
                   const MultiFab& div_coarse,
                   const Geometry& c_geom, const Geometry& f_geom,
                   IntVect ratio)
{
    const int ncomp = div_coarse.nComp();
    Vector<BCRec> bcrec(ncomp);
    for (int n = 0; n < ncomp; ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            bcrec[n].setLo(idim, BCType::int_dir);
            bcrec[n].setHi(idim, BCType::int_dir);
        }
    }

    Interpolater* mapper = &pc_interp;
    PhysBCFunctNoOp phys_bc;

    const int time = 0; // No time extrapolation, so just use 0.

    InterpFromCoarseLevel(div_refined_coarse, time,
                          div_coarse, 0, 0, ncomp,
                          c_geom, f_geom,
                          phys_bc, 0, phys_bc, 0,
                          ratio, mapper, bcrec, 0);
}


Real MFdiff(const MultiFab& lhs, const MultiFab& rhs,
            int strt_comp, int num_comp, int nghost, const std::string name = "")
{
    MultiFab temp(lhs.boxArray(), lhs.DistributionMap(), lhs.nComp(), nghost);
    Copy(temp, lhs, strt_comp, strt_comp, num_comp, nghost);
    temp.minus(rhs, strt_comp, num_comp, nghost);

    if (name != "")
        { amrex::VisMF::Write(temp, std::string("pltfiles/" + name)); }

    Real max_diff = 0;
    for (int i=0; i<num_comp; ++i)
    {
        Real max_i = std::abs(temp.max(i));
        max_diff = (max_diff > max_i) ? max_diff : max_i;
    }

    return max_diff;
}

// ================================================

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
}

void main_main ()
{
    BL_PROFILE("main");

    int n_cell = 0;
    int f_offset = 4;
    int nghost_c = 1;
    int nghost_f = 2;

    amrex::Vector<int> c_lo(AMREX_SPACEDIM,  0);
    amrex::Vector<int> c_hi(AMREX_SPACEDIM, 32);
    amrex::Vector<int> f_lo(AMREX_SPACEDIM, 28);
    amrex::Vector<int> f_hi(AMREX_SPACEDIM,  4);
    int max_grid_size = 64;

    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("f_offset", f_offset);
        pp.query("max_grid_size", max_grid_size);
        pp.query("nghost_c", nghost_c);
        pp.query("nghost_f", nghost_f);

        pp.queryarr("c_hi",  c_hi, 0, AMREX_SPACEDIM);
        pp.queryarr("f_lo",  f_lo, 0, AMREX_SPACEDIM);
        pp.queryarr("f_hi",  f_hi, 0, AMREX_SPACEDIM);

        if (n_cell != 0)
        {
            for (int i=0; i < AMREX_SPACEDIM; ++i)
            { c_lo[i] = 0;
              c_hi[i] = n_cell-1;

              f_lo[i] = f_offset;
              f_hi[i] = n_cell-f_offset-1; }
        }
    }

    int ncomp = 1;
    IntVect ratio{AMREX_D_DECL(2,2,2)};    // For this stencil (octree), always 2.
    IntVect ghost_c{AMREX_D_DECL(nghost_c, nghost_c, nghost_c)};  // For this stencil (octree), need 1 coarse ghost.
    IntVect ghost_f{AMREX_D_DECL(nghost_f, nghost_f, nghost_f)};  // For this stencil (octree), need 1 fine ghost.
    Geometry c_geom, f_geom, f_geom_wghost, f_geom_all, f_geom_partial;

    Array<MultiFab, AMREX_SPACEDIM> c_mf_faces;
    Array<MultiFab, AMREX_SPACEDIM> f_mf_faces;
    Array<MultiFab, AMREX_SPACEDIM> f_mf_copy;
    MultiFab div_coarse, div_refined_coarse, div_fine;

    // For outputting ghost cells for debugging.
    Array<MultiFab, AMREX_SPACEDIM> f_mf_faces_wg;
    MultiFab div_fine_wg;

    AMREX_D_TERM( IntVect x_face{AMREX_D_DECL(1,0,0)};,
                  IntVect y_face{AMREX_D_DECL(0,1,0)};,
                  IntVect z_face{AMREX_D_DECL(0,0,1)};  );

// ***************************************************************
    // Build the Multifabs and Geometries.
    {
        Box domain   (IntVect{c_lo}, IntVect{c_hi});
        Box domain_f (IntVect{f_lo}, IntVect{f_hi});
        Box domain_fg(domain_f);

        amrex::Print() << " Testing on coarse: " << domain << std::endl;
        amrex::Print() << "  w/ fine area covering: " << domain_f << std::endl;

        domain_f.refine(ratio);
        domain_fg.refine(ratio);
        domain_fg.grow(ghost_f);

        const IntVect& crse_hi = domain.bigEnd()*ratio + (ratio-1);
        const IntVect& fine_lo = domain_f.smallEnd();
        const IntVect& fine_hi = domain_f.bigEnd();

        const IntVect& fine_len = domain_f.size();
        const IntVect& fine_lo_partial = fine_lo+(fine_len/6);
        const IntVect& fine_hi_partial = fine_hi-(fine_len/3);

        Box domain_p(fine_lo_partial, fine_hi_partial);
        amrex::Print() << "Partial region: " << domain_p << std::endl;

        RealBox realbox_c    ({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
        RealBox realbox_f_all({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
        RealBox realbox_f({AMREX_D_DECL( double(fine_lo[0])   / double(crse_hi[0]+1),
                                         double(fine_lo[1])   / double(crse_hi[1]+1),
                                         double(fine_lo[2])   / double(crse_hi[2]+1) )},
                          {AMREX_D_DECL( double(fine_hi[0]+1) / double(crse_hi[0]+1),
                                         double(fine_hi[1]+1) / double(crse_hi[1]+1),
                                         double(fine_hi[2]+1) / double(crse_hi[2]+1) )} );
        RealBox realbox_fg({AMREX_D_DECL( double(fine_lo[0]-ghost_f[0])   / double(crse_hi[0]+1),
                                          double(fine_lo[1]-ghost_f[1])   / double(crse_hi[1]+1),
                                          double(fine_lo[2]-ghost_f[2])   / double(crse_hi[2]+1) )},
                           {AMREX_D_DECL( double(fine_hi[0]+ghost_f[0]+1) / double(crse_hi[0]+1),
                                          double(fine_hi[1]+ghost_f[1]+1) / double(crse_hi[1]+1),
                                          double(fine_hi[2]+ghost_f[1]+1) / double(crse_hi[2]+1) )} );
        RealBox realbox_fp({AMREX_D_DECL( double(fine_lo_partial[0])   / double (crse_hi[0]+1),
                                          double(fine_lo_partial[1])   / double (crse_hi[1]+1),
                                          double(fine_lo_partial[2])   / double (crse_hi[2]+1) )},
                           {AMREX_D_DECL( double(fine_hi_partial[0]+1) / double (crse_hi[0]+1),
                                          double(fine_hi_partial[1]+1) / double (crse_hi[1]+1),
                                          double(fine_hi_partial[2]+1) / double (crse_hi[2]+1) )} );

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        // Build coarse and fine boxArrays and DistributionMappings.
        BoxArray ba_c(domain);
        ba_c.maxSize(max_grid_size);

        BoxArray ba_f(domain_f);
        ba_f.maxSize(max_grid_size);

        BoxArray ba_fg(domain_fg);
        ba_fg.maxSize(max_grid_size);

        DistributionMapping dm_c(ba_c);
        DistributionMapping dm_f(ba_f);
        DistributionMapping dm_fg(ba_fg);

        c_geom.define        (domain,    realbox_c,  CoordSys::cartesian, is_periodic);
        f_geom.define        (domain_f,  realbox_f,  CoordSys::cartesian, is_periodic);
        f_geom_wghost.define (domain_fg, realbox_fg, CoordSys::cartesian, is_periodic);
        f_geom_all.define    (amrex::refine(domain, ratio), realbox_f_all, CoordSys::cartesian, is_periodic);
        f_geom_partial.define(domain_p,  realbox_fp, CoordSys::cartesian, is_periodic);

        AMREX_D_TERM( c_mf_faces[0].define( amrex::convert( ba_c,x_face ), dm_c, ncomp, ghost_c);,
                      c_mf_faces[1].define( amrex::convert( ba_c,y_face ), dm_c, ncomp, ghost_c);,
                      c_mf_faces[2].define( amrex::convert( ba_c,z_face ), dm_c, ncomp, ghost_c); );

        AMREX_D_TERM( f_mf_faces[0].define( amrex::convert( ba_f,x_face ), dm_f, ncomp, ghost_f);,
                      f_mf_faces[1].define( amrex::convert( ba_f,y_face ), dm_f, ncomp, ghost_f);,
                      f_mf_faces[2].define( amrex::convert( ba_f,z_face ), dm_f, ncomp, ghost_f); );

        AMREX_D_TERM( f_mf_copy[0].define( amrex::convert( ba_f,x_face ), dm_f, ncomp, ghost_f);,
                      f_mf_copy[1].define( amrex::convert( ba_f,y_face ), dm_f, ncomp, ghost_f);,
                      f_mf_copy[2].define( amrex::convert( ba_f,z_face ), dm_f, ncomp, ghost_f); );

        AMREX_D_TERM( f_mf_faces[0].setVal(0.0);,
                      f_mf_faces[1].setVal(0.0);,
                      f_mf_faces[2].setVal(0.0);  );

        div_coarse.define(ba_c, dm_c, ncomp, ghost_c);
        div_coarse.setVal(0.0);

        div_refined_coarse.define(ba_f, dm_f, ncomp, ghost_f);
        div_refined_coarse.setVal(0.0);

        div_fine.define(ba_f, dm_f, ncomp, ghost_f);
        div_fine.setVal(0.0);

        AMREX_D_TERM( f_mf_faces_wg[0].define( amrex::convert( ba_fg,x_face ), dm_fg, ncomp, 0 );,
                      f_mf_faces_wg[1].define( amrex::convert( ba_fg,y_face ), dm_fg, ncomp, 0 );,
                      f_mf_faces_wg[2].define( amrex::convert( ba_fg,z_face ), dm_fg, ncomp, 0 ); );

        AMREX_D_TERM( f_mf_faces_wg[0].setVal(0.0);,
                      f_mf_faces_wg[1].setVal(0.0);,
                      f_mf_faces_wg[2].setVal(0.0);  );

        div_fine_wg.define(ba_fg, dm_fg, ncomp, 0);
        div_fine_wg.setVal(0.0);

        amrex::Print() << " Testing Face FreeDiv FillPatch with: "
                       << " \n  dimensions = "    << ba_c.minimalBox()
                       << " \n  max_grid_size = " << max_grid_size
                       << " \n  boxes = "         << ba_c.size()
                       << " \n  and ratio = "     << ratio << std::endl;

        amrex::Print() << " Coarse box array: " << ba_c << std::endl;
        amrex::Print() << " Fine box array: " << ba_f << std::endl;
        amrex::Print() << " Fine box w/ ghosts array: " << ba_fg << std::endl;
        amrex::Print() << "============================" << std::endl;
    }

// ***************************************************************
//  Setup initial value on the coarse faces.
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        setupMF(c_mf_faces[i], 1);
    }

    amrex::UtilCreateDirectoryDestructive("pltfiles");
    AMREX_D_TERM( amrex::VisMF::Write(c_mf_faces[0], std::string("pltfiles/cx"));,
                  amrex::VisMF::Write(c_mf_faces[1], std::string("pltfiles/cy"));,
                  amrex::VisMF::Write(c_mf_faces[2], std::string("pltfiles/cz"));  );

// ***************************************************************
//  Calculate divergence on the coarse grid and copy it to a fine grid.
//  This is the target divergence for the final result.

    amrex::Print() << " Calculating coarse divergence. " << std::endl;
    calcDiv(c_mf_faces, div_coarse, c_geom.CellSizeArray());
    amrex::VisMF::Write(div_coarse, std::string("pltfiles/coarse"));

    amrex::Print() << " Copying coarse divergence to fine grid. " << std::endl;
    CoarsenToFine(div_refined_coarse, div_coarse, c_geom, f_geom_all, ratio);
    amrex::VisMF::Write(div_refined_coarse, std::string("pltfiles/coarsetofine"));

    div_fine_wg.ParallelCopy(div_refined_coarse, 0, 0, 1, ghost_f, IntVect::TheZeroVector());
    amrex::VisMF::Write(div_fine_wg, std::string("pltfiles/coarsetofine_wg"));

// ***************************************************************
//  Interp initial coarse values to the fine grid.

    amrex::Print() << " Starting InterpFromCoarse. " << std::endl;
    {
        double time = 1;
        Vector<Real> time_v;
        time_v.push_back(time);

        Array<Vector<BCRec>, AMREX_SPACEDIM> bcrec;
        for (int odim=0; odim < AMREX_SPACEDIM; ++odim)
        {
            bcrec[odim].resize(ncomp);
            for (int n = 0; n < ncomp; ++n)
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    bcrec[odim][n].setLo(idim, BCType::foextrap);
                    bcrec[odim][n].setHi(idim, BCType::foextrap);
                }
            }
        }

        Interpolater* mapper = &face_divfree_interp;
        Array<MultiFab*, AMREX_SPACEDIM> fine_faces;
        Array<MultiFab*, AMREX_SPACEDIM> coarse_faces;

        for (int i=0; i<AMREX_SPACEDIM; ++i)
        {
            fine_faces[i] = &(f_mf_faces[i]);
            coarse_faces[i] = &(c_mf_faces[i]);
        }

        Vector<Array<MultiFab*, AMREX_SPACEDIM> > fine_v;
        Vector<Array<MultiFab*, AMREX_SPACEDIM> > coarse_v;
        fine_v.push_back(fine_faces);
        coarse_v.push_back(coarse_faces);

        Array<PhysBCFunctNoOp, AMREX_SPACEDIM> phys_bc;

        InterpFromCoarseLevel(fine_faces, time,
                              coarse_faces, 0, 0, 1,
//                              c_geom, f_geom_partial,
                              c_geom, f_geom_all,
                              phys_bc, 0, phys_bc, 0,
                              ratio, mapper, bcrec, 0);
    }

    // Check for errors
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (f_mf_faces[i].contains_nan()) {
            amrex::Print() << "******** Nans present in fine velocity in dimension " << i << std::endl;
        }
    }

    AMREX_D_TERM( amrex::VisMF::Write(f_mf_faces[0], std::string("pltfiles/fx"));,
                  amrex::VisMF::Write(f_mf_faces[1], std::string("pltfiles/fy"));,
                  amrex::VisMF::Write(f_mf_faces[2], std::string("pltfiles/fz"));  );

// ***************************************************************
//  Check divergence on the fine grid, subtract the target
//      and report maximum value.

    amrex::Print() << " Calculating Fine Divergence. " << std::endl;
    calcDiv(f_mf_faces, div_fine, f_geom.CellSizeArray());

    div_fine_wg.ParallelCopy(div_fine, 0, 0, 1, ghost_f, IntVect::TheZeroVector());
    amrex::VisMF::Write(div_fine, std::string("pltfiles/fine"));
    amrex::VisMF::Write(div_fine_wg, std::string("pltfiles/finewg"));

    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
       f_mf_faces_wg[i].ParallelCopy(f_mf_faces[i], 0, 0, 1, ghost_f, IntVect::TheZeroVector());
    }

    AMREX_D_TERM( amrex::VisMF::Write(f_mf_faces_wg[0], std::string("pltfiles/fwgx"));,
                  amrex::VisMF::Write(f_mf_faces_wg[1], std::string("pltfiles/fwgy"));,
                  amrex::VisMF::Write(f_mf_faces_wg[2], std::string("pltfiles/fwgz"));  );

    amrex::Print() << " Max InterpFromCoarse divergence error: "
                   << MFdiff(div_fine, div_refined_coarse, 0, 1, nghost_f, "diff") << std::endl;

// ***************************************************************

    // Change coarse values, save current fine values for checking
    //   the final result.
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        setupMF(c_mf_faces[i], 2, BoxArray(amrex::coarsen(f_geom.Domain(), ratio).convert(c_mf_faces[i].ixType())));
        Copy(f_mf_copy[i], f_mf_faces[i], 0, 0, 1, 0);
    }

    AMREX_D_TERM( amrex::VisMF::Write(c_mf_faces[0], std::string("pltfiles/cxFP"));,
                  amrex::VisMF::Write(c_mf_faces[1], std::string("pltfiles/cyFP"));,
                  amrex::VisMF::Write(c_mf_faces[2], std::string("pltfiles/czFP"));  );

    calcDiv(c_mf_faces, div_coarse, c_geom.CellSizeArray());
    CoarsenToFine(div_refined_coarse, div_coarse, c_geom, f_geom_all, ratio);
    amrex::VisMF::Write(div_coarse, std::string("pltfiles/coarsetofineB"));

    amrex::Print() << " Checking new adjustment hasn't changed solution on fine region: "
                   << MFdiff(div_fine, div_refined_coarse, 0, 1, nghost_f, "diffFP") << std::endl;


// ***************************************************************

    // Call FillPatchTwoLevels to update fine ghost cells.
    amrex::Print() << std::endl;
    amrex::Print() << " ********************** " << std::endl;
    amrex::Print() << " Performing DivFree FillPatch. " << std::endl;
    {
        Real time = 1;
        Vector<Real> time_v(1,1);

        Array<Vector<BCRec>, AMREX_SPACEDIM> bcrec;
        for (int odim=0; odim < AMREX_SPACEDIM; ++odim)
        {
            bcrec[odim].resize(ncomp);
            for (int n = 0; n < ncomp; ++n)
            {
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    bcrec[odim][n].setLo(idim, BCType::int_dir);
                    bcrec[odim][n].setHi(idim, BCType::int_dir);
                }
            }
        }

        Interpolater* mapper = &face_divfree_interp;
        Array<MultiFab*, AMREX_SPACEDIM> fine_faces;
        Array<MultiFab*, AMREX_SPACEDIM> coarse_faces;

        for (int i=0; i<AMREX_SPACEDIM; ++i)
        {
            fine_faces[i] = &(f_mf_faces[i]);
            coarse_faces[i] = &(c_mf_faces[i]);
        }

        // Make sure coarse & fine are properly setup for the interpolation stencil.
        amrex::average_down_faces( {AMREX_D_DECL(&f_mf_faces[0], &f_mf_faces[1], &f_mf_faces[2])}, coarse_faces, ratio, 0);

        Vector<Array<MultiFab*, AMREX_SPACEDIM> > fine_v;
        Vector<Array<MultiFab*, AMREX_SPACEDIM> > coarse_v;
        fine_v.push_back(fine_faces);
        coarse_v.push_back(coarse_faces);

        Array<PhysBCFunctNoOp, AMREX_SPACEDIM> phys_bc;

        amrex::Print() << " Starting FillPatch. " << std::endl;

        FillPatchTwoLevels(fine_faces, time,
                           coarse_v, time_v,
                           fine_v, time_v,
                           0, 0, 1, c_geom, f_geom_all,
                           phys_bc, 0, phys_bc, 0,
                           ratio, mapper, bcrec, 0);
    }

// ================================================

    // Checking fine valid cells are identical.
    Real max_diff = 0;
    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        Real max_i = std::abs( MFdiff(f_mf_copy[i], f_mf_faces[i], 0, 1, 0) );
        max_diff = (max_diff > max_i) ? max_diff : max_i;
    }
    amrex::Print() << " Fine values maximum change: " << max_diff << std::endl;

    // Check fine divergence = coarse divergence in ghost cells.
    calcDiv(f_mf_faces, div_fine, f_geom.CellSizeArray());
    amrex::VisMF::Write(div_fine, std::string("pltfiles/fineFP"));

    amrex::Print() << " Max FillPatchTwoLevels divergence error: "
                   << MFdiff(div_fine, div_refined_coarse, 0, 1, nghost_f, "diffFP") << std::endl;

    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        f_mf_faces_wg[i].ParallelCopy(f_mf_faces[i], 0, 0, 1, ghost_f, IntVect::TheZeroVector());
    }

    // Check for errors
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (f_mf_faces_wg[i].contains_nan()) {
            amrex::Print() << "******** Nans present in fine velocity after FillPatch (including ghosts)  in dimension " << i << std::endl;
        }
    }

    AMREX_D_TERM( amrex::VisMF::Write(f_mf_faces_wg[0], std::string("pltfiles/fwgxFP"));,
                  amrex::VisMF::Write(f_mf_faces_wg[1], std::string("pltfiles/fwgyFP"));,
                  amrex::VisMF::Write(f_mf_faces_wg[2], std::string("pltfiles/fwgzFP"));  );

// ***************************************************************

}
