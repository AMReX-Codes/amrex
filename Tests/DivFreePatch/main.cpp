#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

#include <AMReX_FillPatchUtil.H>

using namespace amrex;
void main_main ();

void setupMF(MultiFab* mf)
{
    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        const Box& bx_x = mfi.fabbox();
        Array4<Real> c_x = mf->array(mfi);

        amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
//            c_x(i,j,k) = amrex::Random()*10;
            c_x(i,j,k) = double(i)+double(j)+double(k);
//            c_x(i,j,k) = double(i)*double(i)+double(j)*double(j)+double(k)*double(k);
        });
    }
}

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
    amrex::Vector<int> c_cell_3d(AMREX_SPACEDIM, 32);
    amrex::Vector<int> f_cell_3d(AMREX_SPACEDIM, 28);
    int max_grid_size = 64;

    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("f_offset", f_offset);
        pp.queryarr("c_cell_3d", c_cell_3d, 0, AMREX_SPACEDIM);
        pp.queryarr("f_cell_3d", f_cell_3d, 0, AMREX_SPACEDIM);
        pp.query("max_grid_size", max_grid_size);

        // inputs hierarchy:
        // n_cell > c_cell_3d & f_cell_3d

        if (n_cell != 0)
        {
            for (int i=0; i < AMREX_SPACEDIM; ++i)
            { c_cell_3d[i] = n_cell-1;
              f_cell_3d[i] = n_cell-f_offset-1; }
        }
    }

    int ncomp = 1;
    IntVect ratio{AMREX_D_DECL(2,2,2)};  // For this stencil (octree), always 2.
    IntVect ghost_c{AMREX_D_DECL(1,1,1)};  // For this stencil (octree), need 1 coarse ghost.
    IntVect ghost_f{AMREX_D_DECL(2,2,2)};  // For this stencil (octree), need 2 fine ghost.
    Geometry c_geom, f_geom;

    Array<MultiFab, AMREX_SPACEDIM> c_mf_faces;
    Array<MultiFab, AMREX_SPACEDIM> f_mf_faces;
    MultiFab div_coarse, div_refined_coarse, div_fine;
 
    AMREX_D_TERM( IntVect x_face{AMREX_D_DECL(1,0,0)};,
                  IntVect y_face{AMREX_D_DECL(0,1,0)};,
                  IntVect z_face{AMREX_D_DECL(0,0,1)};  );

    //  Create multifabs.
    {
        Box domain  (IntVect{0}, IntVect{c_cell_3d});
        Box domain_f(IntVect{f_offset}, IntVect{f_cell_3d});
        domain_f.refine(ratio);

        double scale = double(f_cell_3d[0]+1)/double(c_cell_3d[0]+1);
        RealBox realbox_c({AMREX_D_DECL(-1.0,-1.0,-1.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
        RealBox realbox_f({AMREX_D_DECL(-scale,-scale,-scale)}, {AMREX_D_DECL(scale,scale,scale)});
        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        // Build coarse and fine boxArrays and DistributionMappings.
        BoxArray ba_c(domain);
        ba_c.maxSize(max_grid_size);

        BoxArray ba_f(domain_f);
        ba_f.maxSize(max_grid_size);

        DistributionMapping dm_c(ba_c);
        DistributionMapping dm_f(ba_f);

        c_geom.define(domain,   realbox_c, CoordSys::cartesian, is_periodic);
        f_geom.define(domain_f, realbox_f, CoordSys::cartesian, is_periodic);

        AMREX_D_TERM( c_mf_faces[0].define( amrex::convert( ba_c,x_face ), dm_c, ncomp, ghost_c);,
                      c_mf_faces[1].define( amrex::convert( ba_c,y_face ), dm_c, ncomp, ghost_c);,
                      c_mf_faces[2].define( amrex::convert( ba_c,z_face ), dm_c, ncomp, ghost_c); );

        AMREX_D_TERM( f_mf_faces[0].define( amrex::convert( ba_f,x_face ), dm_f, ncomp, ghost_f);,
                      f_mf_faces[1].define( amrex::convert( ba_f,y_face ), dm_f, ncomp, ghost_f);,
                      f_mf_faces[2].define( amrex::convert( ba_f,z_face ), dm_f, ncomp, ghost_f); );

        AMREX_D_TERM( f_mf_faces[0].setVal(0.0);,
                      f_mf_faces[1].setVal(0.0);,
                      f_mf_faces[2].setVal(0.0);  );

        div_coarse.define(ba_c, dm_c, ncomp, ghost_c);
        div_coarse.setVal(0.0);

        div_refined_coarse.define(ba_f, dm_f, ncomp, ghost_f);
        div_refined_coarse.setVal(0.0);

        div_fine.define(ba_f, dm_f, ncomp, ghost_f);
        div_fine.setVal(0.0);

        amrex::Print() << "Testing Face FreeDiv FillPatch with: "
                       << "\n  dimensions = "    << ba_c.minimalBox()
                       << "\n  max_grid_size = " << max_grid_size
                       << "\n  boxes = "         << ba_c.size()
                       << "\n  and ratio = "     << ratio << std::endl;

        amrex::Print() << "Coarse box array: " << ba_c << std::endl;
        amrex::Print() << "Fine box array: " << ba_f << std::endl;
        amrex::Print() << "============================" << std::endl;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i)
    {
        setupMF(&c_mf_faces[i]);
    }

    amrex::UtilCreateDirectoryDestructive("pltfiles");
    AMREX_D_TERM( amrex::VisMF::Write(c_mf_faces[0], std::string("pltfiles/cx"));,
                  amrex::VisMF::Write(c_mf_faces[1], std::string("pltfiles/cy"));,
                  amrex::VisMF::Write(c_mf_faces[2], std::string("pltfiles/cz"));  );

// ***************************************************************

    amrex::Print() << " Calculating coarse divergence. " << std::endl;
    for (MFIter mfi(div_coarse); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> div = div_coarse.array(mfi);

        AMREX_D_TERM( const Array4<Real>& face_x = c_mf_faces[0].array(mfi);,
                      const Array4<Real>& face_y = c_mf_faces[1].array(mfi);,
                      const Array4<Real>& face_z = c_mf_faces[2].array(mfi);  );

        AMREX_D_TERM( const Real& dx = c_geom.CellSize(0);,
                      const Real& dy = c_geom.CellSize(1);,
                      const Real& dz = c_geom.CellSize(2);  );

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            div(i,j,k) = AMREX_D_TERM(   (face_x(i+1,j  ,k  ) - face_x(i,j,k))/dx,
                                       + (face_y(i  ,j+1,k  ) - face_y(i,j,k))/dy,
                                       + (face_z(i  ,j  ,k+1) - face_z(i,j,k))/dz  );
        });
    }
    amrex::VisMF::Write(div_coarse, std::string("pltfiles/coarse"));

    amrex::Print() << " Copying coarse divergence to fine grid. " << std::endl;
    {
        double time = 0;

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

        InterpFromCoarseLevel(div_refined_coarse, time,
                              div_coarse, 0, 0, 1,
                              c_geom, f_geom,
                              phys_bc, 0, phys_bc, 0,
                              ratio, mapper, bcrec, 0);
    }
    amrex::VisMF::Write(div_refined_coarse, std::string("pltfiles/coarsetofine"));

// ***************************************************************

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

        Vector<Array<MultiFab*, AMREX_SPACEDIM> > fine_v;
        Vector<Array<MultiFab*, AMREX_SPACEDIM> > coarse_v;
        fine_v.push_back(fine_faces);
        coarse_v.push_back(coarse_faces);

        Array<PhysBCFunctNoOp, AMREX_SPACEDIM> phys_bc;

        InterpFromCoarseLevel(fine_faces, time,
                              coarse_faces, 0, 0, 1,
                              c_geom, f_geom,
                              phys_bc, 0, phys_bc, 0,
                              ratio, mapper, bcrec, 0);
    }

    AMREX_D_TERM( amrex::VisMF::Write(f_mf_faces[0], std::string("pltfiles/fx"));,
                  amrex::VisMF::Write(f_mf_faces[1], std::string("pltfiles/fy"));,
                  amrex::VisMF::Write(f_mf_faces[2], std::string("pltfiles/fz"));  );

// ***************************************************************

    amrex::Print() << " Calculating Fine Divergence. " << std::endl;
    for (MFIter mfi(div_fine); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> div = div_fine.array(mfi);

        AMREX_D_TERM( Array4<Real> face_x = f_mf_faces[0].array(mfi);,
                      Array4<Real> face_y = f_mf_faces[1].array(mfi);,
                      Array4<Real> face_z = f_mf_faces[2].array(mfi);  );

        AMREX_D_TERM( const Real& dx = f_geom.CellSize(0);,
                      const Real& dy = f_geom.CellSize(1);,
                      const Real& dz = f_geom.CellSize(2);  );

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            div(i,j,k) = AMREX_D_TERM(   (face_x(i+1,j  ,k  ) - face_x(i,j,k))/dx,
                                       + (face_y(i  ,j+1,k  ) - face_y(i,j,k))/dy,
                                       + (face_z(i  ,j  ,k+1) - face_z(i,j,k))/dz  );
        });
    }
    amrex::VisMF::Write(div_fine, std::string("pltfiles/fine"));

    amrex::Print() << " Checking Coarse vs. Fine Divergence. " << std::endl;
    {
        div_fine.minus(div_refined_coarse, 0, 1, 0);
        amrex::VisMF::Write(div_fine, std::string("pltfiles/diff"));
        Real max_error = div_fine.max(0);

        amrex::Print() << " Divergence error = " << max_error << std::endl;
    }


/*
    amrex::Print() << " Performing DivFree FillPatch. " << std::endl;
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

        Vector<Array<MultiFab*, AMREX_SPACEDIM> > fine_v;
        Vector<Array<MultiFab*, AMREX_SPACEDIM> > coarse_v;
        fine_v.push_back(fine_faces);
        coarse_v.push_back(coarse_faces);

        Array<PhysBCFunctNoOp, AMREX_SPACEDIM> phys_bc;

        amrex::Print() << " Starting FillPatch. " << std::endl;

        Geometry f_total_geom = c_geom;
        f_total_geom.refine(ratio);

        FillPatchTwoLevels(fine_faces, ghost, time,
                           coarse_v, time_v,
                           fine_v, time_v,
                           0, 0, 1, c_geom, f_total_geom,
                           phys_bc, 0, phys_bc, 0,
                           ratio, mapper, bcrec, 0);

    }
*/

// ***************************************************************



}
