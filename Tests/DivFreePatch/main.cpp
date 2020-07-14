#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

#include <AMReX_FillPatchUtil.H>

using namespace amrex;

void main_main ();

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
    int n_boxes_per_rank = 0;
    amrex::Vector<int> n_cell_3d (AMREX_SPACEDIM, 512);
    int max_grid_size = 64;
    int nwork = 10;
    int nwrites = 4;

    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.queryarr("n_cell_3d", n_cell_3d, 0, AMREX_SPACEDIM);
        pp.query("n_boxes_per_rank", n_boxes_per_rank);
        pp.query("max_grid_size", max_grid_size);
        pp.query("nwork", nwork);
        pp.query("nwrites", nwrites);

        // inputs hierarchy: 
        // n_cell > n_boxes_per_rank > n_cell_3d

        if (n_cell != 0)
        {
            for (int i=0; i < AMREX_SPACEDIM; ++i)
            { n_cell_3d[i] = n_cell; }
        }
        else if (n_boxes_per_rank != 0)
        {
           n_cell_3d[0] = (max_grid_size) - 1;
           n_cell_3d[1] = (max_grid_size * n_boxes_per_rank) - 1;
           n_cell_3d[2] = (max_grid_size * ParallelDescriptor::NProcs()) - 1;
        }
    }

    BoxArray ba(Box(IntVect(0),IntVect(n_cell_3d)));
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    Array<MultiFab, AMREX_SPACEDIM> mf_faces;
    Gpu::ManagedVector< Array4<Real> > arrs(nwrites);

    // Build x, y and z face MFs.

    IntVect x_face(AMREX_D_DECL(1,0,0));
    IntVect y_face(AMREX_D_DECL(0,1,0));
#if (AMREX_SPACEDIM==3)
    IntVect z_face(AMREX_D_DECL(0,0,1));
#endif

    mf_faces[0].define( amrex::convert( ba,x_face ), dm, 1, 1);
    mf_faces[1].define( amrex::convert( ba,y_face ), dm, 1, 1);
#if (AMREX_SPACEDIM==3)
    mf_faces[2].define( amrex::convert( ba,z_face ), dm, 1, 1);
#endif

    for (MFIter mfi(mf_faces[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.fabbox();
        for (int m = 0; m < AMREX_SPACEDIM; ++m) {
            arrs[m] = mf_faces[m].array(mfi);
        }

        auto arrs_ptr = arrs.dataPtr();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int m = 0; m < AMREX_SPACEDIM; ++m) {
               arrs_ptr[m](i,j,k) = amrex::Random()*10;
            }
        });
        Gpu::streamSynchronize(); // because of arrs
    }

    amrex::Print() << "Testing Face FreeDiv FillPatch with: "
                   << "\n  dimensions = "    << ba.minimalBox() 
                   << "\n  max_grid_size = " << max_grid_size
                   << "\n  boxes = "         << ba.size()
                   << "\n  and nwork = "     << nwork << std::endl;

// ***************************************************************

//    Calc div
//    Do 2LevelFillPatch
//    Compare/check div at fine level

// ***************************************************************

    amrex::Print() << " Calculate divergence. " << std::endl;

    MultiFab div_coarse(ba, dm, 1, 0);
    div_coarse.setVal(0.0);

    for (MFIter mfi(div_coarse); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> div = div_coarse.array(mfi);
        Array4<Real> face_x = mf_faces[0].array(mfi);
        Array4<Real> face_y = mf_faces[1].array(mfi);
#if (AMREX_SPACEDIM==3)
        Array4<Real> face_z = mf_faces[2].array(mfi);
#endif

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            div(i,j,k) = AMREX_D_TERM(   face_x(i+1,j  ,k  ) - face_x(i,j,k),
                                       + face_y(i  ,j+1,k  ) - face_y(i,j,k),
                                       + face_z(i  ,j  ,k+1) - face_z(i,j,k)  );
        });
    }

// ***************************************************************

/*
    {
            const auto& bcrec = get_velocity_bcrec();
            PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > cphysbc
                (geom[lev-1], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
            PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > fphysbc
                (geom[lev], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
#else
            Interpolater* mapper = &cell_cons_interp;
#endif
            FillPatchTwoLevels(vel, IntVect(ng), time,
                               {&(m_leveldata[lev-1]->velocity_o),
                                &(m_leveldata[lev-1]->velocity)},
                               {m_t_old[lev-1], m_t_new[lev-1]},
                               {&(m_leveldata[lev]->velocity_o),
                                &(m_leveldata[lev]->velocity)},
                               {m_t_old[lev], m_t_new[lev]},
                               0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                               cphysbc, 0, fphysbc, 0,
                               refRatio(lev-1), mapper, bcrec, 0);
    }


    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM>& mf, IntVect const& nghost, Real time,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});


    Compute div in fine case & compare to coarse value in matching cell

    for (MFIter mfi(div_coarse); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> div = div_coarse.array(mfi);
        Array4<Real> face_x = mf_faces[0].array(mfi);
        Array4<Real> face_y = mf_faces[1].array(mfi);
        Array4<Real> face_z = mf_faces[2].array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            div(i,j,k) = ( face_x(i+1,j  ,k  ) - face_x(i,j,k) ) +
                         ( face_y(i  ,j+1,k  ) - face_y(i,j,k) ) +
                         ( face_z(i  ,j  ,k+1) - face_z(i,j,k) );
        });
    }

*/


/*

void incflo::fillpatch_velocity (int lev, Real time, MultiFab& vel, int ng)
{
    if (lev == 0) {
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > physbc
            (geom[lev], get_velocity_bcrec(),
             IncfloVelFill{m_probtype, m_bc_velocity});
        FillPatchSingleLevel(vel, IntVect(ng), time,
                             {&(m_leveldata[lev]->velocity_o),
                              &(m_leveldata[lev]->velocity)},
                             {m_t_old[lev], m_t_new[lev]}, 0, 0, AMREX_SPACEDIM, geom[lev],
                             physbc, 0);
    } else {
        const auto& bcrec = get_velocity_bcrec();
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > cphysbc
            (geom[lev-1], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
        PhysBCFunct<GpuBndryFuncFab<IncfloVelFill> > fphysbc
            (geom[lev], bcrec, IncfloVelFill{m_probtype, m_bc_velocity});
#ifdef AMREX_USE_EB
        Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
            (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
        Interpolater* mapper = &cell_cons_interp;
#endif
        FillPatchTwoLevels(vel, IntVect(ng), time,
                           {&(m_leveldata[lev-1]->velocity_o),
                            &(m_leveldata[lev-1]->velocity)},
                           {m_t_old[lev-1], m_t_new[lev-1]},
                           {&(m_leveldata[lev]->velocity_o),
                            &(m_leveldata[lev]->velocity)},
                           {m_t_old[lev], m_t_new[lev]},
                           0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                           cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), mapper, bcrec, 0);
    }
}

*/

}
