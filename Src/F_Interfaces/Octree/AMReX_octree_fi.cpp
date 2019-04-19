
#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FAmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_C.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

namespace
{
    struct treenode {
        int level, grid;
    };
}

extern "C" {

    void amrex_fi_init_octree ()
    {
        ParmParse pp("amr");
        int cnt = pp.countval("max_grid_size");
        int max_grid_size;
        if (cnt == 0) {
            max_grid_size = 8;
            pp.add("max_grid_size", max_grid_size);
        } else if (cnt == 1) {
            pp.get("max_grid_size", max_grid_size);
        } else {
            amrex::Abort("amrex_fi_init_octree: must use the same max_grid_size for all levels");
        }

        int blocking_factor = 2*max_grid_size;
        pp.add("blocking_factor", blocking_factor);

        int max_grid_size_x = max_grid_size;
        pp.query("max_grid_size_x", max_grid_size_x);
        int blocking_factor_x = 2*max_grid_size_x;
        pp.add("blocking_factor_x", blocking_factor_x);

        int max_grid_size_y = max_grid_size;
        pp.query("max_grid_size_y", max_grid_size_y);
        int blocking_factor_y = 2*max_grid_size_y;
        pp.add("blocking_factor_y", blocking_factor_y);

        int max_grid_size_z = max_grid_size;
        pp.query("max_grid_size_z", max_grid_size_z);
        int blocking_factor_z = 2*max_grid_size_z;
        pp.add("blocking_factor_z", blocking_factor_z);

        pp.add("grid_eff", 1.0);

        int max_level;
        pp.get("max_level", max_level);

        Vector<int> ref_ratio(max_level, 2);
        pp.addarr("ref_ratio", ref_ratio);
    }

    void amrex_fi_build_octree_leaves (AmrCore* const amrcore, int* n, Vector<treenode>*& leaves)
    {
        leaves = new Vector<treenode>;
        const int finest_level = amrcore->finestLevel();
        const int myproc = ParallelDescriptor::MyProc();

        FAmrCore* famrcore = dynamic_cast<FAmrCore*>(amrcore);

        famrcore->octree_leaf_grids.resize(finest_level+1);
        famrcore->octree_leaf_dmap.resize(finest_level+1);
        famrcore->octree_leaf_dummy_mf.resize(finest_level+1);
        famrcore->octree_li_full_to_leaf.resize(finest_level+1);
        famrcore->octree_li_leaf_to_full.resize(finest_level+1);

        BoxList bl;
        Vector<int> iproc;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            famrcore->octree_li_full_to_leaf[lev].clear();
            famrcore->octree_li_leaf_to_full[lev].clear();

            const BoxArray& ba = amrcore->boxArray(lev);
            const DistributionMapping& dm = amrcore->DistributionMap(lev);
            const int ngrids = ba.size();
            BL_ASSERT(ba.size() < std::numeric_limits<int>::max());

            if (lev == finest_level)
            {
                famrcore->octree_leaf_grids[lev] = ba;
                famrcore->octree_leaf_dmap[lev] = dm;
                famrcore->octree_leaf_dummy_mf[lev].reset
                    (new MultiFab(ba,dm,1,0,MFInfo().SetAlloc(false)));

                int ilocal = 0;
                for (int i = 0; i < ngrids; ++i) {
                    if (dm[i] == myproc) {
                        leaves->push_back({lev, i});
                        famrcore->octree_li_full_to_leaf[lev].push_back(ilocal++);
                    }
                }
                famrcore->octree_li_leaf_to_full[lev] = famrcore->octree_li_full_to_leaf[lev];
            }
            else
            {
                const BoxArray& fba = amrcore->boxArray(lev+1);
                const IntVect& rr = amrcore->refRatio(lev);
                bl.clear();
                iproc.clear();
                int ilocal_full = 0;
                int ilocal_leaf = 0;
                for (int i = 0; i < ngrids; ++i) {
                    Box bx = ba[i];
                    Box fbx = amrex::refine(bx,rr);
                    if (fba.intersects(fbx)) { // non-leaf
                        if (dm[i] == myproc) {
                            famrcore->octree_li_full_to_leaf[lev].push_back(-1);
                            ++ilocal_full;
                        }
                    } else { // leaves
                        bl.push_back(bx);
                        iproc.push_back(dm[i]);
                        if (dm[i] == myproc) {
                            leaves->push_back({lev, i});
                            famrcore->octree_li_full_to_leaf[lev].push_back(ilocal_leaf++);
                            famrcore->octree_li_leaf_to_full[lev].push_back(ilocal_full++);
                        }
                    }
                }

                if (bl.size() == 0) {
                    famrcore->octree_leaf_grids[lev] = BoxArray();
                    famrcore->octree_leaf_dmap[lev] = DistributionMapping();
                    famrcore->octree_leaf_dummy_mf[lev].reset(new MultiFab());
                } else {
                    bool update_dummy_mf = false;
                    if (famrcore->octree_leaf_grids[lev] != bl.data()) {
                        famrcore->octree_leaf_grids[lev] = BoxArray(bl);
                        update_dummy_mf = true;
                    }
                    if (famrcore->octree_leaf_dmap[lev].ProcessorMap() != iproc) {
                        famrcore->octree_leaf_dmap[lev] = DistributionMapping(iproc);
                        update_dummy_mf = true;
                    }
                    famrcore->octree_leaf_dummy_mf[lev].reset
                        (new MultiFab(famrcore->octree_leaf_grids[lev],
                                      famrcore->octree_leaf_dmap[lev],
                                      1,0,MFInfo().SetAlloc(false)));
                }
            }
        }
        *n = leaves->size();
    }

    void amrex_fi_copy_octree_leaves (Vector<treenode>* leaves, treenode a_copy[])
    {
        const int n = leaves->size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; ++i) {
            a_copy[i] = (*leaves)[i];
        }
        
        delete leaves;
    }

    void amrex_fi_octree_average_down_level (AmrCore* const amrcore, int flev,
                                             MultiFab const * const fine,
                                             MultiFab * const crse,
                                             int scomp, int ncomp)
    {
        FAmrCore* famrcore = dynamic_cast<FAmrCore*>(amrcore);
        const BoxArray& lba = famrcore->octree_leaf_grids[flev];
        const DistributionMapping& ldm = famrcore->octree_leaf_dmap[flev];
        const Vector<int>& li_leaf_to_full = famrcore->octree_li_leaf_to_full[flev];
        const IntVect& rr = famrcore->refRatio(flev-1);

        MultiFab lmf(amrex::coarsen(lba,rr),ldm,ncomp,0);

        MultiFab fvolume;
        if (!Geometry::IsCartesian())
        {
            const Geometry& fgeom = famrcore->Geom(flev);
            fgeom.GetVolume(fvolume, lba, ldm, 0);
        }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(lmf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            FArrayBox * crsefab = lmf.fabPtr(mfi);
            const int li = li_leaf_to_full[mfi.LocalIndex()];
            FArrayBox const* finefab = fine->fabPtrAtLocalIdx(li);
            if (Geometry::IsCartesian()) {
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown(tbx,*crsefab,*finefab,0,scomp,ncomp,rr);
                });
            } else {
                FArrayBox const* finevolfab = fvolume.fabPtr(mfi);
                AMREX_LAUNCH_HOST_DEVICE_LAMBDA ( bx, tbx,
                {
                    amrex_avgdown_with_vol(tbx,*crsefab,*finefab,*finevolfab,
                                           0,scomp,ncomp,rr);
                });
            }
        }

        crse->ParallelCopy(lmf,0,scomp,ncomp);
    }
}
