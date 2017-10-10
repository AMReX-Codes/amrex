
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_F.H>
#include <sstream>
#include <iostream>

namespace amrex
{
    void average_node_to_cellcenter (MultiFab& cc, int dcomp, const MultiFab& nd, int scomp, int ncomp)
    {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
            const Box& bx = mfi.tilebox();
            amrex_fort_avg_nd_to_cc(bx.loVect(), bx.hiVect(), &ncomp,
                                    BL_TO_FORTRAN_N(cc[mfi],dcomp),
                                    BL_TO_FORTRAN_N(nd[mfi],scomp));
        }
    }

    void average_edge_to_cellcenter (MultiFab& cc, int dcomp, const Vector<const MultiFab*>& edge)
    {
	BL_ASSERT(cc.nComp() >= dcomp + BL_SPACEDIM);
	BL_ASSERT(edge.size() == BL_SPACEDIM);
	BL_ASSERT(edge[0]->nComp() == 1);
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVG_EG_TO_CC,bl_avg_eg_to_cc)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN_N(cc[mfi],dcomp),
		 AMREX_D_DECL(BL_TO_FORTRAN((*edge[0])[mfi]),
			BL_TO_FORTRAN((*edge[1])[mfi]),
			BL_TO_FORTRAN((*edge[2])[mfi])));
	}	
    }

    void average_face_to_cellcenter (MultiFab& cc, int dcomp, const Vector<const MultiFab*>& fc)
    {
	BL_ASSERT(cc.nComp() >= dcomp + BL_SPACEDIM);
	BL_ASSERT(fc.size() == BL_SPACEDIM);
	BL_ASSERT(fc[0]->nComp() == 1);

	Real dx[3] = {1.0,1.0,1.0};
	Real problo[3] = {0.,0.,0.};
	int coord_type = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVG_FC_TO_CC,bl_avg_fc_to_cc)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN_N(cc[mfi],dcomp),
		 AMREX_D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 dx, problo, coord_type);
	}
    }

    void average_face_to_cellcenter (MultiFab& cc, const Vector<const MultiFab*>& fc,
				     const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() >= BL_SPACEDIM);
	BL_ASSERT(fc.size() == BL_SPACEDIM);
	BL_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

	const Real* dx     = geom.CellSize();
	const Real* problo = geom.ProbLo();
	int coord_type = Geometry::Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& bx = mfi.tilebox();

	    BL_FORT_PROC_CALL(BL_AVG_FC_TO_CC,bl_avg_fc_to_cc)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(cc[mfi]),
		 AMREX_D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 dx, problo, coord_type);
	}
    }

    void average_cellcenter_to_face (const Vector<MultiFab*>& fc, const MultiFab& cc,
				     const Geometry& geom)
    {
	BL_ASSERT(cc.nComp() == 1);
	BL_ASSERT(cc.nGrow() >= 1);
	BL_ASSERT(fc.size() == BL_SPACEDIM);
	BL_ASSERT(fc[0]->nComp() == 1); // We only expect fc to have the gradient perpendicular to the face

	const Real* dx     = geom.CellSize();
	const Real* problo = geom.ProbLo();
	int coord_type = Geometry::Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
	{
	    const Box& xbx = mfi.nodaltilebox(0);
#if (BL_SPACEDIM > 1)
	    const Box& ybx = mfi.nodaltilebox(1);
#endif
#if (BL_SPACEDIM == 3)
	    const Box& zbx = mfi.nodaltilebox(2);
#endif
	    
	    BL_FORT_PROC_CALL(BL_AVG_CC_TO_FC,bl_avg_cc_to_fc)
		(xbx.loVect(), xbx.hiVect(),
#if (BL_SPACEDIM > 1)
		 ybx.loVect(), ybx.hiVect(),
#endif
#if (BL_SPACEDIM == 3)
		 zbx.loVect(), zbx.hiVect(),
#endif
		 AMREX_D_DECL(BL_TO_FORTRAN((*fc[0])[mfi]),
			BL_TO_FORTRAN((*fc[1])[mfi]),
			BL_TO_FORTRAN((*fc[2])[mfi])),
		 BL_TO_FORTRAN(cc[mfi]),
		 dx, problo, coord_type);
	}
    }

// *************************************************************************************************************

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab.
    // We do NOT assume that the coarse layout is a coarsened version of the fine layout.
    // This version DOES use volume-weighting.
    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
		       const Geometry& fgeom, const Geometry& cgeom, 
                       int scomp, int ncomp, int rr)
    {
         average_down(S_fine,S_crse,fgeom,cgeom,scomp,ncomp,rr*IntVect::TheUnitVector());
    }

    void average_down (const MultiFab& S_fine, MultiFab& S_crse,
		       const Geometry& fgeom, const Geometry& cgeom, 
                       int scomp, int ncomp, const IntVect& ratio)
    {
  
        if (S_fine.is_nodal() || S_crse.is_nodal())
        {
            amrex::Error("Can't use amrex::average_down for nodal MultiFab!");
        }

#if (BL_SPACEDIM == 3)
	amrex::average_down(S_fine, S_crse, scomp, ncomp, ratio);
	return;
#else

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
	const BoxArray& fine_BA = S_fine.boxArray();
	const DistributionMapping& fine_dm = S_fine.DistributionMap();
        BoxArray crse_S_fine_BA = fine_BA; 
	crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA,fine_dm,ncomp,0);

	MultiFab fvolume;
	fgeom.GetVolume(fvolume, fine_BA, fine_dm, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& tbx = mfi.tilebox();

            //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
            //        because the crse fab is a temporary which was made starting at comp 0, it is
            //        not part of the actual crse multifab which came in.
            BL_FORT_PROC_CALL(BL_AVGDOWN_WITH_VOL,bl_avgdown_with_vol)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi],0),
                 BL_TO_FORTRAN(fvolume[mfi]),
                 ratio.getVect(),&ncomp);
	}

        S_crse.copy(crse_S_fine,0,scomp,ncomp);
#endif
   }

// *************************************************************************************************************

    // Average fine cell-based MultiFab onto crse cell-centered MultiFab.
    // We do NOT assume that the coarse layout is a coarsened version of the fine layout.
    // This version does NOT use volume-weighting
    void average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp, int rr)
    {
         average_down(S_fine,S_crse,scomp,ncomp,rr*IntVect::TheUnitVector());
    }


    void sum_fine_to_coarse(const MultiFab& S_fine, MultiFab& S_crse,
                            int scomp, int ncomp, const IntVect& ratio,
                            const Geometry& cgeom, const Geometry& fgeom) 
    {
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());
        BL_ASSERT(ratio == ratio[0]);
        BL_ASSERT(S_fine.nGrow() % ratio[0] == 0);

        const int nGrow = S_fine.nGrow() / ratio[0];

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, nGrow);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse_S_fine, true); mfi.isValid(); ++mfi)
        {
            //  NOTE: The tilebox is defined at the coarse level.
            const Box& tbx = mfi.growntilebox(nGrow);
            
            BL_FORT_PROC_CALL(BL_AVGDOWN, bl_avgdown)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(S_fine[mfi] , scomp),
                 BL_TO_FORTRAN_N(crse_S_fine[mfi], 0),
                 ratio.getVect(),&ncomp);
        }
        
        S_crse.copy(crse_S_fine, 0, scomp, ncomp, nGrow, 0,
                    cgeom.periodicity(), FabArrayBase::ADD);
    }

    void average_down (const MultiFab& S_fine, MultiFab& S_crse, 
                       int scomp, int ncomp, const IntVect& ratio)
    {
        BL_ASSERT(S_crse.nComp() == S_fine.nComp());

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);
        
        if (crse_S_fine_BA == S_crse.boxArray() and S_fine.DistributionMap() == S_crse.DistributionMap())
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(S_crse,true); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& tbx = mfi.tilebox();
                
                BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                    (tbx.loVect(), tbx.hiVect(),
                     BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                     BL_TO_FORTRAN_N(S_crse[mfi],scomp),
                     ratio.getVect(),&ncomp);
            }
        }
        else
        {
            MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp,0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
            {
                //  NOTE: The tilebox is defined at the coarse level.
                const Box& tbx = mfi.tilebox();
                
                //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
                //        because the crse fab is a temporary which was made starting at comp 0, it is
                //        not part of the actual crse multifab which came in.
                
                BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                    (tbx.loVect(), tbx.hiVect(),
                     BL_TO_FORTRAN_N(S_fine[mfi],scomp),
                     BL_TO_FORTRAN_N(crse_S_fine[mfi],0),
                     ratio.getVect(),&ncomp);
            }
            
            S_crse.copy(crse_S_fine,0,scomp,ncomp);
        }
   }

// *************************************************************************************************************

    // Average fine face-based MultiFab onto crse face-based MultiFab.
    // This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_faces (const Vector<const MultiFab*>& fine, const Vector<MultiFab*>& crse,
			     const IntVect& ratio, int ngcrse)
    {
	BL_ASSERT(crse.size()  == BL_SPACEDIM);
	BL_ASSERT(fine.size()  == BL_SPACEDIM);
	BL_ASSERT(crse[0]->nComp() == fine[0]->nComp());

	int ncomp = crse[0]->nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (int n=0; n<BL_SPACEDIM; ++n) {
            for (MFIter mfi(*crse[n],true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.growntilebox(ngcrse);

                BL_FORT_PROC_CALL(BL_AVGDOWN_FACES,bl_avgdown_faces)
                    (tbx.loVect(),tbx.hiVect(),
                     BL_TO_FORTRAN((*fine[n])[mfi]),
                     BL_TO_FORTRAN((*crse[n])[mfi]),
                     ratio.getVect(),n,ncomp);
            }
        }
    }

    //! Average fine edge-based MultiFab onto crse edge-based MultiFab.
    //! This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_edges (const Vector<const MultiFab*>& fine, const Vector<MultiFab*>& crse,
                             const IntVect& ratio, int ngcrse)
    {
	BL_ASSERT(crse.size()  == BL_SPACEDIM);
	BL_ASSERT(fine.size()  == BL_SPACEDIM);
	BL_ASSERT(crse[0]->nComp() == fine[0]->nComp());

	int ncomp = crse[0]->nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (int n=0; n<BL_SPACEDIM; ++n) {
            for (MFIter mfi(*crse[n],true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.growntilebox(ngcrse);

                BL_FORT_PROC_CALL(BL_AVGDOWN_EDGES,bl_avgdown_edges)
                    (tbx.loVect(),tbx.hiVect(),
                     BL_TO_FORTRAN((*fine[n])[mfi]),
                     BL_TO_FORTRAN((*crse[n])[mfi]),
                     ratio.getVect(),n,ncomp);
            }
        }
    }

    //! Average fine node-based MultiFab onto crse node-based MultiFab.
    //! This routine assumes that the crse BoxArray is a coarsened version of the fine BoxArray.
    void average_down_nodal (const MultiFab& fine, MultiFab& crse, const IntVect& ratio, int ngcrse)
    {
        BL_ASSERT(fine.is_nodal());
        BL_ASSERT(crse.is_nodal());
	BL_ASSERT(crse.nComp() == fine.nComp());

	int ncomp = crse.nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(crse,true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.growntilebox(ngcrse);
                
                BL_FORT_PROC_CALL(BL_AVGDOWN_NODES,bl_avgdown_nodes)
                    (tbx.loVect(),tbx.hiVect(),
                     BL_TO_FORTRAN(fine[mfi]),
                     BL_TO_FORTRAN(crse[mfi]),
                     ratio.getVect(),&ncomp);
            }
    }

// *************************************************************************************************************

    void fill_boundary(MultiFab& mf, const Geometry& geom, bool cross)
    {
	amrex::fill_boundary(mf, 0, mf.nComp(), geom, cross);
    }

// *************************************************************************************************************

    void fill_boundary(MultiFab& mf, int scomp, int ncomp, const Geometry& geom, bool cross)
    {
	mf.FillBoundary(scomp,ncomp,geom.periodicity(),cross);
    }


    void print_state(const MultiFab& mf, const IntVect& cell, const int n)
    {


#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mf); mfi.isValid(); ++mfi)
	{
	  if (mfi.validbox().contains(cell)) {
	    if (n >= 0) {
	      amrex::AllPrint().SetPrecision(17) << mf[mfi](cell, n) << std::endl;
	    } else {
	      std::ostringstream ss;
	      ss.precision(17);
	      const int ncomp = mf.nComp();
	      for (int i = 0; i < ncomp-1; ++i)
		{
		  ss << mf[mfi](cell,i) << ", ";
		}
	      ss << mf[mfi](cell,ncomp-1);
	      amrex::AllPrint() << ss.str() << std::endl;	    
	    }
	  }
	}
      
    }

    MultiFab ToMultiFab (const iMultiFab& imf)
    {
        MultiFab mf(imf.boxArray(), imf.DistributionMap(), imf.nComp(), imf.nGrow());

        const int ncomp = imf.nComp();
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.growntilebox();
            amrex_fort_int_to_real(BL_TO_FORTRAN_BOX(tbx), &ncomp,
                                   BL_TO_FORTRAN_ANYD(mf[mfi]),
                                   BL_TO_FORTRAN_ANYD(imf[mfi]));
        }

        return mf;
    }
}
