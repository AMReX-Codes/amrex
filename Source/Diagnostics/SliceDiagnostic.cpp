#include "SliceDiagnostic.H"
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil_F.H>

#include <WarpX.H>

using namespace amrex;
std::unique_ptr<MultiFab> CreateSlice( const amrex::MultiFab *mf, const amrex::Vector<Geometry> dom_geom, amrex::RealBox &slice_realbox, amrex::IntVect &slice_cr_ratio )
{
    std::unique_ptr<MultiFab> smf;
    std::unique_ptr<MultiFab> cs_mf;

    Vector<int> slice_ncells(AMREX_SPACEDIM);
    int nghost = 1;
    int nlevels = dom_geom.size();
    int ncomp = (*mf).nComp();
    
    IntVect cr_ratio(AMREX_D_DECL(0,0,0));
    const auto conversionType = (*mf).ixType();
    IntVect SliceType(AMREX_D_DECL(0,0,0));
 
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim )
    {
        if ( conversionType.nodeCentered(idim))
        {
           SliceType[idim] = 1;
        }
    }

    const RealBox& real_box = dom_geom[0].ProbDomain();
    RealBox slice_cc_nd_box;
    int slice_grid_size = 32;
    int max_ratio = 1;
         
    bool interpolate = false;
    bool coarsen = false;

    // same index space as domain //
    IntVect slice_lo(AMREX_D_DECL(0,0,0));
    IntVect slice_hi(AMREX_D_DECL(1,1,1));
    IntVect slice_lo2(AMREX_D_DECL(0,0,0));
    IntVect interp_lo(AMREX_D_DECL(0,0,0));

    // If inheriting data type //
    CheckSliceInput(real_box, slice_cc_nd_box, slice_realbox, cr_ratio,                                          slice_cr_ratio, dom_geom, SliceType, slice_lo, slice_hi, interp_lo);

    // Determine if interpolation is required and number of cells in slice //
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {    
       // Flag for interpolation if required //
       if ( interp_lo[idim] == 1) {
          interpolate = 1;
       }

       // For the case when a dimension is reduced //
       if ( ( slice_hi[idim] - slice_lo[idim]) == 1) {       
          slice_ncells[idim] = 1;
       }
       else {       
          slice_ncells[idim] = ( slice_hi[idim] - slice_lo[idim]                                                           + 1 )/cr_ratio[idim];

          int refined_ncells = (slice_hi[idim] - slice_lo[idim]) + 1 ;
          if ( cr_ratio[idim] > 1) {          
             coarsen = true;
             if (max_ratio < cr_ratio[idim] ) {
                if (slice_grid_size < cr_ratio[idim] ) {                
                    slice_grid_size  = cr_ratio[idim] ;
                }
                max_ratio = cr_ratio[idim];
             }

             // modify slice_grid_size if >= refines_cells //
             if ( slice_grid_size >= refined_ncells ) {
                slice_grid_size = refined_ncells - 1;
             }

          }
       }
    }

    // Slice generation with index type inheritance //
    Box slice(slice_lo, slice_hi);
   
    // Convert from cc to index type of parent multifab //
    bool slicetypeToBeConverted = 0;
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim ) {    
       if ( SliceType[idim] == 1 ) {
          slicetypeToBeConverted = 1;
       }
    }

    Vector<BoxArray> sba(1);
    sba[0].define(slice);
    sba[0].maxSize(slice_grid_size);

    // Distribution mapping for slice can be different from that of domain
    Vector<DistributionMapping> sdmap(1);
    sdmap[0] = DistributionMapping{sba[0]};
    
    if ( slicetypeToBeConverted==1 ) {
       smf.reset(new MultiFab(amrex::convert(sba[0],SliceType),sdmap[0],ncomp,nghost));
    }
    else {
       smf.reset(new MultiFab(sba[0],sdmap[0],ncomp,nghost));
    }

    // Copy data from domain to slice that has same cell size as that of the domain mf.  
    // src and dst have the same number of ghost cells //
    smf->ParallelCopy(*mf, 0, 0, ncomp,nghost,nghost);

    // inteprolate if required on refined slice //
    if (interpolate == 1 ) {    
       InterpolateSliceValues( *smf, interp_lo, slice_cc_nd_box, dom_geom, ncomp,                                            nghost, slice_lo, slice_hi, SliceType, real_box);
    }

    if ( coarsen == false ) {
       amrex::Print() << " Cell sizes are equal. No averaging required for slice data " << "\n";
       return smf;
    }
    else if ( coarsen == true ) {
       Vector<BoxArray> crse_ba(1);
       crse_ba[0] = sba[0];
       crse_ba[0].coarsen(cr_ratio);

       AMREX_ALWAYS_ASSERT(crse_ba[0].size() == sba[0].size());

       if(slicetypeToBeConverted==1) {
          cs_mf.reset(new MultiFab(amrex::convert(crse_ba[0],SliceType), sdmap[0],                                ncomp,nghost));
       }
       else {
          cs_mf.reset(new MultiFab(crse_ba[0], sdmap[0], ncomp, nghost));
       }

       MultiFab& mfSrc = *smf;
       MultiFab& mfDst = *cs_mf;

       MFIter mfi_dst(mfDst);
       for (MFIter mfi(mfSrc); mfi.isValid(); ++mfi) {       

           FArrayBox& Src_fabox = mfSrc[mfi];

           const Box& Dst_bx = mfi_dst.validbox();
           FArrayBox& Dst_fabox = mfDst[mfi_dst];

           int scomp = 0;
           int dcomp = 0;

           amrex::IntVect cctype(AMREX_D_DECL(0,0,0));
           if( SliceType==cctype ) {
              amrex::amrex_avgdown(Dst_bx, Dst_fabox, Src_fabox, dcomp, scomp,                                                 ncomp, cr_ratio);
           }
           amrex::IntVect ndtype(AMREX_D_DECL(1,1,1));
           if( SliceType == ndtype ) {
              amrex::amrex_avgdown_nodes(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio);
           }
           if( SliceType == WarpX::Ex_nodal_flag  ) {
              amrex::amrex_avgdown_edges(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 0);
           }
           if( SliceType == WarpX::Ey_nodal_flag) {
              amrex::amrex_avgdown_edges(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 1);
           }
           if( SliceType == WarpX::Ez_nodal_flag ) {
              amrex::amrex_avgdown_edges(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 2);
           }
           if( SliceType == WarpX::Bx_nodal_flag) {
              amrex::amrex_avgdown_faces(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 0);
           }
           if( SliceType == WarpX::By_nodal_flag ) {
              amrex::amrex_avgdown_faces(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 1);
           }
           if( SliceType == WarpX::Bz_nodal_flag ) {
              amrex::amrex_avgdown_faces(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 2);
           }

           if ( mfi_dst.isValid() ) {           
              ++mfi_dst;
           }

       }
       return cs_mf;

    }
    
    return smf;
}



void CheckSliceInput(const RealBox real_box, RealBox &slice_cc_nd_box, RealBox &slice_realbox, IntVect &cr_ratio, IntVect slice_cr_ratio, Vector<Geometry> dom_geom, IntVect const SliceType, IntVect &slice_lo, IntVect &slice_hi, IntVect &interp_lo)
{
    IntVect slice_lo2(AMREX_D_DECL(0,0,0));
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
        // Modify coarsening ratio if the input valie is not an exponnt of 2 for AMR //
        if ( slice_cr_ratio[idim] > 0 ) {        
            int log_cr_ratio = floor ( log2( slice_cr_ratio[idim]));
            slice_cr_ratio[idim] = exp2( log_cr_ratio );
        }
   
        // Default coarsening ratio i 1 //
        if ( slice_cr_ratio[idim] > 0 ) cr_ratio[idim] = slice_cr_ratio[idim];

        // Modify lo if input is out of bounds //
        if ( slice_realbox.lo(idim) < real_box.lo(idim) ) {        
            slice_realbox.setLo( idim, real_box.lo(idim));
            amrex::Print() << " slice lo is out of bounds. Modified it in dimension " << idim << " to be aligned with the domain box\n";
        }      
        
        // Modify hi if input in out od bounds //
        if ( slice_realbox.hi(idim) > real_box.hi(idim) ) {        
            slice_realbox.setHi( idim, real_box.hi(idim));
            amrex::Print() << " slice hi is out of bounds. Modified it in dimension " << idim << " to be aligned with the domain box\n";
        }
   
        // Factor to ensure index values computation depending on index type //
        double fac = ( 1.0 - SliceType[idim] ) * dom_geom[0].CellSize(idim) * 0.5;
 
        // if dimension is reduced to one cell length //
        if ( slice_realbox.hi(idim) - slice_realbox.lo(idim) <= 0)
        {
            slice_cc_nd_box.setLo( idim, slice_realbox.lo(idim) );
            slice_cc_nd_box.setHi( idim, slice_realbox.hi(idim) );

            if ( slice_cr_ratio[idim] > 1) cr_ratio[idim] = 1;
    
            // check for interpolation -- compute index lo with floor and ceil
            if ( slice_cc_nd_box.lo(idim) - real_box.lo(idim) >= fac ) {            
                slice_lo[idim] = floor( ( (slice_cc_nd_box.lo(idim) - (real_box.lo(idim)                                          + fac ) ) / dom_geom[0].CellSize(idim)) + fac * 1E-10);
                slice_lo2[idim] = ceil( ( (slice_cc_nd_box.lo(idim) - (real_box.lo(idim)                                          + fac) ) / dom_geom[0].CellSize(idim)) - fac * 1E-10 );    
            }            
            else {            
                slice_lo[idim] =  round((slice_cc_nd_box.lo(idim) - (real_box.lo(idim) ) )                                    / dom_geom[0].CellSize(idim));
                slice_lo2[idim] =  ceil((slice_cc_nd_box.lo(idim) - (real_box.lo(idim) ) )                                    / dom_geom[0].CellSize(idim) );
            }
 
            // flag for interpolation -- if reduced dimension location does not align with data point 
            if ( slice_lo[idim] == slice_lo2[idim]) {            
               if ( slice_cc_nd_box.lo(idim) - real_box.lo(idim) < fac ) {               
                  interp_lo[idim] = 1;
               }
            }
            else {            
               interp_lo[idim] = 1;
            }

            // ncells = 1 if dimension is reduced //
            slice_hi[idim] = slice_lo[idim] + 1;
        }
        else
        {
            // moving realbox.lo and reabox.hi to nearest coarsenable grid point //
            int index_lo = floor( ( (slice_realbox.lo(idim) - (real_box.lo(idim)) )                                         / dom_geom[0].CellSize(idim) ) +fac * 1E-10);
            int index_hi = ceil( ( (slice_realbox.hi(idim) - (real_box.lo(idim)) )                                         / dom_geom[0].CellSize(idim) ) - fac * 1E-10);
            bool modify_cr = true;

            while ( modify_cr == true) {            
                int lo_new = index_lo;
                int hi_new = index_hi;
                int mod_lo = index_lo % cr_ratio[idim];
                int mod_hi = index_hi % cr_ratio[idim];
                modify_cr = false;

                // To ensure that the index.lo is coarsenable //
                if ( mod_lo > 0) {
                   lo_new = index_lo - mod_lo;
                }
                // To ensure that the index.hi is coarsenable //
                if ( mod_hi > 0) {
                   hi_new = index_hi + (cr_ratio[idim] - mod_hi);
                }

            //    //If modified index.hi is > baselinebox.hi, reduce coarsening ratio,                          and provide more points that asked for //               
                if ( (hi_new * dom_geom[0].CellSize(idim)) > real_box.hi(idim) - real_box.lo(idim) )
                {
                   cr_ratio[idim] = cr_ratio[idim]/2;
                   modify_cr = true;
                }

                int ncells = (hi_new - lo_new);
                  
                // If refined cells is not an integer multiple of coarsening ratio,                            then reduce coarsening ratio by factor of 2 // 

                if ( ( ncells % cr_ratio[idim] ) != 0 ) {                
                    cr_ratio[idim] = cr_ratio[idim]/2;
                    modify_cr = true;
                }

                if ( modify_cr == false ) {                
                   index_lo = lo_new;
                   index_hi = hi_new;
                }

                slice_lo[idim] = index_lo;
                slice_hi[idim] = index_hi - 1; // since default is cell-centered    
            }
            slice_realbox.setLo( idim, index_lo * dom_geom[0].CellSize(idim) + real_box.lo(idim) );
            slice_realbox.setHi( idim, index_hi * dom_geom[0].CellSize(idim) + real_box.lo(idim) );
            slice_cc_nd_box.setLo( idim, slice_realbox.lo(idim) + fac );
            slice_cc_nd_box.setHi( idim, slice_realbox.hi(idim) - fac );
        }
    }

}



//void InterpolateSliceValues( Vector<MultiFab*>smf, IntVect interp_lo, RealBox slice_realbox, Vector<Geometry> geom, int ncomp, int nghost, IntVect slice_lo, IntVect slice_hi, IntVect SliceType, const RealBox real_box)
void InterpolateSliceValues( MultiFab& smf, IntVect interp_lo, RealBox slice_realbox, Vector<Geometry> geom, int ncomp, int nghost, IntVect slice_lo, IntVect slice_hi, IntVect SliceType, const RealBox real_box)
{
    for (MFIter mfi(smf); mfi.isValid(); ++mfi)
    {
         //MultiFab &mfSrc= *smf[0];
//         MultiFab &mfSrc= smf;
         const Box& bx = mfi.tilebox();
         const auto IndType = smf.ixType();
         const auto lo = amrex::lbound(bx);
         const auto hi = amrex::ubound(bx);
         FArrayBox& fabox = smf[mfi];

         for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {         
             if ( interp_lo[idim] == 1 ) {             
                InterpolateLo( bx, fabox, slice_lo, geom, idim, SliceType,                                                 slice_realbox, 0, ncomp, nghost, real_box);
             }
         }
    }

}

void InterpolateLo(const Box& bx, FArrayBox &fabox, IntVect slice_lo, Vector<Geometry> geom, int idir, IntVect IndType, RealBox slice_realbox, int srccomp, int ncomp, int nghost, const RealBox real_box )
{
    auto fabarr = fabox.array();
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    double fac = ( 1.0-IndType[idir] )*geom[0].CellSize(idir) * 0.5;
    int imin = slice_lo[idir];
    double minpos = imin*geom[0].CellSize(idir) + fac + real_box.lo(idir);
    double maxpos = (imin+1)*geom[0].CellSize(idir) + fac + real_box.lo(idir);
    double slice_minpos = slice_realbox.lo(idir) ;

    switch (idir) {
    case 0:
    {
        if ( imin >= lo.x && imin <= lo.x) {
           for (int n = srccomp; n < srccomp + ncomp; ++n) {
              for (int k = lo.z; k <= hi.z; ++k) {
                 for (int j = lo.y; j <= hi.y; ++j) {
                     for (int i = lo.x; i <= hi.x; ++i) {
                           double minval = fabarr(i,j,k,n);
                           double maxval = fabarr(i+1,j,k,n);
                           double ratio  = (maxval - minval) / (maxpos - minpos);
                           double xdiff  = slice_minpos - minpos;
                           double newval = minval + xdiff * ratio;
                           fabarr(i,j,k,n) = newval;
                     }
                 }
              }
           }
        }
        break;
    }
    case 1:
    {   
        if ( imin >= lo.y && imin <= lo.y) {
           for (int n = srccomp; n < srccomp+ncomp; ++n) {
              for (int k = lo.z; k <= hi.z; ++k) {
                 for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        double minval = fabarr(i,j,k,n);
                        double maxval = fabarr(i,j+1,k,n);
                        double ratio  = (maxval - minval) / (maxpos - minpos);
                        double xdiff  = slice_minpos - minpos;
                        double newval = minval + xdiff * ratio;
                        fabarr(i,j,k,n) = newval;
                    }
                 }
              }
           }
        }
        break;
    }
    case 2:
    {
        if ( imin >= lo.z && imin <= lo.z) {
           for (int n = srccomp; n < srccomp+ncomp; ++n) {
              for (int k = lo.z; k <= hi.z; ++k) {
                 for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        double minval = fabarr(i,j,k,n);
                        double maxval = fabarr(i,j,k+1,n);
                        double ratio  = (maxval - minval) / (maxpos - minpos);
                        double xdiff  = slice_minpos - minpos;
                        double newval = minval + xdiff * ratio;
                        fabarr(i,j,k,n) = newval;
                    }
                 }
              }
           }
        }
        break;
    }

    }

}







