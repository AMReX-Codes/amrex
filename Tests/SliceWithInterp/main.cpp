#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_VisMF.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Vector.H>
#include <main.H>
#include <AMReX_BLassert.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFabUtil_C.H>


using namespace amrex;

void InitializeVariables(Vector<MultiFab*> mf, Vector<Geometry> geom, int srccomp, int ncomp);

void CreateSlice(Vector<MultiFab*> mf, Vector<Geometry> geom, RealBox slice_realbox, Vector<int> slice_cr_ratio);

void CheckSliceInput(const RealBox real_box, RealBox &slice_cc_nd_box, RealBox &slice_realbox,IntVect &cr_ratio, Vector<int> slice_cr_ratio,Vector<Geometry> const geom, IntVect const SliceType, IntVect& slice_lo, IntVect &slice_hi, IntVect &interp_lo);

void InterpolateLo(const Box& bx, FArrayBox &fabox, IntVect slice_lo, Vector<Geometry> geom, int idir, IntVect IndType, RealBox slice_realbox, int srccomp, int ncomp, int nghost);

void InterpolateSliceValues( Vector<MultiFab*>smf, IntVect interp_lo, RealBox slice_realbox, Vector<Geometry> geom, int ncomp, int nghost, IntVect slice_lo, IntVect slice_hi, IntVect ConvertTypeofSlice);


int main(int argc, char* argv[])
{    
    amrex::Initialize(argc,argv);
    { 
    const int nghost = 1;
    int ncells, max_grid_size, ncomp, nlevs, nppc;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("ncomp", ncomp);
    pp.get("nlevs", nlevs);
    pp.get("nppc", nppc);
    
    
    AMREX_ALWAYS_ASSERT(nlevs < 2); // relax this later
    AMREX_ALWAYS_ASSERT(nghost > 0); 

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(ncells-1, ncells-1, ncells-1)); 
    const Box domain(domain_lo, domain_hi);

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    // Define the refinement ratio
    Vector<IntVect> ref_ratio(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        ref_ratio[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object for each level
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), ref_ratio[lev-1]),
			 &real_box, CoordSys::cartesian, is_per);
    }
    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);

    // default cell-centered multifab //
    Vector<std::unique_ptr<MultiFab> > mf(nlevs);
    // node-based multifab //
    Vector<std::unique_ptr<MultiFab> > rho_mf(nlevs);
    // Ex - Edge based //
    Vector<std::unique_ptr<MultiFab> > Ex_mf(nlevs);
    // Ey - Edge based //
    Vector<std::unique_ptr<MultiFab> > Ey_mf(nlevs);
    // Ez - Edge based //
    Vector<std::unique_ptr<MultiFab> > Ez_mf(nlevs);
    // Bx - face centered //
    Vector<std::unique_ptr<MultiFab> > Bx_mf(nlevs);
    // By - face centered //
    Vector<std::unique_ptr<MultiFab> > By_mf(nlevs);
    // Bz - face centered //
    Vector<std::unique_ptr<MultiFab> > Bz_mf(nlevs);

    // note that current density is edge centered similar to E-Field and therefore for testing purposes it is not explicitly included here. 
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        mf[lev].reset(new MultiFab(ba[lev], dmap[lev], ncomp, nghost));
        rho_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{1,1,1}), dmap[lev], ncomp, nghost));
        Ex_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{0,1,1}), dmap[lev], ncomp, nghost));
        Ey_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{1,0,1}), dmap[lev], ncomp, nghost));
        Ez_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{1,1,0}), dmap[lev], ncomp, nghost));
        Bx_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{1,0,0}), dmap[lev], ncomp, nghost));
        By_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{0,1,0}), dmap[lev], ncomp, nghost));
        Bz_mf[lev].reset(new MultiFab(amrex::convert(ba[lev], IntVect{0,0,1}), dmap[lev], ncomp, nghost));
    }

    InitializeVariables( GetVecOfPtrs(mf), geom , 0, ncomp);
    InitializeVariables( GetVecOfPtrs(rho_mf), geom, 0, ncomp);
    InitializeVariables( GetVecOfPtrs(Ex_mf), geom, 0, ncomp );
    InitializeVariables( GetVecOfPtrs(Ey_mf), geom, 0, ncomp );
    InitializeVariables( GetVecOfPtrs(Ez_mf), geom, 0, ncomp );
    InitializeVariables( GetVecOfPtrs(Bx_mf), geom, 0, ncomp );
    InitializeVariables( GetVecOfPtrs(By_mf), geom, 0, ncomp );
    InitializeVariables( GetVecOfPtrs(Bz_mf), geom, 0, ncomp );

    VisMF::Write((*mf[0]),"vismf_orig_cc");
    VisMF::Write((*rho_mf[0]),"vismf_orig_node");
    VisMF::Write((*Ex_mf[0]),"vismf_orig_Ex");
    VisMF::Write((*Ey_mf[0]),"vismf_orig_Ey");
    VisMF::Write((*Ez_mf[0]),"vismf_orig_Ez");
    VisMF::Write((*Bx_mf[0]),"vismf_orig_Bx");
    VisMF::Write((*By_mf[0]),"vismf_orig_By");
    VisMF::Write((*Bz_mf[0]),"vismf_orig_Bz");

    // Slice generation starts here //
    Vector<Real> slo(AMREX_SPACEDIM);
    Vector<Real> shi(AMREX_SPACEDIM);
    Vector<int> slice_cells(AMREX_SPACEDIM);
    Vector<int> slice_crse_ratio(AMREX_SPACEDIM);
   
    //Read input for slice //
    ParmParse ppg("slice");
    ppg.queryarr("dom_lo",slo,0,AMREX_SPACEDIM);
    ppg.queryarr("dom_hi",shi,0,AMREX_SPACEDIM);
    ppg.queryarr("coarsening_ratio",slice_crse_ratio,0,AMREX_SPACEDIM);

    // Set lo and hi for slice 
    amrex::RealBox slice_realbox;
    slice_realbox.setLo(slo); 
    slice_realbox.setHi(shi); 

    CreateSlice( GetVecOfPtrs(mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(rho_mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(Ex_mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(Ey_mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(Ez_mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(Bx_mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(By_mf), geom, slice_realbox, slice_crse_ratio);
    CreateSlice( GetVecOfPtrs(Bz_mf), geom, slice_realbox, slice_crse_ratio);
    }

    amrex::Finalize();

}

//*********************************************************************************//
// This function generates 1D,2D, or 3D multifab that .lo and .hi contained        //
// within the domain. The slice is generated in two steps :                        //
// 1. Slice multifab (smf) is generated with same cell size as the parent domain   //
//    Then Parallel Copy is used to copy data from domain -> slice                 //
// 2. Then, based on the user-defined cell size for the slice, the slice           //
//    multifab (cs_mf) is coarsened and the data is averaged refined->coarse       //
//    using the in-built average functions in AMReX_MultiFabUtil_3D_C.H            //
// Note : If the user-defined cell size for the slice is same as that of the       //
//        domain, no averaging/coarsening is performed.                            //
//        If the user-defined cell size for the slice is smaller than the domain   //
//        cell size, no slice is generated, since interpolation of data            //
//        from coarse->refine is not performed in this function.                   //
// Also the slice generation currently assumes only 1 level for the amr structure  //
//*********************************************************************************//

void CreateSlice(Vector<MultiFab*> mf, Vector<Geometry> geom, RealBox slice_realbox, Vector<int> slice_cr_ratio)
{

    Vector<int> slice_ncells(AMREX_SPACEDIM);
    int nghost = 1;
    int is_per[AMREX_SPACEDIM];

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is_per[i] = 1;
  
    int nlevs = geom.size();
    int ncomp = (*mf[0]).nComp();
    int max_ratio = 1;
    bool coarsen = false; 
    IntVect cr_ratio(AMREX_D_DECL(1,1,1));
    const auto conversionType = (*mf[0]).ixType();

    // Obtain index type of source multifab //
    IntVect SliceType(AMREX_D_DECL(0,0,0));
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim ) 
    {
        if ( conversionType.nodeCentered(idim) ) 
        {
           SliceType[idim] = 1;
        }
    }

    const RealBox& real_box = geom[0].ProbDomain(); 
 
    RealBox slice_cc_nd_box;    
    // Default max_grid_size for slice //
    int slice_grid_size = 32;

    bool interpolate = false;
    // ensuring that index space for slice is same as domain // 
    IntVect slice_lo(AMREX_D_DECL(0,0,0));
    IntVect slice_hi(AMREX_D_DECL(1,1,1)); 
    IntVect slice_lo2(AMREX_D_DECL(0,0,0));
    IntVect interp_lo(AMREX_D_DECL(0,0,0)); 

    // If inheriting data type //
    CheckSliceInput(real_box, slice_cc_nd_box, slice_realbox, cr_ratio,                                         slice_cr_ratio, geom, SliceType, slice_lo, slice_hi, interp_lo); 

    // Determine if interpolatiojnj is required and number of cells in slice //
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {

        if ( interp_lo[idim] == 1) {
           interpolate = 1;
        }

       // For the case when a dimension is reduced //
       if ( ( slice_hi[idim] - slice_lo[idim]) == 1) {
          slice_ncells[idim] = 1;
       }
       else {
            
            slice_ncells[idim] = ( slice_hi[idim] - slice_lo[idim]                                                           + (1 - SliceType[idim]))/cr_ratio[idim]; 

            int refined_ncells = (slice_hi[idim] - slice_lo[idim]) + ( 1 - SliceType[idim]);
            
            if ( cr_ratio[idim] > 1) 
            {
               coarsen = true; 
               
               if (max_ratio < cr_ratio[idim] ) 
               {
                  if (slice_grid_size < cr_ratio[idim] ) 
                  {
                      slice_grid_size  = cr_ratio[idim] ;
                  }
                  max_ratio = cr_ratio[idim];
               }

               // modify slice_grid_size if >= refines_cells //
               if ( slice_grid_size >= refined_ncells ) {
                   slice_grid_size = refined_ncells;
               }
            }
       }
    }

    // Slice generation with Index Type inheritance //

    // Default index type for slice is cell-centered //
    Box slice(slice_lo, slice_hi);
    
    // Convert from cc to index type of parent multifab //
    bool slicetypeToBeConverted = 0;
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim )
    {
       if ( SliceType[idim] == 1 ) {
          slicetypeToBeConverted = 1;
       }
    }
    
    Vector<BoxArray> sba(nlevs);
    sba[0].define(slice);
    sba[0].maxSize(slice_grid_size);

    // Distribution mapping for slice can be different from that of domain
    Vector<DistributionMapping> sdmap(nlevs);
    sdmap[0] = DistributionMapping{sba[0]};

    // multifab for slice  
    Vector<std::unique_ptr<MultiFab> > smf(nlevs);    
    if ( slicetypeToBeConverted==1 ) {
       smf[0].reset(new MultiFab(amrex::convert(sba[0],SliceType),sdmap[0],ncomp,nghost));
    }
    else {
       smf[0].reset(new MultiFab(sba[0],sdmap[0],ncomp,nghost));
    }

    // Copy data from domain to slice that has same cell size as that of the domain mf.  
    (*smf[0]).ParallelCopy((*mf[0]), 0, 0, ncomp,nghost,nghost);  
    VisMF::Write((*smf[0]),"vismf_init");

    // Interpolation of data on refined slice //
    if (interpolate == 1) {
       InterpolateSliceValues( GetVecOfPtrs(smf), interp_lo, slice_cc_nd_box, geom, ncomp,                                 nghost, slice_lo, slice_hi, SliceType);
       VisMF::Write((*smf[0]),"vismf_output_interpolatedfine");
    }

    if ( coarsen == false ) {
       amrex::Print() << " Cell sizes are equal. No averaging required. " << "\n";
       VisMF::Write((*smf[0]),"vismf_output");
    }
    else if ( coarsen == true ) {
 
       amrex::Print() << " Calling in-built amrex average functions "<< cr_ratio[0] ;
       amrex::Print() << " " << cr_ratio[1] << " " << cr_ratio[2] << "\n";

       Vector<BoxArray> crse_ba(nlevs);
       crse_ba[0] = sba[0];
       crse_ba[0].coarsen(cr_ratio);

       // The input values of max_grid_size is factored by ratio in the coarsened slice //
       int cs_grid_size = double(sba[0].size())*max_ratio ; 
       crse_ba[0].maxSize(cs_grid_size);
       AMREX_ALWAYS_ASSERT(crse_ba[0].size() == sba[0].size());

       // constructing coarsened slice as per user-input if s_cells<ncells // 
       Vector<std::unique_ptr<MultiFab> > cs_mf(nlevs); 

       if(slicetypeToBeConverted==1) {
          cs_mf[0].reset(new MultiFab(amrex::convert(crse_ba[0],SliceType), sdmap[0],                                ncomp,nghost));
       }
       else {
          cs_mf[0].reset(new MultiFab(crse_ba[0], sdmap[0], ncomp, nghost));
       }

       // Assumption :: Currently works only for 1 level //
       for (int ilev = 0; ilev < nlevs; ilev++)
       {
           MultiFab& mfSrc = *smf[ilev];
           MultiFab& mfDst = *cs_mf[ilev];
          
           MFIter mfi_dst(mfDst);
           for (MFIter mfi(mfSrc); mfi.isValid(); ++mfi)
           {
               FArrayBox& Src_fabox = mfSrc[mfi];

               const Box& Dst_bx = mfi_dst.validbox();
               FArrayBox& Dst_fabox = mfDst[mfi_dst];

               int scomp = 0;
               int dcomp = 0;

               if( SliceType==cctype ) {
                  amrex::amrex_avgdown(Dst_bx, Dst_fabox, Src_fabox, dcomp, scomp,                                                 ncomp, cr_ratio);
               }
               if( SliceType == ndtype ) {
                  amrex::amrex_avgdown_nodes(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio);
               }
               if( SliceType == xetype ) {
                  amrex::amrex_avgdown_edges(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 0);
               }
               if( SliceType == yetype) {
                  amrex::amrex_avgdown_edges(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 1);
               }
               if( SliceType == zetype ) {
                  amrex::amrex_avgdown_edges(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 2);
               }
               if (AMREX_SPACEDIM==3) {
                  if( SliceType == xftype) {
                     amrex::amrex_avgdown_faces(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 0);
                  }
                  if( SliceType == yftype ) {
                     amrex::amrex_avgdown_faces(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 1);
                  }
                  if( SliceType == zftype ) {
                     amrex::amrex_avgdown_faces(Dst_bx, Dst_fabox, Src_fabox, dcomp,                                                        scomp, ncomp, cr_ratio, 2);
                  }
               }

               if ( mfi_dst.isValid() ) 
               {
                  ++mfi_dst;
               }
           }
       }
       VisMF::Write((*cs_mf[0]),"vismf_output_coarse");
    }

}

void InitializeVariables(Vector<MultiFab*> mf, Vector<Geometry> geom, int srccomp, int ncomp)
{
    int nlevs = geom.size();
    const Real *dom_dx = geom[0].CellSize();
    for (int lev = 0; lev < nlevs; ++lev) {
        for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
            
            MultiFab &mfSrc= *mf[lev];
            auto const &mf_arr = mfSrc.array(mfi);
            const Box& bx = mfi.growntilebox();
            const auto IndType = (*mf[lev]).ixType();
            const auto lo = amrex::lbound(bx); 
            const auto hi = amrex::ubound(bx); 
            for (int n = srccomp; n < srccomp + ncomp; ++n ) {
               for (int k = lo.z; k<=hi.z; ++k) {
                  for (int j = lo.y; j<=hi.y; ++j) {
                     for (int i = lo.x; i<=hi.x; ++i) {
                         int idim = n;
                         if (n >=3 )
                         {
                            idim = n % 3;
                         }
                         int fac = i;
                         if (idim == 1) fac = j;
                         if (idim == 2) fac = k;
                         mf_arr(i,j,k,n) = fac * dom_dx[idim]                                                                          + ( 1.0 - IndType[idim] ) * dom_dx[idim] *0.5;
                     }
                  }
               }
            }
        }
    }

}

void InterpolateLo(const Box& bx, FArrayBox &fabox, IntVect slice_lo, Vector<Geometry> geom, int idir, IntVect IndType, RealBox slice_realbox, int srccomp, int ncomp, int nghost)
{
    auto fabarr = fabox.array();
    const auto lo = amrex::lbound(bx); 
    const auto hi = amrex::ubound(bx);            
    double fac = ( 1.0-IndType[idir] )*geom[0].CellSize(idir) * 0.5;
    int imin = slice_lo[idir];
    double minpos = imin*geom[0].CellSize(idir) + fac ;
    double maxpos = (imin+1)*geom[0].CellSize(idir) + fac;
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


void InterpolateSliceValues( Vector<MultiFab*>smf, IntVect interp_lo, RealBox slice_realbox, Vector<Geometry> geom, int ncomp, int nghost, IntVect slice_lo, IntVect slice_hi, IntVect SliceType)
{
    for (MFIter mfi(*smf[0]); mfi.isValid(); ++mfi) 
    {
         MultiFab &mfSrc= *smf[0];
         const Box& bx = mfi.tilebox();
         const auto IndType = (*smf[0]).ixType();
         const auto lo = amrex::lbound(bx); 
         const auto hi = amrex::ubound(bx);            
         FArrayBox& fabox = mfSrc[mfi];
 
         for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
         {
             if ( interp_lo[idim] == 1 ) 
             {
                InterpolateLo( bx, fabox, slice_lo, geom, idim, SliceType,                                                 slice_realbox, 0, ncomp, nghost);
             }
         }
    }

}
    

     
void CheckSliceInput(const RealBox real_box, RealBox &slice_cc_nd_box, RealBox &slice_realbox,IntVect &cr_ratio, Vector<int> slice_cr_ratio,Vector<Geometry> const geom, IntVect const SliceType, IntVect& slice_lo, IntVect &slice_hi, IntVect &interp_lo)
{

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) 
    {
       // Modify coarsening ratio if the input value is not an exponent of 2 for AMR //
       if (slice_cr_ratio[idim] > 0 ) 
       {
           int log_cr_ratio = floor( log2(slice_cr_ratio[idim])) ;
           slice_cr_ratio[idim] = exp2( log_cr_ratio);  
       }

       // Default coarsening ratio is 1 //
       if ( slice_cr_ratio[idim] > 0 ) cr_ratio[idim] = slice_cr_ratio[idim];

       // Modify lo if input is out of bounds //
       if ( slice_realbox.lo(idim) < real_box.lo(idim) ) 
       {
            slice_realbox.setLo( idim, real_box.lo(idim) ) ;
            amrex::Print() << " Slice lo out of bounds. Modified slice lo along dim : " <<                                 idim << " to be aligned with domain box\n";
       }

       // Modify hi if input is out of bounds //
       if ( slice_realbox.hi(idim) > real_box.hi(idim) ) 
       {
            slice_realbox.setHi( idim, real_box.hi(idim) );
            amrex::Print() << " Slice hi out of bounds. Modified slice lo along dim : " <<                                 idim << " to be aligned with domain box\n";
       }

       double fac = ( 1.0-SliceType[idim] )*geom[0].CellSize(idim) * 0.5;

       // If dimension is reduced to one cell length //
       if ( slice_realbox.hi(idim) - slice_realbox.lo(idim) <= 0) 
       {
          slice_cc_nd_box.setLo( idim, slice_realbox.lo(idim) );
          slice_cc_nd_box.setHi( idim, slice_realbox.hi(idim) );

          if ( slice_cr_ratio[idim] > 1 )  cr_ratio[idim] = 1; 

          // check for interpolation -- compute index lo with floor and ceil //
          IntVect slice_lo2(AMREX_D_DECL(0,0,0));
          if ( slice_cc_nd_box.lo(idim) - real_box.lo(idim) >= fac ) 
          {
             slice_lo[idim] = floor( ( (slice_cc_nd_box.lo(idim) - (real_box.lo(idim)                                          + fac ) ) / geom[0].CellSize(idim)) + fac * 1E-10);
             slice_lo2[idim] = ceil( ( (slice_cc_nd_box.lo(idim) - (real_box.lo(idim)                                          + fac) ) / geom[0].CellSize(idim)) - fac * 1E-10 );
          }
          else  
          {
              slice_lo[idim] =  round((slice_cc_nd_box.lo(idim) - (real_box.lo(idim) ) )                                    / geom[0].CellSize(idim));
              slice_lo2[idim] =  ceil((slice_cc_nd_box.lo(idim) - (real_box.lo(idim) ) )                                    / geom[0].CellSize(idim) );
          }

          if ( slice_lo[idim] == slice_lo2[idim]) 
          {
             if ( slice_cc_nd_box.lo(idim) - real_box.lo(idim) < fac ) 
             {
                interp_lo[idim] = 1;
             }
          }
          else 
          {
             interp_lo[idim] = 1;
          }

          // ncells = 1 if dimension is reduced //
          slice_hi[idim] = slice_lo[idim] + 1;
 
       }
       else 
       {
          // moving realbox.lo and reabox.hi to nearest coarsenable grid point //
          int index_lo = floor( ( (slice_realbox.lo(idim) - (real_box.lo(idim)) )                                         / geom[0].CellSize(idim) ) +fac * 1E-10);
          int index_hi = ceil( ( (slice_realbox.hi(idim) - (real_box.lo(idim)) )                                         / geom[0].CellSize(idim) ) - fac * 1E-10);
          bool modify_cr = true;
         
          // The input coarsening ratio and lo and hi may require reduction in the ratio //
          while ( modify_cr == true ) 
          {
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

              //If modified index.hi is > baselinebox.hi, reduce coarsening ratio,                          and provide more points that asked for //               
              if ( hi_new * geom[0].CellSize(idim) > real_box.hi(idim) ) 
              {
                 cr_ratio[idim] = cr_ratio[idim]/2;
                 modify_cr = true;
              } 
              
              int ncells = (hi_new - lo_new);

              // If refined cells is not an integer multiple of coarsening ratio,                            then reduce coarsening ratio by factor of 2 // 
 
              if( ( ncells % cr_ratio[idim] ) != 0 )
              {
                  cr_ratio[idim] = cr_ratio[idim]/2;
                  modify_cr = true;
              }
                 
              if ( modify_cr == false ) 
              {
                 index_lo = lo_new;
                 index_hi = hi_new;
              }

              slice_lo[idim] = index_lo;   
              slice_hi[idim] = index_hi - 1; // since default is cell-centered    
          }
          slice_realbox.setLo( idim, index_lo * geom[0].CellSize(idim) );
          slice_realbox.setHi( idim, index_hi * geom[0].CellSize(idim) );
          slice_cc_nd_box.setLo( idim, slice_realbox.lo(idim) + fac );
          slice_cc_nd_box.setHi( idim, slice_realbox.hi(idim) - fac );
       }
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
       amrex::Print() << " Input parameters for the Slice dim : " << idim ;
       amrex::Print() << " slice real box lo = " << slice_realbox.lo(idim) << " hi = " << slice_realbox.hi(idim)  ;   
       amrex::Print() << " cr " << cr_ratio[idim] << "\n";
    }


}
