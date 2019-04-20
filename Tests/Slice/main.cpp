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

void CreateSlice(Vector<MultiFab*> mf, Vector<Geometry> geom, RealBox slice_realbox, Vector<int> slice_ncells, int slice_grid_size);
void InitializeVariables(Vector<MultiFab*> mf, Vector<Geometry> geom);

int main(int argc, char* argv[])
{    
    amrex::Initialize(argc,argv);
    { 
    const int nghost = 0;
    int ncells, max_grid_size, ncomp, nlevs, nppc;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("ncomp", ncomp);
    pp.get("nlevs", nlevs);
    pp.get("nppc", nppc);
    
    
    AMREX_ALWAYS_ASSERT(nlevs < 2); // relax this later

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

    InitializeVariables( GetVecOfPtrs(mf), geom );
    InitializeVariables( GetVecOfPtrs(rho_mf), geom );
    InitializeVariables( GetVecOfPtrs(Ex_mf), geom );
    InitializeVariables( GetVecOfPtrs(Ey_mf), geom );
    InitializeVariables( GetVecOfPtrs(Ez_mf), geom );
    InitializeVariables( GetVecOfPtrs(Bx_mf), geom );
    InitializeVariables( GetVecOfPtrs(By_mf), geom );
    InitializeVariables( GetVecOfPtrs(Bz_mf), geom );

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
    
    int slice_grid_size;

    ParmParse ppg("slice");
    ppg.queryarr("dom_lo",slo,0,AMREX_SPACEDIM);
    ppg.queryarr("dom_hi",shi,0,AMREX_SPACEDIM);
    ppg.queryarr("n_cells",slice_cells,0,AMREX_SPACEDIM);
    ppg.query("max_grid_size",slice_grid_size);
    amrex::RealBox slice_realbox;
    slice_realbox.setLo(slo); 
    slice_realbox.setHi(shi); 

    CreateSlice( GetVecOfPtrs(mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(rho_mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(Ex_mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(Ey_mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(Ez_mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(Bx_mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(By_mf),geom,slice_realbox,slice_cells,slice_grid_size);
    CreateSlice( GetVecOfPtrs(Bz_mf),geom,slice_realbox,slice_cells,slice_grid_size);
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

void CreateSlice(Vector<MultiFab*> mf, Vector<Geometry> geom, RealBox slice_realbox, Vector<int> slice_ncells, int slice_grid_size)
{

    // assumption 1 : nghost = 0; 
    const int nghost = 0;
    int is_per[AMREX_SPACEDIM];

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is_per[i] = 1;
  
    int nlevs = geom.size();
    int ncomp = (*mf[0]).nComp();
    bool genslice = true; 
    int max_ratio = 1;
    bool coarsen = false; 
    IntVect cr_ratio(AMREX_D_DECL(1,1,1));
    const auto conversionType = (*mf[0]).ixType();

    const RealBox& real_box = geom[0].ProbDomain(); 

    // ensuring that index space for slice is same as domain // 
    IntVect slice_lo(AMREX_D_DECL(0,0,0));
    IntVect slice_hi(AMREX_D_DECL(1,1,1)); 

    // Define slice box and error checks //
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
       
       double fac = ( 1.0-conversionType[idim] )*geom[0].CellSize(idim) * 0.5;
//       double compare = real_box.lo(idim) + geom[0].CellSize(idim) *0.5;
       
       slice_lo[idim] = round((slice_realbox.lo(idim) - real_box.lo(idim))/geom[0].CellSize(idim));
       slice_hi[idim] = round((slice_realbox.hi(idim) - real_box.lo(idim))/geom[0].CellSize(idim)); 

       if ( slice_lo[idim] > fac && fac > 0 ) {
            slice_lo[idim] =  round((slice_realbox.lo(idim) - (real_box.lo(idim) + fac))/geom[0].CellSize(idim));
            slice_hi[idim] =  round((slice_realbox.hi(idim) - (real_box.lo(idim) + fac))/geom[0].CellSize(idim));
       }
        
       if ( ( slice_hi[idim] - slice_lo[idim]) == 0) {
          // Only 1 cell is required if hi and lo are equal //
          slice_hi[idim] = slice_lo[idim] + 2;
          slice_ncells[idim] = 1;
       }
       else {

          int refinedcells =  ((slice_realbox.hi(idim) - slice_realbox.lo(idim)) / 
                                                         geom[0].CellSize(idim));

          // compare cell sizes of slice and computational domain //
          if (refinedcells >= slice_ncells[idim]) 
          {
             if (refinedcells % slice_ncells[idim] != 0) 
             {
                amrex::Abort( " SLICEERROR :: cell size of slice is not an integer                                           multiple of the computaitonal domain's cell size" );
             } 
             if ( refinedcells > slice_ncells[idim] && slice_ncells[idim]>0) 
             {
                coarsen = true; 
                cr_ratio[idim] =  refinedcells / slice_ncells[idim] ;
                if ( slice_grid_size % cr_ratio[idim] != 0) 
                {
                   amrex::Abort("SLICE_ERROR :: Max grid size of slice is not an                                              integer multiple of coarsening ratio ");
                }
                if (max_ratio < cr_ratio[idim] ) 
                {
                   max_ratio = cr_ratio[idim];
                }
             }
          }
          else 
          {
             genslice = false; // no slice generation in this case
             amrex::Abort(" SLICEERROR : The cell size for the required diagnostic                                       slice is more refined than the simulation domain cell size.                                  Please change input such that cell size >= domain cell size");
          }
       }
       --slice_hi[idim]; // since default index type is cc  // -= 1; 
    }


    // Slice generation only if slice cell size >= domain cell size //
    if (genslice == true)
    {

    // Default index type for slice is cell-centered //
    Box slice(slice_lo, slice_hi);
    


    // Convert from cc to index type of parent multifab //
    IntVect convertTypeofSlice(AMREX_D_DECL(0,0,0));
    bool slicetypeToBeConverted = 0;
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim )
    {
       if ( conversionType.nodeCentered(idim) ) {
          convertTypeofSlice[idim] = 1;
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
       smf[0].reset(new MultiFab(amrex::convert(sba[0],convertTypeofSlice),sdmap[0],ncomp,nghost));
    }
    else {
       smf[0].reset(new MultiFab(sba[0],sdmap[0],ncomp,nghost));
    }
    smf[0]->setVal(1);
    VisMF::Write((*smf[0]),"vismf_init");
    // Copy data from domain to slice that has same cell size as that of the domain mf.  
    (*smf[0]).ParallelCopy((*mf[0]), 0, 0, ncomp);  

    if ( coarsen == false ) {
       amrex::Print() << " Cell sizes are equal. No averaging required. " << "\n";
       VisMF::Write((*smf[0]),"vismf_output");
    }
    else if ( coarsen == true ) {
       
       Vector<BoxArray> crse_ba(nlevs);

       //Ensures that box arrays are the same and then coarsen using IntVect cr_ratio
       crse_ba[0] = sba[0];
       crse_ba[0].coarsen(cr_ratio);

       // The input values of max_grid_size is factored by ratio in the coarsened slice //
       int cs_grid_size = double(sba[0].size())*max_ratio ; 
       crse_ba[0].maxSize(cs_grid_size);
       AMREX_ALWAYS_ASSERT(crse_ba[0].size() == sba[0].size());

       // constructing coarsened slice as per user-input if s_cells<ncells // 
       Vector<std::unique_ptr<MultiFab> > cs_mf(nlevs); 

       if(slicetypeToBeConverted==1) {
          cs_mf[0].reset(new MultiFab(amrex::convert(crse_ba[0],convertTypeofSlice),sdmap[0],ncomp,nghost));
       }
       else {
          cs_mf[0].reset(new MultiFab(crse_ba[0],sdmap[0],ncomp,nghost));
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

               if( convertTypeofSlice==cctype ) {
                  amrex::amrex_avgdown(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio);
               }
               if( convertTypeofSlice == ndtype ) {
                  amrex::amrex_avgdown_nodes(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio);
               }
               if( convertTypeofSlice == xetype ) {
                  amrex::amrex_avgdown_edges(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio,0);
               }
               if( convertTypeofSlice == yetype) {
                  amrex::amrex_avgdown_edges(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio,1);
               }
               if( convertTypeofSlice == zetype ) {
                  amrex::amrex_avgdown_edges(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio,2);
               }
               if (AMREX_SPACEDIM==3) {
                  if( convertTypeofSlice == xftype) {
                     amrex::amrex_avgdown_faces(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio,0);
                  }
                  if( convertTypeofSlice == yftype ) {
                     amrex::amrex_avgdown_faces(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio,1);
                  }
                  if( convertTypeofSlice == zftype ) {
                     amrex::amrex_avgdown_faces(Dst_bx,Dst_fabox,Src_fabox,dcomp,scomp,ncomp,cr_ratio,2);
                  }
               }

               if ( mfi_dst.isValid() ) {
                  ++mfi_dst;
               }
           }
       }
       VisMF::Write((*cs_mf[0]),"vismf_output");
    }

    }
}

// assuming 6 components for all multifabs just for testing purposes //
void InitializeVariables(Vector<MultiFab*> mf, Vector<Geometry> geom)
{
    int nlevs = geom.size();
    const Real *dom_dx = geom[0].CellSize();
    for (int lev = 0; lev < nlevs; ++lev) {
        for (MFIter mfi(*mf[lev]); mfi.isValid(); ++mfi) {
            
            MultiFab &mfSrc= *mf[lev];
            auto const &mf_arr = mfSrc.array(mfi);
            const Box& bx = mfi.tilebox();
            const auto IndType = (*mf[lev]).ixType();
            const auto lo = amrex::lbound(bx); 
            const auto hi = amrex::ubound(bx); 
            for (int k = lo.z; k<=hi.z; ++k) {
               for (int j = lo.y; j<=hi.y; ++j) {
                  for (int i = lo.x; i<=hi.x; ++i) {
                      int icomp= 0;
		      mf_arr(i,j,k,icomp) = i * dom_dx[0] + (1.0-IndType[0])*dom_dx[0]*0.5; 
                      ++icomp;      
                      mf_arr(i,j,k,icomp) = j * dom_dx[1] + (1.0-IndType[1])*dom_dx[1]*0.5; 
                      ++icomp;      
                      mf_arr(i,j,k,icomp) = k * dom_dx[2] + (1.0-IndType[2])*dom_dx[2]*0.5; 
                      ++icomp;      
                      mf_arr(i,j,k,icomp) = i * dom_dx[0] + (1.0-IndType[0])*dom_dx[0]*0.5;
                      ++icomp;      
                      mf_arr(i,j,k,icomp) = j * dom_dx[1] + (1.0-IndType[1])*dom_dx[1]*0.5; 
                      ++icomp;      
                      mf_arr(i,j,k,icomp) = k * dom_dx[2] + (1.0-IndType[2])*dom_dx[2]*0.5;
                  }
               }
            }
        }
    }

}

