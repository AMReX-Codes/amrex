#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;

using std::string;

void
Nyx::compute_hydro_sources(amrex::Real time, amrex::Real dt, amrex::Real a_old, amrex::Real a_new,
                           MultiFab& S_border, MultiFab& D_border, 
                           MultiFab& ext_src_old, MultiFab& hydro_src, 
                           MultiFab& grav_vector, MultiFab& divu_cc,
                           bool init_flux_register, bool add_to_flux_register) 
{
    amrex::Print() << "Computing the hydro sources ... " << std::endl;
    const int finest_level = parent->finestLevel();
    const Real* dx = geom.CellSize();

    MultiFab fluxes[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        fluxes[j].define(getEdgeBoxArray(j), dmap, NUM_STATE, 0);
        fluxes[j].setVal(0.0);
    }

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister* fine    = 0;
    FluxRegister* current = 0;

    if (do_reflux)
    {
      if (level < finest_level)
      {
         fine = &get_flux_reg(level+1);
         if (init_flux_register)
             fine->setVal(0);

       } 
       if (level > 0) {
         current = &get_flux_reg(level);
       }
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    FArrayBox flux[BL_SPACEDIM], u_gdnv[BL_SPACEDIM];

    for (MFIter mfi(S_border,true); mfi.isValid(); ++mfi)
    {
        const Box& bx        = mfi.tilebox();

        FArrayBox& state     = S_border[mfi];
        FArrayBox& dstate    = D_border[mfi];

        // Allocate fabs for fluxes.
        for (int i = 0; i < BL_SPACEDIM ; i++) {
            const Box &bxtmp = amrex::surroundingNodes(bx, i);
            flux[i].resize(bxtmp, NUM_STATE);
            u_gdnv[i].resize(amrex::grow(bxtmp, 1), 1);
            u_gdnv[i].setVal(1.e200);
        }
        fort_make_hydro_sources
            (&time, bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(state),
             BL_TO_FORTRAN(u_gdnv[0]),
             BL_TO_FORTRAN(u_gdnv[1]),
             BL_TO_FORTRAN(u_gdnv[2]),
             BL_TO_FORTRAN(ext_src_old[mfi]),
             BL_TO_FORTRAN(hydro_src[mfi]),
             BL_TO_FORTRAN(divu_cc[mfi]),
             BL_TO_FORTRAN(grav_vector[mfi]),
             dx, &dt,
             BL_TO_FORTRAN(flux[0]),
             BL_TO_FORTRAN(flux[1]),
             BL_TO_FORTRAN(flux[2]),
             &a_old, &a_new, 
             &print_fortran_warnings);

        for (int i = 0; i < BL_SPACEDIM; ++i) 
          fluxes[i][mfi].copy(flux[i], mfi.nodaltilebox(i));
        
    } // end of MFIter loop

    } // end of parallel

    if (add_to_flux_register)
    {
       if (do_reflux) {
         if (current) {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
             current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, 1);
           }
         }
         if (fine) {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
	         fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.,FluxRegister::ADD);
           }
         }
       }
    }
}
