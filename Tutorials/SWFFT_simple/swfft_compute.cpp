
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BaseFab.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;

void
swfft_compute(MultiFab& phi_spatial, MultiFab& phi_dft, Geometry& geom, int verbose)
{
    const BoxArray& ba = phi_dft.boxArray();
    amrex::Print() << "BA " << ba << std::endl;
    const DistributionMapping& dm = phi_dft.DistributionMap();

    if (ba.size() != ParallelDescriptor::NProcs()) {
      amrex::Error("Need same number of MPI processes as grids");
      exit(0);
    }

    if (phi_spatial.nGrow() != 0 || phi_dft.nGrow() != 0) 
       amrex::Error("Current implementation requires that both phi_spatial and phi_dft have no ghost cells");

    // We assume that all grids have the same size hence 
    // we have the same nx,ny,nz on all ranks
    int nx = ba[0].size()[0];
    int ny = ba[0].size()[1];
#if (AMREX_SPACEDIM == 2)
    int nz = 1;
#elif (AMREX_SPACEDIM == 3)
    int nz = ba[0].size()[2];
#endif

    Box domain(geom.Domain());

    int nbx = domain.length(0) / nx;
    int nby = domain.length(1) / ny;
#if (AMREX_SPACEDIM == 2)
    int nbz = 1;
#elif (AMREX_SPACEDIM == 3)
    int nbz = domain.length(2) / nz;
#endif
    int nboxes = nbx * nby * nbz;
    if (nboxes != ba.size()) 
       amrex::Error("NBOXES NOT COMPUTED CORRECTLY");
    amrex::Print() << "Number of boxes:\t" << nboxes << std::endl;

    Vector<int> rank_mapping;
    rank_mapping.resize(nboxes);

    DistributionMapping dmap = phi_spatial.DistributionMap();

    for (int ib = 0; ib < nboxes; ++ib)
    {
        int i = ba[ib].smallEnd(0) / nx;
        int j = ba[ib].smallEnd(1) / ny;
#if (AMREX_SPACEDIM == 2)
	int k = 0;
#elif (AMREX_SPACEDIM == 3)
	int k = ba[ib].smallEnd(2) / nz;
#endif

        // This would be the "correct" local index if the data wasn't being transformed
        int local_index = k*nbx*nby + j*nbx + i;

        // This is what we pass to dfft to compensate for the Fortran ordering
        //      of amrex data in MultiFabs.
        // int local_index = i*nby*nbz + j*nbz + k;

        rank_mapping[local_index] = dmap[ib];

        if (verbose)
          amrex::Print() << "LOADING RANK NUMBER " << dmap[ib] << " FOR GRID NUMBER " << ib 
                         << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }

    Real start_time = amrex::second();

    // Assume for now that nx = ny = nz
#if (AMREX_SPACEDIM == 2)
    int Ndims[3] = { 1, nby, nbx};
    int     n[3] = { 1, domain.length(1), domain.length(0)};
#elif (AMREX_SPACEDIM == 3)
    int Ndims[3] = { nbz, nby, nbx };
    int     n[3] = { domain.length(2), domain.length(1), domain.length(0)};
#endif
    hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
    hacc::Dfft dfft(d);
    
    for (MFIter mfi(phi_spatial,false); mfi.isValid(); ++mfi)
    {
       int gid = mfi.index();

       size_t local_size  = dfft.local_size();
   
       std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
       std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

       a.resize(nx*ny*nz);
       b.resize(nx*ny*nz);

       dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

       // *******************************************
       // Copy real data from Rhs into real part of a -- no ghost cells and
       // put into C++ ordering (not Fortran)
       // *******************************************
       complex_t zero(0.0, 0.0);
       size_t local_indx = 0;
       for(size_t k=0; k<(size_t)nz; k++) {
        for(size_t j=0; j<(size_t)ny; j++) {
         for(size_t i=0; i<(size_t)nx; i++) {

           complex_t temp(phi_spatial[mfi].dataPtr()[local_indx],0.);
           a[local_indx] = temp;
      	   local_indx++;

         }
       }
      }

//  *******************************************
//  Compute the forward transform
//  *******************************************
       dfft.forward(&a[0]);
    
//  *******************************************
//  Redistribute from z-pencils to blocks
//  *******************************************
       d.redistribute_2_to_3(&a[0],&b[0],2);
       
       size_t global_size  = dfft.global_size();
       double fac;

       // fac = sqrt(1.0 / (double)global_size);       
       fac = 1.0; // Overwrite fac

       local_indx = 0;
       // int local_indx_p = 0;
       for(size_t k=0; k<(size_t)nz; k++) {
        for(size_t j=0; j<(size_t)ny; j++) {
         for(size_t i=0; i<(size_t)nx; i++) {

           // Divide by 2 pi N
           phi_dft[mfi].dataPtr()[local_indx] = fac * std::abs(b[local_indx]);
       	   local_indx++;

         }
        }
       }
    }

}
