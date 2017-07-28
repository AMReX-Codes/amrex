

#include <fstream>

#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BaseUmap.H>
#include <AMReX_BaseUmap_f.H>
#include <AMReX_ArrayLim.H>

//BL_FORT_PROC_DECL(FILLFAB,fillfab)(Real* d, const int* nx, const int* ny);

using namespace amrex;

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    const int NY = 4;
    const int NX = 3;
    const int NZ = 2;

    Box bx(IntVect::TheZeroVector(),IntVect(D_DECL(NX-1,NY-1,NZ-1)));

    FArrayBox fab(bx,2);
    fab.setVal(1.e200);

    BaseUmap<Real> umap(bx, 2);


    Real v;
    v = 2.0;

#pragma omp parallel
    {
        std::vector<double> localdata;
        std::vector<int> localkeys;
        localdata.resize(0);
        localkeys.resize(0);

#pragma omp for
        for(int k=0;k<NZ;k++)
        {
            for(int j=0;j<NY;j++)
            {
                for(int i=0;i<NX;i++)
                {
                    int ind=i+NX*j+NX*NY*k;

                    if(ind%2 == 0)
                    {
                        umap.setVal_Threadsafe( v, IntVect(i,j,k), 0, 1,localdata,localkeys);
                        //umap.setVal( v, IntVect(i,j,k), 0, 1);
                    }
                    if(ind%2 == 1)
                    {
                        umap.setVal_Threadsafe( v+1, IntVect(i,j,k), 0, 1,localdata,localkeys);
                        //umap.setVal( v+1, IntVect(i,j,k), 0, 1);
                    }
                }
            }
        }

#pragma omp critical
        {
            umap.stitch_data_critical(localdata,localkeys);
        }
    }

    umap.set_key_indices();

    for (IntVect iv=bx.smallEnd(); iv<=bx.bigEnd(); bx.next(iv)) {
        v = umap.getVal(iv, 0, 1);
        std::cout << iv << " has value " << v <<std::endl;
        //
    }

    // Convert to references.
    // umap.nPts() returns a long, fort_umap_norm expects an int
    // Keytable location may not be quite right
    int npts = umap.nPts();
    int max_mv = umap.MaxMV();
    int ncomp = umap.nComp();
    Real norm =  fort_umap_norm(ARLIM_3D(umap.box().loVect()), ARLIM_3D(umap.box().hiVect()),
            umap.dataPtr(),&npts, 
            umap.keyTablePtr(), ARLIM_3D(umap.box().loVect()), ARLIM_3D(umap.box().hiVect()),
            &max_mv, &ncomp, 0);


    std::cout << "Norm (from fotran): " << norm << std::endl;


    amrex::Finalize();

    return 0;
}
