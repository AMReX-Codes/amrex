#include <AMReX.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Gpu.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        constexpr int ncomp = 36;

        Box box(IntVect(0), IntVect(31));
        IArrayBox fab(box, ncomp);

        Array<Array4<int>,ncomp> aa;
        for (int n = 0; n < ncomp; ++n) {
            aa[n] = fab.array(n);
        }

        // Each Array4 has a size of 64. aa with 36 Array4s has a size of
        // 2304 bytes. So the following kernle has a size greater than 2048
        // bytes, the parameter size limit of oneAPI. In AMReX, we manually
        // copy the kernel to the device to work around the limit.
        ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            for (int n = 0; n < ncomp; ++n) {
                aa[n](i,j,k) = n;
            }
        });
        Gpu::streamSynchronize();

        for (int n = 0; n < ncomp; ++n) {
            auto s = fab.template sum<RunOn::Device>(n);
            AMREX_ALWAYS_ASSERT(Long(s) == box.numPts()*n);
        }
    }

    {
        constexpr int ncomp = 18;

        Box ccbox(IntVect(0), IntVect(31));
        Box ndbox = amrex::convert(ccbox, IntVect(1));
        IArrayBox ccfab(ccbox, ncomp);
        IArrayBox ndfab(ndbox, ncomp);

        Array<Array4<int>,ncomp> ccaa;
        Array<Array4<int>,ncomp> ndaa;
        for (int n = 0; n < ncomp; ++n) {
            ccaa[n] = ccfab.array(n);
            ndaa[n] = ndfab.array(n);
        }

        ParallelFor(ccbox, ndbox,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        for (int n = 0; n < ncomp; ++n) {
                            ccaa[n](i,j,k) = n;
                        }
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        for (int n = 0; n < ncomp; ++n) {
                            ndaa[n](i,j,k) = n;
                        }
                    });
        Gpu::streamSynchronize();

        for (int n = 0; n < ncomp; ++n) {
            auto ccs = ccfab.template sum<RunOn::Device>(n);
            auto nds = ndfab.template sum<RunOn::Device>(n);
            AMREX_ALWAYS_ASSERT(Long(ccs) == ccbox.numPts()*n &&
                                Long(nds) == ndbox.numPts()*n);
        }
    }

    {
        constexpr int ncomp = 12;

        Box ccbox(IntVect(0), IntVect(31));
        Box fcbox = amrex::convert(ccbox, IntVect::TheDimensionVector(AMREX_SPACEDIM-1));
        Box ndbox = amrex::convert(ccbox, IntVect(1));
        IArrayBox ccfab(ccbox, ncomp);
        IArrayBox fcfab(fcbox, ncomp);
        IArrayBox ndfab(ndbox, ncomp);

        Array<Array4<int>,ncomp> ccaa;
        Array<Array4<int>,ncomp> fcaa;
        Array<Array4<int>,ncomp> ndaa;
        for (int n = 0; n < ncomp; ++n) {
            ccaa[n] = ccfab.array(n);
            fcaa[n] = fcfab.array(n);
            ndaa[n] = ndfab.array(n);
        }

        ParallelFor(ccbox, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        ccaa[n](i,j,k) = n;
                    },
                    fcbox, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        fcaa[n](i,j,k) = n;
                    },
                    ndbox, ncomp,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                    {
                        ndaa[n](i,j,k) = n;
                    });
        Gpu::streamSynchronize();

        for (int n = 0; n < ncomp; ++n) {
            auto ccs = ccfab.template sum<RunOn::Device>(n);
            auto fcs = fcfab.template sum<RunOn::Device>(n);
            auto nds = ndfab.template sum<RunOn::Device>(n);
            AMREX_ALWAYS_ASSERT(Long(ccs) == ccbox.numPts()*n &&
                                Long(fcs) == fcbox.numPts()*n &&
                                Long(nds) == ndbox.numPts()*n);
        }
    }

    amrex::Finalize();
}
