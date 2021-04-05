
#include <AMReX.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BLProfiler.H>

using namespace amrex;

void fused_test (iMultiFab const& mfa, iMultiFab const& mfb, iMultiFab& mfc)
{
    Gpu::FuseSafeGuard fsg(true);
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
        Array4<int const> const& a = mfa.const_array(mfi);
        Array4<int const> const& b = mfb.const_array(mfi);
        Array4<int      > const& c = mfc.array(mfi);
        Box const& vbx = mfi.validbox();
        amrex::Gpu::Register(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            c(i,j,k,0) = a(i,j,k) + b(i,j,k);
        });
        amrex::Gpu::Register(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            c(i,j,k,1) = a(i,j,k) - b(i,j,k);
        });
    }
    amrex::Gpu::LaunchFusedKernels();
}

void unfused_test (iMultiFab const& mfa, iMultiFab const& mfb, iMultiFab& mfc)
{
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
        Array4<int const> const& a = mfa.const_array(mfi);
        Array4<int const> const& b = mfb.const_array(mfi);
        Array4<int      > const& c = mfc.array(mfi);
        Box const& vbx = mfi.validbox();
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            c(i,j,k,0) = a(i,j,k) + b(i,j,k);
        });
        amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            c(i,j,k,1) = a(i,j,k) - b(i,j,k);
        });
    }
}

int verify (iMultiFab& mfc)
{
    int min0 = mfc.min(0);
    int min1 = mfc.min(1);
    int max0 = mfc.max(0);
    int max1 = mfc.max(1);
    bool r;
    if (min0 != 4 || max0 != min0 || min1 != -2 || max1 != min1) {
        amrex::Print() << "Failed!" << std::endl;
        r = false;
    } else {
        amrex::Print() << "Success" << std::endl;
        r = true;
    }
    mfc.setVal(-1);
    return r;
}

int main(int argc, char* argv[])
{
    bool r = true;
    Real tf, tuf;
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE("main()");

        Box bx(IntVect(0),IntVect(127));
        BoxArray ba(bx);
        ba.maxSize(16);
        DistributionMapping dm{ba};

        iMultiFab mfa(ba,dm,1,0);
        iMultiFab mfb(ba,dm,1,0);
        iMultiFab mfc(ba,dm,2,0);
        mfa.setVal(1);
        mfb.setVal(3);
        mfc.setVal(-1);

        {
            BL_PROFILE("fused-test1");
            fused_test(mfa,mfb,mfc);
        }
        r = r && verify(mfc);

        {
            BL_PROFILE("fused-test2");
            Real t = amrex::second();
            fused_test(mfa,mfb,mfc);
            tf = amrex::second() - t;
        }
        r = r && verify(mfc);

        {
            BL_PROFILE("unfused-test");
            Real t = amrex::second();
            unfused_test(mfa,mfb,mfc);
            tuf = amrex::second() - t;
        }
        r = r && verify(mfc);
    }
    amrex::Finalize();

    if (r) {
        std::cout << "\nTest passes.\n";
    } else {
        std::cout << "\nTest fails.\n";
    }
    std::cout << "Unfused kernels take " << tuf << " seconds.\n"
              << "Fused   kernels take " << tf  << " seconds." << std::endl;

    return r ? EXIT_SUCCESS : EXIT_FAILURE;
}
