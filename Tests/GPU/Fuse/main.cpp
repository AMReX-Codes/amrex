
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BLProfiler.H>

using namespace amrex;

void fused_test (MultiFab const& mfa, MultiFab const& mfb, MultiFab& mfc)
{
    Gpu::FuseSafeGuard fsg(true);
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
        Array4<Real const> const& a = mfa.const_array(mfi);
        Array4<Real const> const& b = mfb.const_array(mfi);
        Array4<Real      > const& c = mfc.array(mfi);
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

void unfused_test (MultiFab const& mfa, MultiFab const& mfb, MultiFab& mfc)
{
    for (MFIter mfi(mfc); mfi.isValid(); ++mfi) {
        Array4<Real const> const& a = mfa.const_array(mfi);
        Array4<Real const> const& b = mfb.const_array(mfi);
        Array4<Real      > const& c = mfc.array(mfi);
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

void verify (MultiFab& mfc)
{
    Real min0 = mfc.min(0);
    Real min1 = mfc.min(1);
    Real max0 = mfc.max(0);
    Real max1 = mfc.max(1);
    if (min0 != 4.0 || max0 != min0 || min1 != -2.0 || max1 != min1) {
        amrex::Print() << "Failed!" << std::endl;
    } else {
        amrex::Print() << "Success" << std::endl;
    }
    mfc.setVal(-1.);
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE("main()");

        Box bx(IntVect(0),IntVect(127));
        BoxArray ba(bx);
        ba.maxSize(16);
        DistributionMapping dm{ba};

        MultiFab mfa(ba,dm,1,0);
        MultiFab mfb(ba,dm,1,0);
        MultiFab mfc(ba,dm,2,0);
        mfa.setVal(1.0);
        mfb.setVal(3.0);
        mfc.setVal(-1.0);

        {
            BL_PROFILE("fuesd-test1");
            fused_test(mfa,mfb,mfc);
        }
        verify(mfc);
        
        {
            BL_PROFILE("fused-test2");
            fused_test(mfa,mfb,mfc);
        }
        verify(mfc);

        {
            BL_PROFILE("unfused-test");
            unfused_test(mfa,mfb,mfc);
        }
        verify(mfc);
    }
    amrex::Finalize();
}

