
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";

    {
        amrex::Box domain_box(amrex::IntVect(0), amrex::IntVect(127));
        amrex::BoxArray ba(domain_box);
        ba.maxSize(64);
        amrex::DistributionMapping dm{ba};
        amrex::MultiFab mf(ba, dm, 1, 0);

        mf.setVal(1.0);

        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            const auto& a = mf.array(mfi);
            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                a(i,j,k) = 1.*i + 10.*j + 100.*k;
            });
        }

        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.validbox();
            amrex::FArrayBox* fab = mf.fabPtr(mfi);
            amrex::launch(box,
            [=] AMREX_GPU_DEVICE (amrex::Box const& tbx)
            {
                fab->setVal(2.0);
                const auto a = fab->array();
                const auto lo = amrex::lbound(tbx);
                const auto hi = amrex::ubound(tbx);
                for         (int k = lo.z; k <= hi.z; ++k) {
                    for     (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            a(i,j,k) *= 3.5;
                        }
                    }
                }
            });
        }
    }

    amrex::Finalize();
}

