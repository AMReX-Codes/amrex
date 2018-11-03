
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BLFort.H>

using namespace amrex;

void flux_to_dudt_c (Box const& bx,
                     FArrayBox& dudtfab,
                     FArrayBox const& fxfab,
                     FArrayBox const& fyfab,
         FArrayBox const& fzfab,
                     Array<Real,AMREX_SPACEDIM> const& dxinv);

extern "C" {
    void flux_to_dudt_f (const int* lo, const int* hi,
                         Real* dudt, const int* ulo, const int* uhi,
                         const Real* fx, const int* xlo, const int* xhi,
                         const Real* fy, const int* ylo, const int* yhi,
                         const Real* fz, const int* zlo, const int* zhi,
                         const Real* dxinv, const int* ncomp);
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        const Box bx(IntVect(0), IntVect(63));
        const int ncomp = 5;
        FArrayBox dudtfab(bx,ncomp);
        FArrayBox fxfab(amrex::convert(bx,IntVect::TheDimensionVector(0)),ncomp);
        FArrayBox fyfab(amrex::convert(bx,IntVect::TheDimensionVector(1)),ncomp);
        FArrayBox fzfab(amrex::convert(bx,IntVect::TheDimensionVector(2)),ncomp);

        fxfab.setVal(1.0);
        fyfab.setVal(1.0);
        fzfab.setVal(1.0);

        Array<Real,3> dxinv {10.,10.,10.};

        double t0 = amrex::second();

        for (int i = 0; i < 1000; ++i) {
            asm volatile("" : "+r" (dxinv[0]));
            flux_to_dudt_c(bx, dudtfab, fxfab, fyfab, fzfab, dxinv);
        }

        double t1 = amrex::second();

        for (int i = 0; i < 1000; ++i) {
            asm volatile("" : "+r" (dxinv[0]));
            flux_to_dudt_f(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_ANYD(dudtfab),
                           BL_TO_FORTRAN_ANYD(fxfab),
                           BL_TO_FORTRAN_ANYD(fyfab),
                           BL_TO_FORTRAN_ANYD(fzfab),
                           dxinv.data(),
                           &ncomp);
        }

        double t2 = amrex::second();

        amrex::Print() << " C++ time: " << t1-t0 << ", Fortran time: " << t2-t1 << std::endl;
    }
    amrex::Finalize();
}

