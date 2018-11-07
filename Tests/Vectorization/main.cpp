
#include <AMReX.H>
#include <AMReX_Print.H>

#include <kdecl.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        const Box bx(IntVect(0), IntVect(63));
        const Box& bxg2 = amrex::grow(bx,2);
        const int ncomp = 7;
        FArrayBox ufab(bxg2,ncomp);
        FArrayBox qfab(bxg2,8);
        FArrayBox dudtfab(bx,ncomp);
        FArrayBox fxfab(amrex::convert(bx,IntVect::TheDimensionVector(0)),ncomp);
        FArrayBox fyfab(amrex::convert(bx,IntVect::TheDimensionVector(1)),ncomp);
        FArrayBox fzfab(amrex::convert(bx,IntVect::TheDimensionVector(2)),ncomp);

        ufab.setVal(3.0);
        fxfab.setVal(1.0);
        fyfab.setVal(1.0);
        fzfab.setVal(1.0);

        Array<Real,3> dxinv {10.,10.,10.};

        // ctoprim
        {
            ctoprim_f(BL_TO_FORTRAN_BOX(bxg2),
                      BL_TO_FORTRAN_ANYD(ufab),
                      BL_TO_FORTRAN_ANYD(qfab));

            double t0 = amrex::second();

            for (int i = 0; i < 1000; ++i) {
                __asm__ __volatile__("" : : : "memory");
                ctoprim_f(BL_TO_FORTRAN_BOX(bxg2),
                          BL_TO_FORTRAN_ANYD(ufab),
                          BL_TO_FORTRAN_ANYD(qfab));
            }
            
            double t1 = amrex::second();
            
            for (int i = 0; i < 1000; ++i) {
                __asm__ __volatile__("" : : : "memory");
                ctoprim_c_simd(bxg2, ufab, qfab);
            }
            
            double t2 = amrex::second();
            
            for (int i = 0; i < 1000; ++i) {
                __asm__ __volatile__("" : : : "memory");
                ctoprim_c_nosimd(bxg2, ufab, qfab);
            }
            
            double t3 = amrex::second();
            
            amrex::Print() << "ctoprim: Fortran time: " << t1-t0 << "\n"
                           << "         C++ w/ simd time: " << t2-t1  << "\n"
                           << "         C++ w/o simd time: " << t3-t2 << std::endl;
        }

        // flux_to_dudt 
        {
            flux_to_dudt_f(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_ANYD(dudtfab),
                           BL_TO_FORTRAN_ANYD(fxfab),
                           BL_TO_FORTRAN_ANYD(fyfab),
                           BL_TO_FORTRAN_ANYD(fzfab),
                           dxinv.data(),
                           &ncomp);

            double t0 = amrex::second();

            for (int i = 0; i < 1000; ++i) {
                __asm__ __volatile__("" : : : "memory");
                flux_to_dudt_f(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_ANYD(dudtfab),
                               BL_TO_FORTRAN_ANYD(fxfab),
                               BL_TO_FORTRAN_ANYD(fyfab),
                               BL_TO_FORTRAN_ANYD(fzfab),
                               dxinv.data(),
                               &ncomp);
            }
            
            double t1 = amrex::second();
            
            for (int i = 0; i < 1000; ++i) {
                __asm__ __volatile__("" : : : "memory");
                flux_to_dudt_c_simd(bx, dudtfab, fxfab, fyfab, fzfab, dxinv);
            }
            
            double t2 = amrex::second();
            
            for (int i = 0; i < 1000; ++i) {
                __asm__ __volatile__("" : : : "memory");
                flux_to_dudt_c_nosimd(bx, dudtfab, fxfab, fyfab, fzfab, dxinv);
            }
            
            double t3 = amrex::second();
            
            amrex::Print() << "flux_to_dudt: Fortran time: " << t1-t0 << "\n"
                           << "              C++ w/ simd time: " << t2-t1  << "\n"
                           << "              C++ w/o simd time: " << t3-t2 << std::endl;
        }
    }
    amrex::Finalize();
}

