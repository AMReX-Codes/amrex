#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <random>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    ULong seed;
    {
        std::random_device rd; // non-deterministic random numbers
        std::uniform_int_distribution<ULong> dist(0,std::numeric_limits<ULong>::max());
        seed = dist(rd);
        amrex::ResetRandomSeed(seed);
    }
    for (int icell = 0; icell < 10000; ++icell)
    {
        int ncells = int(amrex::Random_int(102400)) + 4;
        Box domain(IntVect(0),IntVect(ncells-1));

        for (int ieps = 0; ieps < 1000; ++ieps)
        {
            std::array<Real,AMREX_SPACEDIM> rblo{AMREX_D_DECL(Real(0.),Real(-1.),Real(-0.3))};
            std::array<Real,AMREX_SPACEDIM> rbhi{AMREX_D_DECL(Real(1.),Real( 0.),Real( 0.5))};
            if (ieps % 100 != 0) {
                auto eps = (amrex::Random() - Real(0.5)) * Real(1.e-4);
                AMREX_D_TERM(rblo[0] += eps;,
                             rblo[1] -= eps;,
                             rblo[2] += eps);
                AMREX_D_TERM(rbhi[0] -= eps;,
                             rbhi[1] += eps;,
                             rbhi[2] -= eps);
            }

            RealBox rb(rblo, rbhi);
            Geometry geom(domain, rb, 0, {AMREX_D_DECL(0,0,0)});

            auto rlo = geom.ProbLoArrayInParticleReal();
            auto rhi = geom.ProbHiArrayInParticleReal();
            auto plo = geom.ProbLoArray();
            auto dxinv = geom.InvCellSizeArray();
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto index = [&] (ParticleReal x) -> int
                {
                    return int(std::floor((x - plo[idim])*dxinv[idim]));
                };
                auto epsilon = std::numeric_limits<ParticleReal>::epsilon()
                    * std::max(ParticleReal(geom.CellSize(idim)),std::abs(rlo[idim]))
                    * ParticleReal(2.5);
                auto rlom = rlo[idim] - epsilon;
                epsilon = std::numeric_limits<ParticleReal>::epsilon()
                    * std::max(ParticleReal(geom.CellSize(idim)),std::abs(rhi[idim]))
                    * ParticleReal(2.5);
                auto rhip = rhi[idim] + epsilon;
                bool pass = (index(rlom)      == -1)
                    &&      (index(rlo[idim]) == 0 )
                    &&      (index(rhi[idim]) == ncells-1)
                    &&      (index(rhip)      == ncells);
                if (!pass) {
                    amrex::AllPrint().SetPrecision(17)
                        << "Random seed = " << seed << "\n"
                        << "RealBox: " << rb << "\n"
                        << "Geometry: " << geom << "\n"
                        << " rlo[" << idim << "] = " << rlo[idim]
                        << " rhi[" << idim << "] = " << rhi[idim]
                        << " rlo_minus = " << rlom
                        << " rhi_plus  = " << rhip << "\n"
                        << " ilo = " << index(rlo[idim])
                        << " ihi = " << index(rhi[idim])
                        << " ilo-1 = " << index(rlom)
                        << " ihi+1 = " << index(rhip)
                        << "\n";
                    amrex::Abort("Failed");
                }
            }
        }
    }
    amrex::Finalize();
}
