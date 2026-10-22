#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_Print.H>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int failures = 0;
    {
        // Set parameters before the first instance of this specialization.
        amrex::ParmParse pp("particles");
        bool prepost = false;
        bool unlink = true;
        pp.query("use_prepost", prepost);
        pp.query("do_unlink", unlink);

        const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(3));
        amrex::RealBox bounds;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { bounds.setLo(d,0); bounds.setHi(d,1); }
        const int periodic[AMREX_SPACEDIM] = {AMREX_D_DECL(0, 0, 0)};
        const amrex::Geometry geom(domain, &bounds, 0, periodic);
        const amrex::BoxArray ba(domain);
        const amrex::DistributionMapping dm(ba);
        using PC = amrex::ParticleContainerRTSoA<>;
        PC first(geom, dm, ba);
        PC second(geom, dm, ba);

        auto check = [&](const char* name, const PC& pc, bool prepost, bool unlink) {
            amrex::Print() << name << ": use_prepost=" << pc.GetUsePrePost()
                           << " do_unlink=" << pc.GetUseUnlink()
                           << " expected=" << prepost << ',' << unlink << '\n';
            if (pc.GetUsePrePost() != prepost || pc.GetUseUnlink() != unlink) {
                ++failures;
            }
        };
        check("first", first, prepost, unlink);
        check("second", second, prepost, unlink);
        first.SetUsePrePost(!prepost);
        first.SetUseUnlink(!unlink);
        check("first after setters", first, !prepost, !unlink);
        check("second after first setters", second, prepost, unlink);
        PC third(geom, dm, ba);
        check("third", third, prepost, unlink);
        amrex::Print() << "Failed checks: " << failures << '\n';
    }
    amrex::Finalize();
    return failures == 0 ? 0 : 1;
}
