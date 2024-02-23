#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>

// compute the integral of f dV, where f is one of the fields in the plotfile
// this understands axisymmetric geometry.

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile;
    std::string varname_arg;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-v" || name == "--variable") {
            varname_arg = amrex::get_command_argument(++farg);
        } else {
            break;
        }
        ++farg;
    }

    if (pltfile.empty() && farg <= narg) {
        pltfile = amrex::get_command_argument(farg);
    }

    if (pltfile.empty()) {
        amrex::Print()
            << "\n"
            << " Compute the integral of f dV, where f is the field specified via -v.\n"
            << " Works with 1-, 2-, or 3-d datasets, including 2-d axisymmetric.\n"
            << "\n"
            << " Usage:\n"
            << "    fvolumesum [-v variable] plotfile\n"
            << "\n"
            << " args [-v|--variable]  varname      : the field to integrate over\n"
            << "\n"
            << std::endl;
        return;
    }


    PlotFileData pf(pltfile);
    const Vector<std::string>& var_names_pf = pf.varNames();

    std::string var_name;

    if (varname_arg.empty()) {
        var_name = "density";
    } else {
        var_name = varname_arg;
    }

    // make sure that variable name is valid

    bool found = false;
    for (auto const& vpf : var_names_pf) {
        if (var_name == vpf) {
            found = true;
            break;
        }
    }

    if (! found) {
        amrex::Abort("Error: invalid variable name");
    }


    Array<Real,AMREX_SPACEDIM> problo = pf.probLo();

    const int dim = pf.spaceDim();

    int fine_level = pf.finestLevel();

    Vector<Real> pos;

    Real lsum = 0.0;

    int coord = pf.coordSys();

    constexpr Real pi = 3.1415926535897932;

    for (int ilev = 0; ilev <= fine_level; ++ilev) {

        Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(ilev);

        Real volfac = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

        if (coord == 1) {
            AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == 2);
            // axisymmetric V = pi (r_r**2 - r_l**2) * dz
            //                = pi dr * dz * (r_r + r_l)
            //                = 2 pi r dr dz
            volfac *= 2 * pi; // 2 * pi * dr * dz part here
        } else if (coord == 2) {
            AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == 1);
            // 1-d spherical V = 4/3 pi (r_r**3 - r_l**3)
            volfac *= (4.0_rt/3.0_rt) * pi; // 4/3 * pi * dr part here
        }

        auto xlo = problo[0];
        auto dx0 = dx[0];
        AMREX_ASSERT(coord == 0 || coord == 1 || coord == 2);
        auto f_vol = [=] AMREX_GPU_DEVICE (int i) {
                         if (coord == 0) {
                             return volfac;
                         } else if (coord == 1) {
                             return volfac * (xlo + (Real(i)+0.5_rt)*dx0);
                         } else {
                             Real r_r = xlo + Real(i+1)*dx0;
                             Real r_l = xlo + Real(i  )*dx0;
                             return volfac * (r_r*r_r + r_l*r_r + r_l*r_l);
                         }
                     };

        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                pf.boxArray(ilev+1), ratio);
            const MultiFab& mf = pf.get(ilev, var_name);
            auto const& ima = mask.const_arrays();
            auto const& ma = mf.const_arrays();
            lsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mf,
                    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                        -> GpuTuple<Real>
                    {
                        return { (ima[bno](i,j,k) == 0) ? ma[bno](i,j,k)*f_vol(i) : 0._rt };
                    });
        } else {
            const MultiFab& mf = pf.get(ilev, var_name);
            auto const& ma = mf.const_arrays();
            lsum += ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, mf,
                    [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                        -> GpuTuple<Real>
                    {
                        return { ma[bno](i,j,k)*f_vol(i) };
                    });
        }
    }

    ParallelDescriptor::ReduceRealSum(lsum);

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "integral of " << var_name << " = " << lsum << std::endl;

    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
