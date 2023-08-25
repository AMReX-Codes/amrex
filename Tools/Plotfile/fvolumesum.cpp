#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
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

        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                pf.boxArray(ilev+1), ratio);
            const MultiFab& mf = pf.get(ilev, var_name);
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.ok()) {
                    const auto& m = mask.array(mfi);
                    const auto& fab = mf.array(mfi);
                    const auto lo = amrex::lbound(bx);
                    const auto hi = amrex::ubound(bx);
                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {
                                if (m(i,j,k) == 0) { // not covered by fine
                                    Array<Real,AMREX_SPACEDIM> p
                                        = {AMREX_D_DECL(problo[0]+static_cast<Real>(i+0.5)*dx[0],
                                                        problo[1]+static_cast<Real>(j+0.5)*dx[1],
                                                        problo[2]+static_cast<Real>(k+0.5)*dx[2])};

                                    // compute the volume
                                    Real vol = std::numeric_limits<Real>::quiet_NaN();
                                    if (coord == 0) {
                                        // Cartesian
                                        vol = 1.0_rt;
                                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                            vol *= dx[idim];
                                        }
                                    } else if (coord == 1) {
                                        // axisymmetric V = pi (r_r**2 - r_l**2) * dz
                                        //                = pi dr * dz * (r_r + r_l)
                                        //                = 2 pi r dr dz
                                        vol = 2 * pi * p[0] * dx[0] * dx[1];
                                    } else if (coord == 2) {
                                        // 1-d spherical V = 4/3 pi (r_r**3 - r_l**3)
                                        Real r_r = problo[0]+static_cast<Real>(i+1)*dx[0];
                                        Real r_l = problo[0]+static_cast<Real>(i)*dx[0];
                                        vol = (4.0_rt/3.0_rt) * pi * dx[0] * (r_r*r_r + r_l*r_r + r_l*r_l);
                                    }

                                    lsum += fab(i,j,k) * vol;
                                }
                            }
                        }
                    }
                }
            }
        } else {
            const MultiFab& mf = pf.get(ilev, var_name);
            for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.ok()) {
                    const auto& fab = mf.array(mfi);
                    const auto lo = amrex::lbound(bx);
                    const auto hi = amrex::ubound(bx);
                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {
                                Array<Real,AMREX_SPACEDIM> p
                                    = {AMREX_D_DECL(problo[0]+static_cast<Real>(i+0.5)*dx[0],
                                                    problo[1]+static_cast<Real>(j+0.5)*dx[1],
                                                    problo[2]+static_cast<Real>(k+0.5)*dx[2])};

                                    // compute the volume
                                    Real vol = std::numeric_limits<Real>::quiet_NaN();
                                    if (coord == 0) {
                                        // Cartesian
                                        vol = 1.0_rt;
                                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                                            vol *= dx[idim];
                                        }
                                    } else if (coord == 1) {
                                        // axisymmetric V = pi (r_r**2 - r_l**2) * dz
                                        //                = pi dr * dz * (r_r + r_l)
                                        //                = 2 pi r dr dz
                                        vol = 2 * pi * p[0] * dx[0] * dx[1];
                                    } else if (coord == 2) {
                                        // 1-d spherical V = 4/3 pi (r_r**3 - r_l**3)
                                        Real r_r = problo[0]+static_cast<Real>(i+1)*dx[0];
                                        Real r_l = problo[0]+static_cast<Real>(i)*dx[0];
                                        vol = (4.0_rt/3.0_rt) * pi * dx[0] * (r_r*r_r + r_l*r_r + r_l*r_l);
                                    }

                                    lsum += fab(i,j,k) * vol;
                            }
                        }
                    }
                }
            }
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
