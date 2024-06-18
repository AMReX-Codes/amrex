//
// Laterally average a variable with optional density weighting (Favre average)

#include <iostream>
// #include <stringstream>
#include <regex>
#include <string>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

std::string inputs_name = "";

void PrintUsage ();


int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv, false);

    {

        // timer for profiling
        BL_PROFILE_VAR("main()", pmain);

        // Input arguments
        std::string varname{"density"};
        bool do_favre{false};

        const int narg = amrex::command_argument_count();

        int farg = 1;
        while (farg <= narg) {
            const std::string& name = amrex::get_command_argument(farg);
            if (name == "-v" || name == "--variable") {
                varname = amrex::get_command_argument(++farg);
            } else if (name == "-f" || name == "--favre") {
                do_favre = true;
            } else {
                break;
            }
            ++farg;
        }

        if (farg > narg) {
            PrintUsage();
            return -1;
        }

        const std::string& pltfile = amrex::get_command_argument(farg);
        std::string slcfile = pltfile + ".slice";

        PlotFileData pf(pltfile);

        int fine_level = pf.finestLevel();
        const int dim = pf.spaceDim();

        // get the index bounds and dx.
        Box domain = pf.probDomain(fine_level);
        auto dx_fine = pf.cellSize(fine_level);
        auto problo = pf.probLo();
        auto probhi = pf.probHi();

        // compute the size of the radially-binned array -- we'll take the
        // vertical direction to be the dimensionality

        int nbins = static_cast<int>(std::abs(probhi[dim-1] - problo[dim-1]) / dx_fine[dim-1]);

        // height coordinate
        Vector<Real> h(nbins);

        for (auto i = 0; i < nbins; i++) {
            h[i] = (i + 0.5) * dx_fine[dim-1];
        }

        // find variable indices
        const Vector<std::string>& var_names_pf = pf.varNames();

        // density can be call "density" in Castro or "rho" in MAESTROeX
        int dens_comp = std::distance(var_names_pf.cbegin(),
                                  std::find(var_names_pf.cbegin(), var_names_pf.cend(), "density"));

        if (dens_comp == var_names_pf.size()) {
            dens_comp = std::distance(var_names_pf.cbegin(),
                                      std::find(var_names_pf.cbegin(), var_names_pf.cend(), "rho"));
        }

        if (dens_comp == var_names_pf.size() && do_favre) {
            amrex::Error("density not found");
        }

        int var_comp = std::distance(var_names_pf.cbegin(),
                                  std::find(var_names_pf.cbegin(), var_names_pf.cend(), varname));

        if (var_comp == var_names_pf.size()) {
            amrex::Error("variable " + varname + " not found");
        }


        // allocate storage for data
        Vector<Real> var_bin(nbins, 0.);
        Vector<Real> volcount(nbins, 0.);

        // we will use a mask that tells us if a zone on the current level
        // is covered by data on a finer level.

        for (int ilev = 0; ilev <= fine_level; ++ilev) {

            Array<Real, AMREX_SPACEDIM> dx_level = pf.cellSize(ilev);
            Real vol = dx_level[0];
            if (dim >= 2) {
                vol *= dx_level[1];
            }
            if (dim == 3) {
                vol *= dx_level[2];
            }

            if (ilev < fine_level) {
                IntVect ratio{pf.refRatio(ilev)};
                for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                    ratio[idim] = 1;
                }
                const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                    pf.boxArray(ilev+1), ratio);

                const MultiFab& lev_data_mf = pf.get(ilev);

                for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    if (bx.ok()) {
                        const auto& m = mask.array(mfi);
                        const auto& fab = lev_data_mf.array(mfi);
                        const auto lo = amrex::lbound(bx);
                        const auto hi = amrex::ubound(bx);

                        for (int k = lo.z; k <= hi.z; ++k) {
                            for (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {
                                    if (m(i,j,k) == 0) { // not covered by fine

                                        Real height{0.0};
                                        if (dim == 1) {
                                            height = problo[0] + static_cast<Real>(i+0.5)*dx_level[0];
                                        } else if (dim == 2) {
                                            height = problo[1] + static_cast<Real>(j+0.5)*dx_level[1];
                                        } else {
                                            height = problo[2] + static_cast<Real>(k+0.5)*dx_level[2];
                                        }

                                        int index = static_cast<int>(height / dx_fine[dim-1]);

                                        // add to the bin, weighting by the size

                                        if (do_favre) {
                                            var_bin[index] += fab(i,j,k,dens_comp) * fab(i,j,k,var_comp) * vol;
                                            volcount[index] += fab(i,j,k,dens_comp) * vol;
                                        } else {
                                            var_bin[index] += fab(i,j,k,var_comp) * vol;
                                            volcount[index] += vol;
                                        }

                                    } // mask

                                }
                            }
                        }

                    } // bx.ok()

                } // MFIter

            } else {
                // this is the finest level

                const MultiFab& lev_data_mf = pf.get(ilev);

                for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    if (bx.ok()) {
                        const auto& fab = lev_data_mf.array(mfi);
                        const auto lo = amrex::lbound(bx);
                        const auto hi = amrex::ubound(bx);

                        for (int k = lo.z; k <= hi.z; ++k) {
                            for (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {

                                    Real height{0.0};
                                    if (dim == 1) {
                                        height = problo[0] + static_cast<Real>(i+0.5)*dx_level[0];
                                    } else if (dim == 2) {
                                        height = problo[1] + static_cast<Real>(j+0.5)*dx_level[1];
                                    } else {
                                        height = problo[2] + static_cast<Real>(k+0.5)*dx_level[2];
                                    }

                                    int index = static_cast<int>(height / dx_fine[dim-1]);

                                    // add to the bin, weighting by the size

                                    if (do_favre) {
                                        var_bin[index] += fab(i,j,k,dens_comp) * fab(i,j,k,var_comp) * vol;
                                        volcount[index] += fab(i,j,k,dens_comp) * vol;
                                    } else {
                                        var_bin[index] += fab(i,j,k,var_comp) * vol;
                                        volcount[index] += vol;
                                    }

                                }
                            }
                        }

                    } // bx.ok()

                } // MFIter


            }

        } // level loop


        ParallelDescriptor::ReduceRealSum(var_bin.data(), static_cast<int>(var_bin.size()));
        ParallelDescriptor::ReduceRealSum(volcount.data(), static_cast<int>(volcount.size()));

        //normalize

        for (int i = 0; i < nbins; i++) {
            if (volcount[i] != 0.0) {
                var_bin[i] /= volcount[i];
            }
        }

        // now open the slicefile and write out the data
        amrex::Print() << "outputting lateral average of " << varname << " to " << slcfile << '\n';

        std::ofstream slicefile;
        slicefile.open(slcfile);
        slicefile.setf(std::ios::scientific);
        slicefile.precision(12);
        const auto w = 24;

        // write the header
        slicefile << "# " << std::setw(w) << "height" << std::setw(w) << varname << '\n';

        // write the data in columns
        const auto SMALL = 1.e-20;
        for (auto i = 0; i < nbins; i++) {
            if (std::abs(var_bin[i]) < SMALL) var_bin[i] = 0.0;

            slicefile << std::setw(w) << h[i] << std::setw(w) << var_bin[i] << '\n';
        }

        slicefile.close();

        // destroy timer for profiling
        BL_PROFILE_VAR_STOP(pmain);

    }

    amrex::Finalize();
}



//
// Print usage info
//
void PrintUsage ()
{
    amrex::Print() << "\n"
                   << " Horizontally average a variable\n"
                   << " Usage: \n"
                   << "     faverage {[-v|--variable] name} {[-f|--favre]} plotfile\n"
                   << "\n"
                   << "    -v name    : variable to average (default: density)\n"
                   << "    -f         : do Favre (density-weighted) average\n"
                   << '\n';
}
