#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string slicefile;
    std::string pltfile;
    int idir = 0;
    std::string varnames_arg;
    bool center = true;
    int coarse_level = 0;
    int fine_level = -1;  // This will be fixed later
    Real ycoord = std::numeric_limits<Real>::lowest();
    bool scientific = false;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-s" or name == "--slicefile") {
            slicefile = amrex::get_command_argument(++farg);
        } else if (name == "-d" or name == "--direction") {
            idir = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-v" or name == "--variable") {
            varnames_arg = amrex::get_command_argument(++farg);
        } else if (name == "-y") {
            ycoord = std::stod(amrex::get_command_argument(++farg));
        } else if (name == "-l" or name == "--lower_left") {
            center = false;
        } else if (name == "-c" or name == "--coarse_level") {
            coarse_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-f" or name == "--fine_level") {
            fine_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-e" or name == "--scientific") {
            scientific = true;
        } else {
            break;
        }
        ++farg;
    }

    if (pltfile.empty() and farg <= narg) {
        pltfile = amrex::get_command_argument(farg);
    }

    if (pltfile.empty()) {
        amrex::Print()
            << "\n"
            << " Extract at 1D slice through a plotfile in any coordinate direction.\n"
            << " Works with 1-, 2-, or 3-d datasets.\n"
            << "\n"
            << " Usage:\n"
            << "    fextract [-s outfile] [-d dir] [-v variable] [-y ycood] [-c coarse_level] [-f fine_level] plotfile\n"
            << "\n"
            << " args [-s|--slicefile] slice file      : slice file          (optional)\n"
            << "      [-d|--direction] idir            : slice direction {0 (default), 1, or 2}\n"
            << "      [-v|--variable]  varname(s)      : only output the values of variable\n"
            << "                                         varname (space separated string for\n"
            << "                                         multiple variables)\n"
            << "      [-l|--lower_left]                : slice through lower left corner\n"
            << "                                         instead of center\n"
            << "      [-y]                             : y-coordinate to pass through\n"
            << "                                         (overrides center/lower-left)\n"
            << "      [-c|--coarse_level] coarse level : coarsest level to extract from\n"
            << "      [-f|--fine_level]   fine level   : finest level to extract from\n"
            << "      [-e|--scientific]                : output data in scientific notation\n"
            << "\n"
            << " If a job_info file is present in the plotfile, that information is made\n"
            << " available at the end of the slice file (commented out), for reference.\n"
            << std::endl;
        return;
    }

    // slicefile not defined, default to plotfile.slice
    if (slicefile.empty()) {
        slicefile = pltfile;
        if (slicefile.back() == '/') {
            slicefile.pop_back();
        }
        slicefile += ".slice";
    }

    PlotFileData pf(pltfile);
    const int dim = pf.spaceDim();

    if (idir < 0 or idir >= dim) {
        amrex::Print() << " invalid direction\n";
        return;
    } else if (idir == 0) {
        amrex::Print() << " slicing along x-direction and output to " << slicefile << "\n";
    } else if (idir == 1) {
        amrex::Print() << " slicing along y-direction and output to " << slicefile << "\n";
    } else if (idir == 2) {
        amrex::Print() << " slicing along z-direction and output to " << slicefile << "\n";
    }

    const Vector<std::string>& var_names_pf = pf.varNames();

    Vector<std::string> var_names;
    if (varnames_arg.empty()) {
        var_names = var_names_pf;
    } else {
        std::istringstream is(varnames_arg);
        var_names.assign(std::istream_iterator<std::string>{is},
                         std::istream_iterator<std::string>{  });
        var_names.erase(std::remove_if(var_names.begin(), var_names.end(),
                                       [&](std::string const& x) {
                                           return var_names_pf.end() ==
                                               std::find(var_names_pf.begin(), var_names_pf.end(), x);
                                       }),
                        var_names.end());
        if (var_names.empty()) {
            amrex::Abort("ERROR: no valid variable names");
        }
    }

    Array<Real,AMREX_SPACEDIM> problo = pf.probLo();
    Array<Real,AMREX_SPACEDIM> dx0 = pf.cellSize(0);
    Box probdom0 = pf.probDomain(0);
    const auto lo0 = amrex::lbound(probdom0);
    const auto hi0 = amrex::ubound(probdom0);

    // compute the index of the center or lower left of the domain on the
    // coarse grid.  These are used to set the position of the slice in
    // the transverse direction.

    int iloc = 0, jloc = 0, kloc = 0;
    if (center) {
        iloc = (hi0.x-lo0.x+1)/2 + lo0.x;
        jloc = (hi0.y-lo0.y+1)/2 + lo0.y;
        kloc = (hi0.z-lo0.z+1)/2 + lo0.z;
    }

    if (ycoord > -1.e-36 and AMREX_SPACEDIM >= 2) {
        // we specified the y value to pass through
        for (int j = lo0.y; j <= hi0.y; ++j) {
            amrex::Real yc = problo[1] + (j+0.5)*dx0[1];
            if (yc > ycoord) {
                jloc = j;
                break;
            }
        }
    }

    const IntVect ivloc{AMREX_D_DECL(iloc,jloc,kloc)};

    if (fine_level < 0) fine_level = pf.finestLevel();
    // sanity check on valid selected levels
    if (fine_level > pf.finestLevel() or coarse_level < 0 or coarse_level > fine_level) {
        amrex::Abort("Invalid level selection");
    }

    Vector<Real> pos;
    Vector<Vector<Real> > data(var_names.size());

    IntVect rr{1};
    for (int ilev = coarse_level; ilev <= fine_level; ++ilev) {
        Box slice_box(ivloc*rr,ivloc*rr);
        slice_box.setSmall(idir, std::numeric_limits<int>::lowest());
        slice_box.setBig(idir, std::numeric_limits<int>::max());

        Array<Real,AMREX_SPACEDIM> dx = pf.cellSize(ilev);

        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                pf.boxArray(ilev+1), ratio);
            for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                const MultiFab& mf = pf.get(ilev, var_names[ivar]);
                for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox() & slice_box;
                    if (bx.ok()) {
                        const auto& m = mask.array(mfi);
                        const auto& fab = mf.array(mfi);
                        const auto lo = amrex::lbound(bx);
                        const auto hi = amrex::ubound(bx);
                        for         (int k = lo.z; k <= hi.z; ++k) {
                            for     (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {
                                    if (m(i,j,k) == 0) { // not covered by fine
                                        if (pos.size() == data[ivar].size()) {
                                            Array<Real,AMREX_SPACEDIM> p
                                                = {AMREX_D_DECL(problo[0]+(i+0.5)*dx[0],
                                                                problo[1]+(j+0.5)*dx[1],
                                                                problo[2]+(k+0.5)*dx[2])};
                                            pos.push_back(p[idir]);
                                        }
                                        data[ivar].push_back(fab(i,j,k));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            rr *= ratio;
        } else {
            for (int ivar = 0; ivar < var_names.size(); ++ivar) {
                const MultiFab& mf = pf.get(ilev, var_names[ivar]);
                for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox() & slice_box;
                    if (bx.ok()) {
                        const auto& fab = mf.array(mfi);
                        const auto lo = amrex::lbound(bx);
                        const auto hi = amrex::ubound(bx);
                        for         (int k = lo.z; k <= hi.z; ++k) {
                            for     (int j = lo.y; j <= hi.y; ++j) {
                                for (int i = lo.x; i <= hi.x; ++i) {
                                    if (pos.size() == data[ivar].size()) {
                                        Array<Real,AMREX_SPACEDIM> p
                                            = {AMREX_D_DECL(problo[0]+(i+0.5)*dx[0],
                                                            problo[1]+(j+0.5)*dx[1],
                                                            problo[2]+(k+0.5)*dx[2])};
                                        pos.push_back(p[idir]);
                                    }
                                    data[ivar].push_back(fab(i,j,k));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

#ifdef BL_USE_MPI
    {
        const int numpts = pos.size();
        auto numpts_vec = ParallelDescriptor::Gather(numpts,
                                                     ParallelDescriptor::IOProcessorNumber());
        Vector<int> recvcnt, disp;
        Vector<Real> allpos;
        Vector<Vector<Real> > alldata(data.size());
        if (ParallelDescriptor::IOProcessor()) {
            recvcnt.resize(numpts_vec.size());
            disp.resize(numpts_vec.size());
            int ntot = 0;
            disp[0] = 0;
            for (int i = 0; i < numpts_vec.size(); ++i) {
                ntot += numpts_vec[i];
                recvcnt[i] = numpts_vec[i];
                if (i+1 < numpts_vec.size()) {
                    disp[i+1] = disp[i] + numpts_vec[i];
                }
            }
            allpos.resize(ntot);
            alldata.resize(data.size());
            for (auto& v : alldata) {
                v.resize(ntot);
            }
        } else {
            recvcnt.resize(1);
            disp.resize(1);
            allpos.resize(1);
            for (auto& v: alldata) {
                v.resize(1);
            }
        }
        ParallelDescriptor::Gatherv(pos.data(), numpts, allpos.data(), recvcnt, disp,
                                    ParallelDescriptor::IOProcessorNumber());
        for (int i = 0; i < data.size(); ++i) {
            ParallelDescriptor::Gatherv(data[i].data(), numpts, alldata[i].data(), recvcnt, disp,
                                        ParallelDescriptor::IOProcessorNumber());
        }
        if (ParallelDescriptor::IOProcessor()) {
            pos = std::move(allpos);
            data = std::move(alldata);
        }
    }
#endif

    if (ParallelDescriptor::IOProcessor()) {
        Vector<std::pair<Real,int> > posidx;
        posidx.reserve(pos.size());
        for (int i = 0; i < pos.size(); ++i) {
            posidx.emplace_back(pos[i],i);
        }
        std::sort(posidx.begin(), posidx.end());

        const std::string dirstr = (idir == 0) ? "x" : ((idir == 1) ? "y" : "z");

        std::ofstream ofs(slicefile, std::ios::trunc);

        ofs << "# 1-d slice in " << dirstr << "-direction, file: " << pltfile << "\n";
        ofs << "# time = " << std::setprecision(17) << pf.time() << "\n";

        ofs << "#" << std::setw(24) << dirstr;
        for (auto const& vname : var_names) {
            ofs << " " <<  std::setw(24) << std::right << vname;
        }
        ofs << "\n";
        if (scientific) {
            ofs << std::scientific;
        }
        for (int i = 0; i < posidx.size(); ++i) {
            ofs << std::setw(25) << std::right << std::setprecision(12) << posidx[i].first;
            for (int j = 0; j < var_names.size(); ++j) {
                ofs << std::setw(25) << std::right << std::setprecision(12) << data[j][posidx[i].second];
            }
            ofs << "\n";
        }

        // job_info? if so write it out to the slice file end
        std::ifstream jobinfo(pltfile+"/job_info");
        if (jobinfo.good()) {
            ofs << "\n";
            std::string s;
            while (std::getline(jobinfo, s)) {
                ofs << "#" << s << "\n";
            }
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
