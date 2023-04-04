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
    Real xcoord = std::numeric_limits<Real>::lowest();
    Real ycoord = std::numeric_limits<Real>::lowest();
    Real zcoord = std::numeric_limits<Real>::lowest();
    bool scientific = false;
    bool csv = false;
    int  precision = 17;
    Real tolerance = std::numeric_limits<Real>::lowest();
    bool print_info = false;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-s" || name == "--slicefile") {
            slicefile = amrex::get_command_argument(++farg);
        } else if (name == "-d" || name == "--direction") {
            idir = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-v" || name == "--variable") {
            varnames_arg = amrex::get_command_argument(++farg);
        } else if (name == "-x") {
            xcoord = Real(std::stod(amrex::get_command_argument(++farg)));
        } else if (name == "-y") {
            ycoord = Real(std::stod(amrex::get_command_argument(++farg)));
        } else if (name == "-z") {
            zcoord = Real(std::stod(amrex::get_command_argument(++farg)));
        } else if (name == "-l" || name == "--lower_left") {
            center = false;
        } else if (name == "-c" || name == "--coarse_level") {
            coarse_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-f" || name == "--fine_level") {
            fine_level = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-e" || name == "--scientific") {
            scientific = true;
        } else if (name == "-csv" || name == "--csv") {
            csv = true;
        } else if (name == "-p" || name == "--precision") {
            precision = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-t" || name == "--tolerance") {
            tolerance = Real(std::stod(amrex::get_command_argument(++farg)));
        } else if (name == "-i" || name == "--info") {
            print_info = true;
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
            << "      [-x][-y][-z]                     : (x,y,z)-coordinate to pass through\n"
            << "                                         (overrides center/lower-left)\n"
            << "      [-c|--coarse_level] coarse level : coarsest level to extract from\n"
            << "      [-f|--fine_level]   fine level   : finest level to extract from\n"
            << "      [-e|--scientific]                : output data in scientific notation\n"
            << "      [-p|--precision]    precision    : decimal precision {17 (default)}\n"
            << "      [-t|--tolerance]    tolerance    : set to 0 any value lower than tolerance\n"
            << "      [-i|--info]                      : output job info\n"
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

    if (xcoord > -1.e36 && AMREX_SPACEDIM >= 1) {
        // we specified the x value to pass through
        iloc = hi0.x;
        for (int i = lo0.x; i <= hi0.x; ++i) {
            amrex::Real xc = problo[0] + (Real(i)+Real(0.5))*dx0[0];
            if (xc > xcoord) {
                iloc = i;
                break;
            }
        }
    }

    if (ycoord > -1.e36 && AMREX_SPACEDIM >= 2) {
        // we specified the y value to pass through
        jloc = hi0.y;
        for (int j = lo0.y; j <= hi0.y; ++j) {
            amrex::Real yc = problo[1] + (Real(j)+Real(0.5))*dx0[1];
            if (yc > ycoord) {
                jloc = j;
                break;
            }
        }
    }

    if (zcoord > -1.e36 && AMREX_SPACEDIM == 3) {
        // we specified the z value to pass through
        kloc = hi0.z;
        for (int k = lo0.z; k <= hi0.z; ++k) {
            amrex::Real zc = problo[2] + (Real(k)+Real(0.5))*dx0[2];
            if (zc > zcoord) {
                kloc = k;
                break;
            }
        }
    }

    const int dim = pf.spaceDim();

    if (idir < 0 || idir >= dim) {
        amrex::Print() << " invalid direction\n";
        return;
    } else if (idir == 0) {
        amrex::Print() << " slicing along x-direction at coarse grid (j,k)=(" << jloc << "," << kloc << ") and output to " << slicefile << "\n";
    } else if (idir == 1) {
        amrex::Print() << " slicing along y-direction at coarse grid (i,k)=(" << iloc << "," << kloc << ") and output to " << slicefile << "\n";
    } else if (idir == 2) {
        amrex::Print() << " slicing along z-direction at coarse grid (i,j)=(" << iloc << "," << jloc << ") and output to " << slicefile << "\n";
    }


    const IntVect ivloc{AMREX_D_DECL(iloc,jloc,kloc)};

    if (fine_level < 0) fine_level = pf.finestLevel();
    // sanity check on valid selected levels
    if (fine_level > pf.finestLevel() || coarse_level < 0 || coarse_level > fine_level) {
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
                                                = {AMREX_D_DECL(problo[0]+static_cast<Real>(i+0.5)*dx[0],
                                                                problo[1]+static_cast<Real>(j+0.5)*dx[1],
                                                                problo[2]+static_cast<Real>(k+0.5)*dx[2])};
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
                                            = {AMREX_D_DECL(problo[0]+static_cast<Real>(i+0.5)*dx[0],
                                                            problo[1]+static_cast<Real>(j+0.5)*dx[1],
                                                            problo[2]+static_cast<Real>(k+0.5)*dx[2])};
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
        const auto numpts = int(pos.size());
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
            for (int i = 0, N = int(numpts_vec.size()); i < N; ++i) {
                ntot += numpts_vec[i];
                recvcnt[i] = numpts_vec[i];
                if (i+1 < N) {
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

        if (scientific) {
          ofs << std::scientific;
        }

        if (csv)
        {
           ofs << std::setw(24) << std::left << dirstr;
           for (auto const& vname : var_names) {
              ofs << ", " <<  std::setw(24) << std::left << vname;
           }
           ofs << "\n";
           for (int i = 0; i < posidx.size(); ++i) {
              ofs << std::setw(25) << std::left << std::setprecision(17) << posidx[i].first;
              for (int j = 0; j < var_names.size(); ++j) {
                  ofs << ", " << std::setw(25) << std::left << std::setprecision(17) << data[j][posidx[i].second];
              }
              ofs << "\n";
           }
        }
        else
        {
           ofs << "# 1-d slice in " << dirstr << "-direction, file: " << pltfile << "\n";
           ofs << "# time = " << std::setw(20) << std::setprecision(precision) << pf.time() << "\n";

           ofs << "#" << std::setw(24) << dirstr;
           for (auto const& vname : var_names) {
             ofs << " " <<  std::setw(24) << std::right << vname;
           }
           ofs << "\n";

           for (int i = 0; i < posidx.size(); ++i) {
             ofs << std::setw(25) << std::right << std::setprecision(precision) << posidx[i].first;
             for (int j = 0; j < var_names.size(); ++j) {
               if (std::abs(data[j][posidx[i].second])< tolerance ) data[j][posidx[i].second] = 0.;
               ofs << std::setw(25) << std::right << std::setprecision(precision) << data[j][posidx[i].second];
             }
             ofs << "\n";
           }

           if (print_info) {
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
    }
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
