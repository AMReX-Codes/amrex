#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>
#include <cstdlib>
#include <iterator>
#include <sstream>
#include <string>

using namespace amrex;

void PrintUsage ()
{
    amrex::Print()
        << "\n"
        << " Compute gradient of variables in a given plofile. The results are saved in a\n"
        << " new plotfile.\n"
        << "\n"
        << " Usage:\n"
        << "    fgradient [-o outfile] [-v varnames] [-xlo xlo_bc] [-ylo ylo_bc]\n"
        << "              [-zlo zlo_bc] [-xhi xhi_bc] [-yhi yhi_bc] [-zhi zhi_bc]\n"
        << "              [-i interp] plotfile\n"
        << "\n"
        << " optional arguments:\n"
        << "    -o|--output outfile\n"
        << "            Gradient plotfile name. If it's not provided, it is grad.PLOTFILE,\n"
        << "            where PLOTFILE is the original plotfile name.\n"
        << "\n"
        << "    -v|--variables varnames\n"
        << "            Variable names separated by spaces. If there are multiple names,\n"
        << "            they must be inside a pair of double quotes. If not provided, all\n"
        << "            variables in the plotfile are included.\n"
        << "\n"
        << "    -xlo|--xlo_boundary xlo_bc\n"
        << "            Domain boundary conditions at x-lo. Choices are\n"
        << "                reflect_odd  : reflecting odd\n"
        << "                reflect_even : reflecting even\n"
        << "                foextrap     : first-order extrapolation\n"
        << "                hoextrap     : second-order extrapolation.\n"
        << "                periodic     : periodic boundary.\n"
        << "            One could provide either different BC types for different\n"
        << "            variables or just one type for all variables. If not provided, the\n"
        << "            default is hoextrap, equivalent to one-sided derivative.\n"
        << "\n"
        << "    -ylo|--ylo_boundary ylo_bc\n"
        << "            Domain boundary conditions at y-lo.\n"
        << "\n"
        << "    -zlo|--zlo_boundary zlo_bc\n"
        << "            Domain boundary conditions at z-lo.\n"
        << "\n"
        << "    -xhi|--xhi_boundary xhi_bc\n"
        << "            Domain boundary conditions at x-hi.\n"
        << "\n"
        << "    -yhi|--yhi_boundary yhi_bc\n"
        << "            Domain boundary conditions at y-hi.\n"
        << "\n"
        << "    -zhi|--zhi_boundary zhi_bc\n"
        << "            Domain boundary conditions at z-hi.\n"
        << "\n"
        << "    -i|--interpolater interp\n"
        << "            Interpolation from coarse to fine level. Choices are 0 (default)\n"
        << "            or 1, where 0 is amrex::cell_cons_interp (conservative linear)\n"
        << "            and 1 is amrex::pc_interp (piecewise constant).\n"
        << "\n"
        << "    -h|--help\n"
        << "            Print usage.\n"
        << "\n";
}

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile;
    std::string outfile;
    std::string varnames_arg;
    std::string xlobc_arg;
    std::string xhibc_arg;
    std::string ylobc_arg;
    std::string yhibc_arg;
    std::string zlobc_arg;
    std::string zhibc_arg;
    int interpolater = 0;
    bool print_help = false;

    int farg = 1;
    while (farg <= narg) {
        const std::string& name = amrex::get_command_argument(farg);
        if (name == "-o" || name == "--output") {
            outfile = amrex::get_command_argument(++farg);
        } else if (name == "-v" || name == "--variables") {
            varnames_arg = amrex::get_command_argument(++farg);
        } else if (name == "-xlo" || name == "--xlo_boundary") {
            xlobc_arg = amrex::get_command_argument(++farg);
        } else if (name == "-xhi" || name == "--xhi_boundary") {
            xhibc_arg = amrex::get_command_argument(++farg);
        } else if (name == "-ylo" || name == "--ylo_boundary") {
            ylobc_arg = amrex::get_command_argument(++farg);
        } else if (name == "-yhi" || name == "--yhi_boundary") {
            yhibc_arg = amrex::get_command_argument(++farg);
        } else if (name == "-zlo" || name == "--zlo_boundary") {
            zlobc_arg = amrex::get_command_argument(++farg);
        } else if (name == "-zhi" || name == "--zhi_boundary") {
            zhibc_arg = amrex::get_command_argument(++farg);
        } else if (name == "-i" || name == "--interpolater") {
            interpolater = std::stoi(amrex::get_command_argument(++farg));
        } else if (name == "-h" || name == "--help") {
            print_help = true;
        } else {
            break;
        }
        ++farg;
    }

    if (farg <= narg) {
        pltfile = amrex::get_command_argument(farg);
    }

    if (pltfile.empty() || print_help) {
        PrintUsage();
        return;
    }

    if (pltfile.back() == '/') {
        pltfile.pop_back();
    }

    if (outfile.empty()) {
        outfile = "grad."+VisMF::BaseName(pltfile);
    }

    PlotFileData pf(pltfile);

    const int ndims = pf.spaceDim();
    AMREX_ALWAYS_ASSERT(ndims <= AMREX_SPACEDIM);

    const int nlevs = pf.finestLevel() + 1;

    Vector<std::string> varnames;
    if (varnames_arg.empty()) {
        varnames = pf.varNames();
    } else {
        auto const& all_names = pf.varNames();
        std::istringstream is(varnames_arg);
        Vector<std::string> args(std::istream_iterator<std::string>{is},
                                 std::istream_iterator<std::string>{  });
        for (auto const& vn : args) {
            if (std::find(std::begin(all_names),std::end(all_names), vn)
                != std::end(all_names)) {
                varnames.push_back(vn);
            } else {
                amrex::Abort("fgradient: Unknown variable "+vn);
            }
        }
    }

    Vector<std::string> gvarnames;
    for (auto const& name : varnames) {
        gvarnames.push_back(name+"_x");
#if (AMREX_SPACEDIM >= 2)
        if (ndims > 1) {
            gvarnames.push_back(name+"_y");
        }
#endif
#if (AMREX_SPACEDIM == 3)
        if (ndims > 2) {
            gvarnames.push_back(name+"_z");
        }
#endif
    }

    const auto nvars = int(varnames.size());

    BCRec bcr_default;
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    IntVect ng(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (idim < ndims) {
            bcr_default.setLo(idim, BCType::hoextrapcc);
            bcr_default.setHi(idim, BCType::hoextrapcc);
        } else {
            bcr_default.setLo(idim, BCType::int_dir);
            bcr_default.setHi(idim, BCType::int_dir);
            is_periodic[idim] = 1;
            ng[idim] = 0;
        }
    }

    auto f_bcarg = [] (std::string const& bcarg) -> Vector<int>
    {
        Vector<int> r;
        std::istringstream is(bcarg);
        Vector<std::string> args(std::istream_iterator<std::string>{is},
                                 std::istream_iterator<std::string>{  });
        for (auto const& bc : args) {
            if (bc == "reflect_odd") {
                r.push_back(BCType::reflect_odd);
            } else if (bc == "reflect_even") {
                r.push_back(BCType::reflect_even);
            } else if (bc == "foextrap") {
                r.push_back(BCType::foextrap);
            } else if (bc == "hoextrap") {
                r.push_back(BCType::hoextrapcc);
            } else if (bc == "periodic") {
                r.push_back(BCType::int_dir);
            } else {
                amrex::Abort("fgradient: Unknown BC type "+bc);
            }
        }
        return r;
    };

    Array<Vector<int>,3> lobc{f_bcarg(xlobc_arg),
                              f_bcarg(ylobc_arg),
                              f_bcarg(zlobc_arg)};
    Array<Vector<int>,3> hibc{f_bcarg(xhibc_arg),
                              f_bcarg(yhibc_arg),
                              f_bcarg(zhibc_arg)};

    Vector<MultiFab> gmf(nlevs);
    Vector<Geometry> geom;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        gmf[ilev].define(pf.boxArray(ilev), pf.DistributionMap(ilev), nvars*ndims, 0);

        for (int ivar = 0; ivar < nvars; ++ivar)
        {
            Vector<BCRec> bcr{bcr_default};
            auto is_per = is_periodic;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                int bc_type = bcr_default.lo(idim);
                auto nbcs = int(lobc[idim].size());
                if (nbcs == 1) {
                    bc_type = lobc[idim][0];
                } else if (ivar < nbcs) {
                    bc_type = lobc[idim][ivar];
                }
                if (bc_type == BCType::int_dir) {
                    is_per[idim] = 1;
                }
                bcr[0].setLo(idim, bc_type);
                //
                bc_type = bcr_default.hi(idim);
                nbcs = int(hibc[idim].size());
                if (nbcs == 1) {
                    bc_type = hibc[idim][0];
                } else if (ivar < nbcs) {
                    bc_type = hibc[idim][ivar];
                }
                if (bc_type == BCType::int_dir) {
                    AMREX_ALWAYS_ASSERT(is_per[idim] == 1);
                }
                bcr[0].setHi(idim, bc_type);
            }

            Geometry vargeom(pf.probDomain(ilev), RealBox(pf.probLo(),pf.probHi()),
                             pf.coordSys(), is_per);
            if (ivar == 0) {
                geom.push_back(vargeom);
            }

            PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> physbcf
                (vargeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));

            MultiFab mf(pf.boxArray(ilev), pf.DistributionMap(ilev), 1, ng);
            if (ilev == 0) {
                MultiFab smf = pf.get(ilev, varnames[ivar]);
                FillPatchSingleLevel(mf, ng, Real(0.0), {&smf}, {Real(0.0)},
                                     0, 0, 1, vargeom, physbcf, 0);
            } else {
                MultiFab cmf = pf.get(ilev-1, varnames[ivar]);
                MultiFab fmf = pf.get(ilev  , varnames[ivar]);
                Geometry cgeom(pf.probDomain(ilev-1), RealBox(pf.probLo(),pf.probHi()),
                               pf.coordSys(), is_per);
                PhysBCFunct<GpuBndryFuncFab<FabFillNoOp>> cphysbcf
                    (cgeom, bcr, GpuBndryFuncFab<FabFillNoOp>(FabFillNoOp{}));
                auto* mapper = (interpolater == 0)
                                ? (Interpolater*)(&cell_cons_interp)
                                : (Interpolater*)(&pc_interp);
                IntVect ratio(pf.refRatio(ilev-1));
                for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                    ratio[idim] = 1;
                }
                FillPatchTwoLevels(mf, ng, Real(0.0), {&cmf}, {Real(0.0)},
                                   {&fmf}, {Real(0.0)}, 0, 0, 1, cgeom, vargeom,
                                   cphysbcf, 0, physbcf, 0, ratio, mapper, bcr, 0);
            }

            auto const& dx = pf.cellSize(ilev);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& ga = gmf[ilev].array(mfi,ivar*ndims);
                auto const& a = mf.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    ga(i,j,k,0) = (a(i+1,j,k) - a(i-1,j,k)) / (Real(2.0)*dx[0]);
#if (AMREX_SPACEDIM >= 2)
                    if (ndims > 1) {
                        ga(i,j,k,1) = (a(i,j+1,k) - a(i,j-1,k)) / (Real(2.0)*dx[1]);
                    }
#endif
#if (AMREX_SPACEDIM == 3)
                    if (ndims > 2) {
                        ga(i,j,k,2) = (a(i,j,k+1) - a(i,j,k-1)) / (Real(2.0)*dx[2]);
                    }
#endif
                });
            }
        }
    }

    Vector<int> level_steps;
    Vector<IntVect> ref_ratio;
    for (int ilev = 0; ilev < nlevs; ++ilev) {
        level_steps.push_back(pf.levelStep(ilev));
        if (ilev < pf.finestLevel()) {
            ref_ratio.push_back(IntVect(pf.refRatio(ilev)));
            for (int idim = ndims; idim < AMREX_SPACEDIM; ++idim) {
                ref_ratio[ilev][idim] = 1;
            }
        }
    }

    WriteMultiLevelPlotfile(outfile, nlevs, GetVecOfConstPtrs(gmf), gvarnames,
                            geom, pf.time(), level_steps, ref_ratio);
}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
