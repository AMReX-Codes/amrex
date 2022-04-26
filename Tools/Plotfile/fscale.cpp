#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

using namespace amrex;

void main_main() {
  bool b_scale = false;

  const int narg = amrex::command_argument_count();
  Real length_unit = NAN;

  int farg = 1;
  while (farg <= narg) {
    const auto fname = get_command_argument(farg);
    if (fname == "-s" || fname == "--scale") {
      b_scale = true;
      ++farg;
      length_unit = std::stod(amrex::get_command_argument(farg));
    } else {
      break;
    }
    ++farg;
  }

  if (!b_scale) {
    amrex::Abort("ERROR: length scale must be specified!");
  }

  if (std::isnan(length_unit)) {
    amrex::Abort("ERROR: length scale is invalid!");
  }

  const auto &fname = amrex::get_command_argument(farg);
  ++farg;
  const auto &fname_output = amrex::get_command_argument(farg);

  if (farg > narg) {
    amrex::Print()
        << "\n"
        << " Scale the RealBox dimensions of a plotfile.\n"
        << " Works with 1-, 2-, or 3-d datasets.\n"
        << "\n"
        << " usage:\n"
        << "    fboxinfo -s|--scale length_unit plotfile output_plotfile\n"
        << "\n"
        << " args:\n"
        << "    -s|--scale     divide RealBox coordinates by this value\n"
        << std::endl;
    return;
  }

  PlotFileData plotfile(fname);
  amrex::Print() << " reading plotfile: " << fname << "\n";

  const int nlevels = plotfile.finestLevel() + 1;
  const Vector<std::string> &varnames = plotfile.varNames();
  Vector<MultiFab> mf;
  Vector<Geometry> geom;
  Vector<int> level_step;
  Vector<IntVect> ref_ratio;

  for (int ilev = 0; ilev < nlevels; ++ilev) {
    // is this readable from the plotfile??
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};

    // can this be anisotropic??
    const int ref = plotfile.refRatio(ilev);

    mf.push_back(plotfile.get(ilev));
    ref_ratio.push_back(IntVect{AMREX_D_DECL(ref, ref, ref)});
    level_step.push_back(plotfile.levelStep(ilev));

    // rescale geometry
    const Box prob_domain = plotfile.probDomain(ilev);
    auto const &dlo = plotfile.probLo();
    auto const &dhi = plotfile.probHi();
    std::array<amrex::Real, AMREX_SPACEDIM> new_dlo;
    std::array<amrex::Real, AMREX_SPACEDIM> new_dhi;

    for (int k = 0; k < AMREX_SPACEDIM; ++k) {
      new_dlo[k] = dlo[k] / length_unit;
      new_dhi[k] = dhi[k] / length_unit;
    }
    amrex::RealBox rescaledRealBox(new_dlo, new_dhi);
    const Geometry geom_l(prob_domain, rescaledRealBox, plotfile.coordSys(),
                          is_periodic);
    geom.push_back(geom_l);
  }

  // save to new plotfile
  amrex::Print() << " writing to new plotfile: " << fname_output << "\n";

  const Vector<const MultiFab *> vec_mf = GetVecOfConstPtrs(mf);
  WriteMultiLevelPlotfileHeaders(fname_output, nlevels, vec_mf, varnames, geom,
                                 plotfile.time(), level_step, ref_ratio);

  amrex::Print() << std::endl;
}

int main(int argc, char *argv[]) {
  amrex::SetVerbose(0);
  amrex::Initialize(argc, argv, false);
  main_main();
  amrex::Finalize();
}
