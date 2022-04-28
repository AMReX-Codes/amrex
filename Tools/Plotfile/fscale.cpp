#include "AMReX_ParallelDescriptor.H"
#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <fstream>
#include <sstream>

using namespace amrex;

namespace {
void GotoNextLine(std::istream &is) {
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}
} // namespace

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

  amrex::Print() << " reading plotfile header: " << fname << "\n";

  // read header file from disk
  std::string File(fname + "/Header");
  Vector<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
  std::istringstream is(std::string(fileCharPtr.dataPtr()),
                        std::istringstream::in);

  // stringstream for output header
  std::ostringstream os;

  // parse input file line-by-line
  std::string m_file_version;
  is >> m_file_version;
  os << m_file_version << '\n';

  int m_ncomp;
  is >> m_ncomp;
  os << m_ncomp << '\n';

  Vector<std::string> m_var_names;
  m_var_names.resize(m_ncomp);
  for (int i = 0; i < m_ncomp; ++i) {
    is >> m_var_names[i];
    os << m_var_names[i] << '\n';
  }

  int m_spacedim;
  is >> m_spacedim;
  os << m_spacedim << '\n';

  Real m_time;
  is >> m_time;
  os << std::setprecision(30) << m_time << '\n';

  int m_finest_level, m_nlevels;
  is >> m_finest_level;
  os << m_finest_level << '\n';
  m_nlevels = m_finest_level + 1;

  Array<Real, AMREX_SPACEDIM> m_prob_lo{{AMREX_D_DECL(0., 0., 0.)}};
  for (int i = 0; i < m_spacedim; ++i) {
    is >> m_prob_lo[i];
    os << (m_prob_lo[i] / length_unit) << ' ';
  }
  os << '\n';

  Array<Real, AMREX_SPACEDIM> m_prob_hi{{AMREX_D_DECL(1., 1., 1.)}};
  for (int i = 0; i < m_spacedim; ++i) {
    is >> m_prob_hi[i];
    os << (m_prob_hi[i] / length_unit) << ' ';
  }
  os << '\n';

  Vector<int> m_ref_ratio;
  m_ref_ratio.resize(m_nlevels, 0);
  for (int i = 0; i < m_finest_level; ++i) {
    is >> m_ref_ratio[i];
    os << m_ref_ratio[i] << ' ';
  }
  os << '\n';

  GotoNextLine(is);

  Vector<Box> m_prob_domain;
  m_prob_domain.resize(m_nlevels);
  for (int i = 0; i < m_nlevels; ++i) {
    is >> m_prob_domain[i];
    os << m_prob_domain[i] << ' ';
  }
  os << '\n';

  Vector<int> m_level_steps;
  m_level_steps.resize(m_nlevels);
  for (int i = 0; i < m_nlevels; ++i) {
    is >> m_level_steps[i];
    os << m_level_steps[i] << ' ';
  }
  os << '\n';

  Vector<Array<Real, AMREX_SPACEDIM>> m_cell_size;
  m_cell_size.resize(m_nlevels,
                     Array<Real, AMREX_SPACEDIM>{{AMREX_D_DECL(1., 1., 1.)}});
  for (int ilev = 0; ilev < m_nlevels; ++ilev) {
    for (int idim = 0; idim < m_spacedim; ++idim) {
      is >> m_cell_size[ilev][idim];
      os << (m_cell_size[ilev][idim] / length_unit) << ' ';
    }
    os << '\n';
  }

  GotoNextLine(is);

  // process remaining lines and write out until EOF
  for (std::string line; std::getline(is, line);) {
    os << line << '\n';
  }

  // write ostringstream 'os' to file
  if (ParallelDescriptor::IOProcessor()) {
    amrex::Print() << " writing new header to plotfile: " << fname_output << "\n";
    std::ofstream outfile;
    outfile.open(fname_output + "/Header");
    outfile << os.str();
    outfile.close();
  }

  amrex::Print() << std::endl;
}

int main(int argc, char *argv[]) {
  amrex::SetVerbose(0);
  amrex::Initialize(argc, argv, false);
  main_main();
  amrex::Finalize();
  return 0;
}
