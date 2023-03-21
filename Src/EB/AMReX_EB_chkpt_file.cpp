#include <AMReX_EB_chkpt_file.H>

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)
#include <utility>

namespace {

const std::string level_prefix = "Level_";

void gotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

}

namespace amrex::EB2 {

// Header information includes the cut and covered boxes (if any)
// Checkpoint file contains data for cut boxes
void
ChkptFile::writeHeader (const BoxArray& cut_ba, const BoxArray& covered_ba,
                        const Geometry& geom,
                        const IntVect& ngrow, bool extend_domain_face,
                        int max_grid_size) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::string HeaderFileName(m_restart_file + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;

        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                std::ofstream::trunc |
                std::ofstream::binary);

        if ( ! HeaderFile.good() )
            FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << "Checkpoint version: 1\n";

        const int nlevels = 1;
        HeaderFile << nlevels << "\n";

        // Geometry
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            HeaderFile << geom.ProbLo(i) << ' ';
        HeaderFile << '\n';

        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            HeaderFile << geom.ProbHi(i) << ' ';
        HeaderFile << '\n';

        // ngrow
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            HeaderFile << ngrow[i] << ' ';
        HeaderFile << '\n';

        // extend domain face
        HeaderFile << extend_domain_face << "\n";

        // max grid size
        HeaderFile << max_grid_size << "\n";

        // BoxArray
        for (int lev = 0; lev < nlevels; ++lev)
        {
            cut_ba.writeOn(HeaderFile);
            HeaderFile << '\n';

            if (! covered_ba.empty()) {
                covered_ba.writeOn(HeaderFile);
                HeaderFile << '\n';
            }
        }
    }
}

void
ChkptFile::writeToFile (const MultiFab& mf, const std::string& mf_name) const
{
    VisMF::Write(mf, MultiFabFileFullPrefix(0, m_restart_file,
                level_prefix, mf_name));
}


ChkptFile::ChkptFile (std::string fname)
    : m_restart_file(std::move(fname))
{}

void
ChkptFile::read_from_chkpt_file (BoxArray& cut_grids, BoxArray& covered_grids,
                                 DistributionMapping& dmap,
                                 MultiFab& volfrac, MultiFab& centroid,
                                 MultiFab& bndryarea, MultiFab& bndrycent,
                                 MultiFab& bndrynorm, Array<MultiFab,AMREX_SPACEDIM>& areafrac,
                                 Array<MultiFab,AMREX_SPACEDIM>& facecent,
                                 Array<MultiFab,AMREX_SPACEDIM>& edgecent,
                                 MultiFab& levelset, int ng_gfab, const Geometry& geom,
                                 const IntVect& ngrow_finest, bool extend_domain_face,
                                 int max_grid_size) const
{
    Real prob_lo[] = {AMREX_D_DECL(std::numeric_limits<Real>::max(),
                                   std::numeric_limits<Real>::max(),
                                   std::numeric_limits<Real>::max())};
    Real prob_hi[] = {AMREX_D_DECL(std::numeric_limits<Real>::lowest(),
                                   std::numeric_limits<Real>::lowest(),
                                   std::numeric_limits<Real>::lowest())};

    std::string File(m_restart_file + "/Header");

    if (amrex::Verbose()) amrex::Print() << "file=" << File << std::endl;

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    std::getline(is, line);

    int nlevs;
    is >> nlevs;
    gotoNextLine(is);
    AMREX_ASSERT(nlevs == 1);

    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
#ifdef AMREX_USE_FLOAT
            prob_lo[i++] = std::stof(word);
#else
            prob_lo[i++] = std::stod(word);
#endif
        }
    }

    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
#ifdef AMREX_USE_FLOAT
            prob_hi[i++] = std::stof(word);
#else
            prob_hi[i++] = std::stod(word);
#endif
        }
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(prob_lo[idim] - geom.ProbLo()[idim]) < std::numeric_limits<Real>::epsilon(),
                                         "EB2::ChkptFile cannot read from a different problem domain");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(prob_hi[idim] - geom.ProbHi()[idim]) < std::numeric_limits<Real>::epsilon(),
                                         "EB2::ChkptFile cannot read from a different problem domain");
    }

    IntVect ngrow_chkptfile;
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            ngrow_chkptfile[i++] = std::stoi(word);
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngrow_chkptfile == ngrow_finest, "EB2::ChkptFile cannot read from different ngrow");

    bool edf_chkptfile;
    is >> edf_chkptfile;
    gotoNextLine(is);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(extend_domain_face == edf_chkptfile,
                                     "EB2::ChkptFile cannot read from different extend_domain_face");

    int mgs_chkptfile;
    is >> mgs_chkptfile;
    gotoNextLine(is);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_grid_size == mgs_chkptfile,
                                     "EB2::ChkptFile cannot read from different max_grid_size");

    if (amrex::Verbose()) amrex::Print() << "Loading cut_grids\n";
    cut_grids.readFrom(is);
    gotoNextLine(is);

    if (is.peek() != EOF) {
        if (amrex::Verbose()) amrex::Print() << "Loading covered_grids\n";
        covered_grids.readFrom(is);
        gotoNextLine(is);
    }

    dmap.define(cut_grids, ParallelDescriptor::NProcs());

    // volfrac
    {
        if (amrex::Verbose()) amrex::Print() << "  Loading " << m_volfrac_name << std::endl;

        volfrac.define(cut_grids, dmap, 1, ng_gfab);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_volfrac_name);
        VisMF::Read(volfrac, prefix);
    }

    // centroid
    {
        if (amrex::Verbose()) amrex::Print() << "  Loading " << m_centroid_name << std::endl;

        centroid.define(cut_grids, dmap, AMREX_SPACEDIM, ng_gfab);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_centroid_name);
        VisMF::Read(centroid, prefix);
    }

    // bndryarea
    {
        if (amrex::Verbose()) amrex::Print() << "  Loading " << m_bndryarea_name << std::endl;

        bndryarea.define(cut_grids, dmap, 1, ng_gfab);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_bndryarea_name);
        VisMF::Read(bndryarea, prefix);
    }

    // bndrycent
    {
        if (amrex::Verbose()) amrex::Print() << "  Loading " << m_bndrycent_name << std::endl;

        bndrycent.define(cut_grids, dmap, AMREX_SPACEDIM, ng_gfab);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_bndrycent_name);
        VisMF::Read(bndrycent, prefix);
    }

    // bndrynorm
    {
        if (amrex::Verbose()) amrex::Print() << "  Loading " << m_bndrynorm_name << std::endl;

        bndrynorm.define(cut_grids, dmap, AMREX_SPACEDIM, ng_gfab);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_bndrynorm_name);
        VisMF::Read(bndrynorm, prefix);
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // areafrac
        {
            if (amrex::Verbose()) amrex::Print() << "  Loading " << m_areafrac_name[idim] << std::endl;

            areafrac[idim].define(convert(cut_grids, IntVect::TheDimensionVector(idim)), dmap, 1, ng_gfab);

            auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_areafrac_name[idim]);
            VisMF::Read(areafrac[idim], prefix);
        }

        // facecent
        {
            if (amrex::Verbose()) amrex::Print() << "  Loading " << m_facecent_name[idim] << std::endl;

            facecent[idim].define(convert(cut_grids, IntVect::TheDimensionVector(idim)), dmap, AMREX_SPACEDIM-1, ng_gfab);

            auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_facecent_name[idim]);
            VisMF::Read(facecent[idim], prefix);
        }

        // edgecent
        {
            if (amrex::Verbose()) amrex::Print() << "  Loading " << m_edgecent_name[idim] << std::endl;

            IntVect edge_type{1}; edge_type[idim] = 0;
            edgecent[idim].define(convert(cut_grids, edge_type), dmap, 1, ng_gfab);

            auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_edgecent_name[idim]);
            VisMF::Read(edgecent[idim], prefix);
        }
    }

    // levelset
    {
        if (amrex::Verbose()) amrex::Print() << "  Loading " << m_levelset_name << std::endl;

        levelset.define(convert(cut_grids,IntVect::TheNodeVector()), dmap, 1, ng_gfab);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_levelset_name);
        VisMF::Read(levelset, prefix);
    }
}

void
ChkptFile::write_to_chkpt_file (const BoxArray& cut_grids,
                                const BoxArray& covered_grids,
                                const MultiFab& volfrac,
                                const MultiFab& centroid, const MultiFab& bndryarea,
                                const MultiFab& bndrycent, const MultiFab& bndrynorm,
                                const Array<MultiFab,AMREX_SPACEDIM>& areafrac,
                                const Array<MultiFab,AMREX_SPACEDIM>& facecent,
                                const Array<MultiFab,AMREX_SPACEDIM>& edgecent,
                                const MultiFab& levelset, const Geometry& geom,
                                const IntVect& ngrow, bool extend_domain_face,
                                int max_grid_size) const
{

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\n\t Writing checkpoint " << m_restart_file << std::endl;
    }

    const int nlevels = 1;
    PreBuildDirectorHierarchy(m_restart_file, level_prefix, nlevels, true);

    writeHeader(cut_grids, covered_grids, geom, ngrow, extend_domain_face, max_grid_size);

    writeToFile(volfrac, m_volfrac_name);
    writeToFile(centroid, m_centroid_name);
    writeToFile(bndryarea, m_bndryarea_name);
    writeToFile(bndrycent, m_bndrycent_name);
    writeToFile(bndrynorm, m_bndrynorm_name);
    writeToFile(levelset, m_levelset_name);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        writeToFile(areafrac[idim], m_areafrac_name[idim]);
        writeToFile(facecent[idim], m_facecent_name[idim]);
        writeToFile(edgecent[idim], m_edgecent_name[idim]);
    }
}

}
