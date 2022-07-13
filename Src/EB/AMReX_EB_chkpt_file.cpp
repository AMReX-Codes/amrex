#include <AMReX_EB_chkpt_file.H>

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)

#include <AMReX_MultiFabUtil.H>


namespace {

const std::string level_prefix = "Level_";

void gotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

}

namespace amrex { namespace EB2 {

void ChkptFile::writeHeader (const BoxArray& ba) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        std::string HeaderFileName(m_restart_file + "/Header");
        VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;

        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                std::ofstream::trunc |
                std::ofstream::binary);

        if ( ! HeaderFile.good() )
            amrex::FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << "Checkpoint version: 1\n";

        const int nlevels = 1;
        HeaderFile << nlevels << "\n";

        Geometry geometry;

        // Geometry
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            HeaderFile << geometry.ProbLo(i) << ' ';
        HeaderFile << '\n';

        for (int i = 0; i < AMREX_SPACEDIM; ++i)
            HeaderFile << geometry.ProbHi(i) << ' ';
        HeaderFile << '\n';

        // BoxArray
        for (int lev = 0; lev < nlevels; ++lev)
        {
            ba.writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }
}

void ChkptFile::writeToFile (const MultiFab& mf, const std::string& mf_name) const
{
    amrex::VisMF::Write(mf, amrex::MultiFabFileFullPrefix(0, m_restart_file,
                level_prefix, mf_name));
}


ChkptFile::ChkptFile (const std::string &fname)
    : m_restart_file(fname)
{}

void ChkptFile::fill_from_chkpt_file(BoxArray& grids, DistributionMapping& dmap,
        MultiFab& volfrac, MultiFab& centroid, MultiFab& bndryarea, MultiFab& bndrycent,
        MultiFab& bndrynorm, Array<MultiFab,AMREX_SPACEDIM>& areafrac,
        Array<MultiFab,AMREX_SPACEDIM>& facecent,
        Array<MultiFab,AMREX_SPACEDIM>& edgecent,
        MultiFab& levelset, int ng) const
{
    Real prob_lo[AMREX_SPACEDIM];
    Real prob_hi[AMREX_SPACEDIM];

    std::string File(m_restart_file + "/Header");

    Print() << "file=" << File << std::endl;

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    std::getline(is, line);

    int  nlevs;
    int  int_tmp;
    Real real_tmp;

    is >> nlevs;
    gotoNextLine(is);
    AMREX_ASSERT(nlevs == 1);

    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            prob_lo[i++] = std::stod(word);
        }
    }

    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            prob_hi[i++] = std::stod(word);
        }
    }

    BoxArray orig_ba;
    orig_ba.readFrom(is);
    gotoNextLine(is);

    Box orig_domain(orig_ba.minimalBox());

    BoxList bl;
    for (int nb = 0; nb < orig_ba.size(); nb++) {
        Box b(orig_ba[nb]);
        bl.push_back(b);
    }

    grids.define(bl);
    dmap.define(grids, ParallelDescriptor::NProcs());

    // volfrac
    {
        Print() << "  Loading volfrac" << std::endl;

        volfrac.define(grids, dmap, 1, ng);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_volfrac_name);
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        volfrac.ParallelCopy(mf, 0, 0, 1, 0, ng);
    }

    // centroid
    {
        Print() << "  Loading centroid" << std::endl;

        centroid.define(grids, dmap, AMREX_SPACEDIM, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_centroid_name);
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        centroid.ParallelCopy(mf, 0, 0, AMREX_SPACEDIM, 0, ng);
    }

    // bndryarea
    {
        Print() << "  Loading bndryarea" << std::endl;

        bndryarea.define(grids, dmap, 1, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_bndryarea_name);
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        bndryarea.ParallelCopy(mf, 0, 0, 1, 0, ng);
    }

    // bndrycent
    {
        Print() << "  Loading bndrycent" << std::endl;

        bndrycent.define(grids, dmap, AMREX_SPACEDIM, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_bndrycent_name);
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        bndrycent.ParallelCopy(mf, 0, 0, AMREX_SPACEDIM, 0, ng);
    }

    // bndrynorm
    {
        Print() << "  Loading bndrynorm" << std::endl;

        bndrynorm.define(grids, dmap, AMREX_SPACEDIM, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_bndrynorm_name);
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        bndrynorm.ParallelCopy(mf, 0, 0, AMREX_SPACEDIM, 0, ng);
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // areafrac
        {
            Print() << "  Loading areafrac " << idim << std::endl;

            areafrac[idim].define(convert(grids, IntVect::TheDimensionVector(idim)), dmap, 1, ng);

            auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_areafrac_name[idim]);
            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);
            areafrac[idim].ParallelCopy(mf, 0, 0, 1, 0, ng);
        }

        // facecent
        {
            Print() << "  Loading facecent " << idim << std::endl;

            facecent[idim].define(convert(grids, IntVect::TheDimensionVector(idim)), dmap, AMREX_SPACEDIM-1, ng);

            auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_facecent_name[idim]);
            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);
            facecent[idim].ParallelCopy(mf, 0, 0, AMREX_SPACEDIM-1, 0, ng);
        }

        // edgecent
        {
            Print() << "  Loading edgecent " << idim << std::endl;

            IntVect edge_type{1}; edge_type[idim] = 0;
            edgecent[idim].define(convert(grids, edge_type), dmap, 1, ng);

            auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_edgecent_name[idim]);
            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);
            edgecent[idim].ParallelCopy(mf, 0, 0, 1, 0, ng);
        }
    }

    // levelset
    {
        Print() << "  Loading levelset" << std::endl;

        levelset.define(convert(grids,IntVect::TheNodeVector()), dmap, 1, ng);
        levelset.setVal(-1.);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, level_prefix, m_levelset_name);
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        levelset.ParallelCopy(mf, 0, 0, 1, 0, ng);
    }
}

void ChkptFile::write_to_chkpt_file (const BoxArray& ba,
        const MultiFab& volfrac,
        const MultiFab& centroid, const MultiFab& bndryarea,
        const MultiFab& bndrycent, const MultiFab& bndrynorm,
        const Array<MultiFab,AMREX_SPACEDIM>& areafrac,
        const Array<MultiFab,AMREX_SPACEDIM>& facecent,
        const Array<MultiFab,AMREX_SPACEDIM>& edgecent,
        const MultiFab& levelset) const
{

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "\n\t Writing checkpoint " << m_restart_file << std::endl;
    }

    const int nlevels = 1;
    amrex::PreBuildDirectorHierarchy(m_restart_file, level_prefix, nlevels, true);

    writeHeader(ba);

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

}}
