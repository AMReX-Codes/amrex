#include <AMReX_EB_chkpt_file.H>

#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>    // amrex::VisMF::Write(MultiFab)
#include <AMReX_VectorIO.H> // amrex::[read,write]IntData(array_of_ints)

#include <AMReX_MultiFabUtil.H>


namespace {

const int ng_to_copy = 0;

void gotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
}

namespace amrex { namespace EB2 {

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

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, "Level_", "volfrac");
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        volfrac.ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng);
    }

    // centroid
    {
        Print() << "  Loading centroid" << std::endl;

        centroid.define(grids, dmap, AMREX_SPACEDIM, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", "centroid");
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        centroid.ParallelCopy(mf, 0, 0, AMREX_SPACEDIM, ng_to_copy, ng);
    }

    // bndryarea
    {
        Print() << "  Loading bndryarea" << std::endl;

        bndryarea.define(grids, dmap, 1, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", "bndryarea");
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        bndryarea.ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng);
    }

    // bndrycent
    {
        Print() << "  Loading bndrycent" << std::endl;

        bndrycent.define(grids, dmap, AMREX_SPACEDIM, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", "bndrycent");
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        bndrycent.ParallelCopy(mf, 0, 0, AMREX_SPACEDIM, ng_to_copy, ng);
    }

    // bndrynorm
    {
        Print() << "  Loading bndrynorm" << std::endl;

        bndrynorm.define(grids, dmap, AMREX_SPACEDIM, ng);

        auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", "bndrynorm");
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        bndrynorm.ParallelCopy(mf, 0, 0, AMREX_SPACEDIM, ng_to_copy, ng);
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // areafrac
        {
            Print() << "  Loading areafrac " << idim << std::endl;

            areafrac[idim].define(convert(grids, IntVect::TheDimensionVector(idim)), dmap, 1, ng);

            auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", Concatenate("areafrac", idim, 1));
            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);
            areafrac[idim].ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng);
        }

        // facecent
        {
            Print() << "  Loading facecent " << idim << std::endl;

            facecent[idim].define(convert(grids, IntVect::TheDimensionVector(idim)), dmap, AMREX_SPACEDIM-1, ng);

            auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", Concatenate("facecent", idim, 1));
            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);
            facecent[idim].ParallelCopy(mf, 0, 0, AMREX_SPACEDIM-1, ng_to_copy, ng);
        }

        // edgecent
        {
            Print() << "  Loading edgecent " << idim << std::endl;

            IntVect edge_type{1}; edge_type[idim] = 0;
            edgecent[idim].define(convert(grids, edge_type), dmap, 1, ng);

            auto prefix = amrex::MultiFabFileFullPrefix(0, m_restart_file, "Level_", Concatenate("edgecent", idim, 1));
            MultiFab mf(The_Pinned_Arena());
            VisMF::Read(mf, prefix);
            edgecent[idim].ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng);
        }
    }

    // levelset
    {
        Print() << "  Loading levelset" << std::endl;

        levelset.define(convert(grids,IntVect::TheNodeVector()), dmap, 1, ng);
        levelset.setVal(-1.);

        auto prefix = MultiFabFileFullPrefix(0, m_restart_file, "Level_", "levelset");
        MultiFab mf(The_Pinned_Arena());
        VisMF::Read(mf, prefix);
        levelset.ParallelCopy(mf, 0, 0, 1, ng_to_copy, ng);
    }
}

}}
