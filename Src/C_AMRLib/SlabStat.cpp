//BL_COPYRIGHT_NOTICE

//
// $Id: SlabStat.cpp,v 1.9 2000-07-11 21:03:56 sstanley Exp $
//

#include <AmrLevel.H>
#include <ParmParse.H>
#include <SlabStat.H>

SlabStatRec::SlabStatRec (const aString&  name,
                          int             ncomp,
                          Array<aString>& vars,
                          int             ngrow,
                          int             level,
                          const BoxArray& boxes,
                          SlabStatFunc    func)
    :
    m_name(name),
    m_ncomp(ncomp),
    m_vars(vars),
    m_ngrow(ngrow),
    m_level(level),
    m_boxes(boxes),
    m_func(func),
    m_mf(m_boxes,m_ncomp,0),
    m_tmp_mf(m_boxes,m_vars.length(),m_ngrow),
    m_interval(0)
{
    m_mf.setVal(0);
}

//
// This is a helper function for the following constructor.
//
// It opens "file" and looks for the Boxes on which the statistic
// "name" is to be calculated.
//
// The structure of the file:
//
//   #
//   # Lines beginning with a pound sign are comments.
//   #
//   default:
//     N
//     level
//     b1
//     b2
//     ...
//     bN
//   name:
//     M
//     level
//     b1
//     b2
//     ...
//     bM
//

static
void
Boxes (const aString& file,
       const aString& name,
       BoxArray&      boxes,
       int&           boxesLevel)
{
    const aString TheDflt = "default:";
    const aString TheName = name + ":";

    ifstream is(file.c_str(),ios::in);

    if (!is.good())
        Utility::FileOpenFailed(file);

#define STRIP while( is.get() != '\n' )

    BoxArray ba_dflt;
    BoxArray ba_name;
    int bxLvl_dflt;
    int bxLvl_name;

    aString line;

    while (line.getline(is))
    {
        if (line.isNull() || line[0] == '#') continue;

        if (line == TheDflt || line == TheName)
        {
            bool dflt = (line == TheDflt) ? true : false;

            int N; int lvl; Box bx; BoxList bl;

            is >> N; STRIP;
            is >> lvl; STRIP;

            BL_ASSERT(N > 0);

            for (int i = 0; i < N; i++)
            {
                is >> bx; STRIP; bl.append(bx);
            }

            (dflt ? ba_dflt : ba_name) = BoxArray(bl);
            (dflt ? bxLvl_dflt : bxLvl_name) = lvl;
        }
    }

    is.close();

    if (!ba_dflt.ready() && !ba_name.ready())
        BoxLib::Abort("slabstats.boxes doesn't have appropriate structure");

    boxes =  ba_name.ready() ? ba_name : ba_dflt;
    boxesLevel = ba_name.ready() ? bxLvl_name : bxLvl_dflt;
#undef STRIP
}

SlabStatRec::SlabStatRec (const aString&  name,
                          int             ncomp,
                          Array<aString>& vars,
                          int             ngrow,
                          SlabStatFunc    func)
    :
    m_name(name),
    m_ncomp(ncomp),
    m_vars(vars),
    m_ngrow(ngrow),
    m_func(func),
    m_interval(0)
{
    ParmParse pp("slabstat");

    aString file;

    if (!pp.query("boxes", file))
        BoxLib::Abort("SlabStatRec: slabstat.boxes isn't defined");

    Boxes(file, m_name, m_boxes, m_level);

    m_mf.define(m_boxes, m_ncomp, 0, Fab_allocate);

    m_tmp_mf.define(m_boxes, m_vars.length(), m_ngrow, Fab_allocate);

    m_mf.setVal(0);
}

SlabStatRec::~SlabStatRec () {}

SlabStatList::~SlabStatList ()
{
    for (ListIterator<SlabStatRec*> li(m_list); li; ++li)
    {
        delete li();
    }
}

void
SlabStatList::update (AmrLevel& amrlevel,
                      Real      time,
                      Real      dt)
{
    for (ListIterator<SlabStatRec*> li(m_list); li; ++li)
    {
        if (li()->level() == amrlevel.Level())
        {
            li()->m_interval += dt;

            for (int i = 0; i < li()->tmp_mf().nComp(); i++)
            {
                amrlevel.derive(li()->vars()[i],time+dt,li()->tmp_mf(),i);
            }

            for (MultiFabIterator dmfi(li()->mf()); dmfi.isValid(); ++dmfi)
            {
                DependentMultiFabIterator smfi(dmfi,li()->tmp_mf());

                const int nsrc = smfi().nComp();
                const int ndst = dmfi().nComp();

                li()->func()(smfi().dataPtr(),
                             ARLIM(smfi().box().loVect()),
                             ARLIM(smfi().box().hiVect()),
                             &nsrc,
                             dmfi().dataPtr(),
                             ARLIM(dmfi().box().loVect()),
                             ARLIM(dmfi().box().hiVect()),
                             &ndst,
                             &dt,
                             amrlevel.Geom().CellSize());
            }
        }
    }
}

void
SlabStatList::checkPoint (PArray<AmrLevel>& amrLevels,
                          const int level0_step)
{
    if (m_list.isEmpty()) return;
    //
    // We put SlabStats in a subdirectory of the directory, 'slabstats'.
    //
    aString statdir = "slabstats";
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::UtilCreateDirectory(statdir, 0755))
            Utility::CreateDirectoryFailed(statdir);

    statdir = statdir + "/stats";
    statdir = Utility::Concatenate(statdir, level0_step);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!Utility::UtilCreateDirectory(statdir, 0755))
            Utility::CreateDirectoryFailed(statdir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        const aString HeaderFileName = statdir + "/Header";

        ofstream HeaderFile;

        HeaderFile.open(HeaderFileName.c_str(), ios::out|ios::trunc);

        if (!HeaderFile.good())
            Utility::FileOpenFailed(HeaderFileName);

        int prec = HeaderFile.precision(30);

        HeaderFile << m_list.length() << '\n';

        for (ListIterator<SlabStatRec*> li(m_list); li; ++li)
            HeaderFile << li()->name() << " " << li()->level() << '\n';

        HeaderFile << m_list.firstElement()->interval() << '\n';

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            HeaderFile << Geometry::ProbLo(dir) << " ";
        HeaderFile << '\n';

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            HeaderFile << Geometry::ProbHi(dir) << " ";
        HeaderFile << '\n';

        for (int level = 0; level < amrLevels.length(); level++)
        {
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
                HeaderFile << amrLevels[level].Geom().CellSize(dir) << " ";

            HeaderFile << '\n';
        }

        HeaderFile.precision(prec);

        if (!HeaderFile.good())
            BoxLib::Error("SlabStat::checkPoint() failed");
    }
    //
    // Write out the SlabStat MultiFabs.
    //
    const aString path = statdir + "/";

    for (ListIterator<SlabStatRec*> li(m_list); li; ++li)
    {
        RunStats::addBytes(VisMF::Write(li()->mf(),path+li()->name()));

        li()->m_interval = 0;

        li()->mf().setVal(0);
    }
}
