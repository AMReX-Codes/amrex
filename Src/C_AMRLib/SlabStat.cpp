//
// $Id: SlabStat.cpp,v 1.12 2001-08-01 21:50:46 lijewski Exp $
//

#include <AmrLevel.H>
#include <ParmParse.H>
#include <SlabStat.H>

SlabStatRec::SlabStatRec (const std::string&  name,
                          int                 ncomp,
                          Array<std::string>& vars,
                          int                 ngrow,
                          int                 level,
                          const BoxArray&     boxes,
                          SlabStatFunc        func)
    :
    m_name(name),
    m_ncomp(ncomp),
    m_vars(vars),
    m_ngrow(ngrow),
    m_level(level),
    m_boxes(boxes),
    m_func(func),
    m_mf(m_boxes,m_ncomp,0),
    m_tmp_mf(m_boxes,m_vars.size(),m_ngrow),
    m_interval(0)
{
    m_mf.setVal(0);
}

const std::string&
SlabStatRec::name () const
{
    return m_name;
}

int
SlabStatRec::nComp () const
{
    return m_ncomp;
}

const Array<std::string>&
SlabStatRec::vars () const
{
    return m_vars;
}

long
SlabStatRec::nVariables () const
{
    return m_vars.size();
}

int
SlabStatRec::nGrow () const
{
    return m_ngrow;
}

int
SlabStatRec::level () const
{
    return m_level;
}

const BoxArray&
SlabStatRec::boxes () const
{
    return m_boxes;
}

SlabStatFunc
SlabStatRec::func () const
{
    return m_func;
}

MultiFab&
SlabStatRec::mf ()
{
    return m_mf;
}

MultiFab&
SlabStatRec::tmp_mf ()
{
    return m_tmp_mf;
}

Real
SlabStatRec::interval () const
{
    return m_interval;
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
Boxes (const std::string& file,
       const std::string& name,
       BoxArray&          boxes,
       int&               boxesLevel)
{
    const std::string TheDflt = "default:";
    const std::string TheName = name + ":";

    std::ifstream is(file.c_str(),std::ios::in);

    if (!is.good())
        BoxLib::FileOpenFailed(file);

#define STRIP while( is.get() != '\n' )

    BoxArray ba_dflt;
    BoxArray ba_name;
    int bxLvl_dflt;
    int bxLvl_name;

    std::string line;

    while (std::getline(is,line))
    {
        if (line.empty() || line[0] == '#') continue;

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

    if (ba_dflt.size() == 0 && ba_name.size() == 0)
        BoxLib::Abort("slabstats.boxes doesn't have appropriate structure");

    boxes =  ba_name.size() ? ba_name : ba_dflt;
    boxesLevel = ba_name.size() ? bxLvl_name : bxLvl_dflt;
#undef STRIP
}

SlabStatRec::SlabStatRec (const std::string&  name,
                          int                 ncomp,
                          Array<std::string>& vars,
                          int                 ngrow,
                          SlabStatFunc        func)
    :
    m_name(name),
    m_ncomp(ncomp),
    m_vars(vars),
    m_ngrow(ngrow),
    m_func(func),
    m_interval(0)
{
    ParmParse pp("slabstat");

    std::string file;

    if (!pp.query("boxes", file))
        BoxLib::Abort("SlabStatRec: slabstat.boxes isn't defined");

    Boxes(file, m_name, m_boxes, m_level);

    m_mf.define(m_boxes, m_ncomp, 0, Fab_allocate);

    m_tmp_mf.define(m_boxes, m_vars.size(), m_ngrow, Fab_allocate);

    m_mf.setVal(0);
}

SlabStatRec::~SlabStatRec () {}

SlabStatList::SlabStatList () {}

SlabStatList::~SlabStatList ()
{
    for (std::list<SlabStatRec*>::iterator li = m_list.begin();
         li != m_list.end();
         ++li)
    {
        delete *li;
    }
}

void
SlabStatList::add (const std::string&  name,
                   int                 ncomp,
                   Array<std::string>& vars,
                   int                 ngrow,
                   int                 level,
                   const BoxArray&     boxes,
                   SlabStatFunc        func)
{
    m_list.push_back(new SlabStatRec(name,ncomp,vars,ngrow,level,boxes,func));
}

void
SlabStatList::add (const std::string&  name,
                   int                 ncomp,
                   Array<std::string>& vars,
                   int                 ngrow,
                   SlabStatFunc        func)
{
    m_list.push_back(new SlabStatRec(name,ncomp,vars,ngrow,func));
}

std::list<SlabStatRec*>&
SlabStatList::list ()
{
    return m_list;
}

void
SlabStatList::update (AmrLevel& amrlevel,
                      Real      time,
                      Real      dt)
{
    for (std::list<SlabStatRec*>::iterator li = m_list.begin();
         li != m_list.end();
         ++li)
    {
        if ((*li)->level() == amrlevel.Level())
        {
            (*li)->m_interval += dt;

            for (int i = 0; i < (*li)->tmp_mf().nComp(); i++)
            {
                amrlevel.derive((*li)->vars()[i],time+dt,(*li)->tmp_mf(),i);
            }

            for (MFIter dmfi((*li)->mf()); dmfi.isValid(); ++dmfi)
            {
                FArrayBox&       dfab = (*li)->mf()[dmfi];
                const FArrayBox& sfab = (*li)->tmp_mf()[dmfi];

                const int nsrc = sfab.nComp();
                const int ndst = dfab.nComp();

                (*li)->func()(sfab.dataPtr(),
                              ARLIM(sfab.box().loVect()),
                              ARLIM(sfab.box().hiVect()),
                              &nsrc,
                              dfab.dataPtr(),
                              ARLIM(dfab.box().loVect()),
                              ARLIM(dfab.box().hiVect()),
                              &ndst,
                              &dt,
                              amrlevel.Geom().CellSize());
            }
        }
    }
}

void
SlabStatList::checkPoint (PArray<AmrLevel>& amrLevels,
                          int               level0_step)
{
    if (m_list.empty()) return;
    //
    // We put SlabStats in a subdirectory of the directory, 'slabstats'.
    //
    std::string statdir = "slabstats";
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(statdir, 0755))
            BoxLib::CreateDirectoryFailed(statdir);

    statdir = statdir + "/stats";
    statdir = BoxLib::Concatenate(statdir, level0_step);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(statdir, 0755))
            BoxLib::CreateDirectoryFailed(statdir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Only the IOProcessor() writes to the header file.
        //
        const std::string HeaderFileName = statdir + "/Header";

        std::ofstream HeaderFile;

        HeaderFile.open(HeaderFileName.c_str(), std::ios::out|std::ios::trunc);

        if (!HeaderFile.good())
            BoxLib::FileOpenFailed(HeaderFileName);

        int prec = HeaderFile.precision(30);

        HeaderFile << m_list.size() << '\n';

        for (std::list<SlabStatRec*>::const_iterator li = m_list.begin();
             li != m_list.end();
             ++li)
        {
            HeaderFile << (*li)->name() << " " << (*li)->level() << '\n';
        }

        HeaderFile << (*m_list.begin())->interval() << '\n';

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            HeaderFile << Geometry::ProbLo(dir) << " ";
        HeaderFile << '\n';

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            HeaderFile << Geometry::ProbHi(dir) << " ";
        HeaderFile << '\n';

        for (int level = 0; level < amrLevels.size(); level++)
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
    const std::string path = statdir + "/";

    for (std::list<SlabStatRec*>::iterator li = m_list.begin();
         li != m_list.end();
         ++li)
    {
        VisMF::Write((*li)->mf(),path+(*li)->name());

        (*li)->m_interval = 0;

        (*li)->mf().setVal(0);
    }
}
