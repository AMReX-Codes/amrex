//BL_COPYRIGHT_NOTICE

//
// $Id: SlabStat.cpp,v 1.3 2000-03-01 20:19:03 lijewski Exp $
//

#include <AmrLevel.H>
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
SlabStatList::checkPoint (const aString& ckdir)
{
    if (m_list.isEmpty()) return;
    //
    // We put SlabStats in subdirectory of `ckdir'.
    //
    const aString statdir = ckdir + "/stats";
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

        HeaderFile << m_list.length() << '\n';

        for (ListIterator<SlabStatRec*> li(m_list); li; ++li)
            HeaderFile << li()->name() << '\n';

        HeaderFile << m_list.firstElement()->interval() << '\n';
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
