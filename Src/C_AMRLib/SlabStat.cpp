//BL_COPYRIGHT_NOTICE

//
// $Id: SlabStat.cpp,v 1.2 2000-03-01 15:43:35 lijewski Exp $
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
    m_interval(0)
{
    m_mf.setVal(0);
    //
    // Build m_tmp_mf without using ghost cells.
    // This'll facilitate using it in Parallel MF routines.
    //
    BoxArray ba(boxes.length());
        
    for (int j = 0; j < boxes.length(); ++j)
        ba.set(j,::grow(boxes[j],ngrow));

    m_tmp_mf.define(ba,m_vars.length(),0,Fab_allocate);
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
            const int L = li()->level();
            MultiFab& M = li()->tmp_mf();
            const int N = M.nGrow();

            li()->m_interval += dt;
            //
            // Build MultiFab with no ghost cells paralleling the AmrLevel.
            //
            BoxArray ba(amrlevel.boxArray().length());

            for (int j = 0; j < ba.length(); ++j)
                ba.set(j,::grow(amrlevel.boxArray()[j],N));

            MultiFab tmpMF(ba,1,0);

            for (int i = 0; i < M.nComp(); i++)
            {
                MultiFab* mf = amrlevel.derive(li()->vars()[i],time+dt,N);
                //
                // Fill tmpMF from *mf.  They have the same distribution.
                //
                BL_ASSERT(ba.length() == mf->boxArray().length());

                for (MultiFabIterator mfi(*mf); mfi.isValid(); ++mfi)
                    tmpMF[mfi.index()].copy((*mf)[mfi.index()],0,0,1);

                delete mf;
                //
                // Parallel MultiFab -> MultiFab copy.
                // This'll do what we want since neither has ghost cells.
                //
                M.copy(tmpMF,0,i,1);
            }
            //
            // Finally compute the statistics.
            //
            for (MultiFabIterator dmfi(li()->mf()); dmfi.isValid(); ++dmfi)
            {
                DependentMultiFabIterator smfi(dmfi,M);

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
        RunStats::addBytes(VisMF::Write(li()->mf(),path+li()->name()));
}
