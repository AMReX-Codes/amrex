//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.cpp,v 1.8 1998-11-30 21:13:38 lijewski Exp $
//
#ifdef BL_USE_NEW_HFILES
#include <cstdio>
#include <cstring>
#else
#include <stdio.h>
#include <string.h>
#endif

#include <AmrLevel.H>
#include <ParmParse.H>
#include <StationData.H>

StationData::StationData () {}

StationData::~StationData ()
{
    m_ofile.close();
}

void
StationData::init ()
{
    //
    // ParmParse variables:
    //
    //   StationData.vars     -- Names of StateData components to output
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //
    ParmParse pp("StationData");

    if (pp.contains("vars"))
    {
        const int N = pp.countval("vars");

        m_vars.resize(N);
        m_typ.resize(N);
        m_ncomp.resize(N);

        for (int i = 0; i < N; i++)
        {
            pp.get("vars", m_vars[i], i);

            if (!AmrLevel::isStateVariable(m_vars[i],m_typ[i],m_ncomp[i]))
            {
                cout << "StationData::init(): `"
                     << m_vars[i]
                     << "' is not a state variable\n";
                BoxLib::Abort();
            }
        }
    }

    if (m_vars.length() > 0 && pp.contains("coord"))
    {
        static int identifier;

        Array<Real> data(BL_SPACEDIM);

        const int N  = pp.countname("coord");

        m_stn.resize(N);

        for (int k = 0; k < N; k++)
        {
            pp.getktharr("coord", k, data, 0, BL_SPACEDIM);

            D_TERM(m_stn[k].pos[0] = data[0];,
                   m_stn[k].pos[1] = data[1];,
                   m_stn[k].pos[2] = data[2];);

            assert(Geometry::ProbDomain().contains(data.dataPtr()));

            m_stn[k].id = identifier++;
        }
    }

    if (m_vars.length() > 0)
    {
        //
        // Make "Stations/" directory.
        //
        // Only the I/O processor makes the directory if it doesn't exist.
        //
        static const aString Dir("Stations");

        if (ParallelDescriptor::IOProcessor())
            if (!Utility::UtilCreateDirectory(Dir, 0755))
                Utility::CreateDirectoryFailed(Dir);
        //
        // Everyone must wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Open the data file.
        //
        char buf[80];

        sprintf(buf,
                "Stations/stn_CPU_%04d",
                ParallelDescriptor::MyProc());

        m_ofile.open(buf, ios::out|ios::app);

        m_ofile.precision(15);

        assert(!m_ofile.bad());
        //
        // Output the list of stations.
        //
        ofstream os("Stations/Station.List", ios::out);

        for (int i = 0; i < m_stn.length(); i++)
        {
            os << m_stn[i].id;

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                os << '\t' << m_stn[i].pos[k];
            }

            os << '\n';
        }
    }
}

void
StationData::report (Real            time,
                     int             level,
                     const AmrLevel& amrlevel)
{
    if (m_stn.length() <= 0)
        return;

    static Array<Real> data;

    const int N      = m_vars.length();
    const int MyProc = ParallelDescriptor::MyProc();

    data.resize(N);

    for (int i = 0; i < m_stn.length(); i++)
    {
        if (m_stn[i].own && m_stn[i].level == level)
        {
            //
            // Fill in the `data' vector.
            //
            for (int j = 0; j < N; j++)
            {
                const MultiFab& mf = amrlevel.get_new_data(m_typ[j]);

                assert(mf.nComp() > m_ncomp[j]);
                assert(mf.DistributionMap()[m_stn[i].grd] == MyProc);
                //
                // Find IntVect so we can index into FAB.
                // We want to use Geometry::CellIndex().
                // Must adjust the position to account for NodeCentered-ness.
                //
                IndexType ityp = amrlevel.get_desc_lst()[m_typ[j]].getType();

                Real pos[BL_SPACEDIM];

                D_TERM(pos[0] = m_stn[i].pos[0] + .5 * ityp[0];,
                       pos[1] = m_stn[i].pos[1] + .5 * ityp[1];,
                       pos[2] = m_stn[i].pos[2] + .5 * ityp[2];);

                IntVect idx = amrlevel.Geom().CellIndex(&pos[0]);

                data[j] = mf[m_stn[i].grd](idx,m_ncomp[j]);
            }
            //
            // Write data to output stream.
            //
            m_ofile << m_stn[i].id << ' ' << time << ' ';

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                m_ofile << m_stn[i].pos[k] << ' ';
            }
            for (int k = 0; k < N; k++)
            {
                m_ofile << data[k] << ' ';
            }
            m_ofile << '\n';
        }
    }

    m_ofile.flush();

    assert(!m_ofile.bad());
}

void
StationData::findGrid (const PArray<AmrLevel>& levels,
                       const Array<Geometry>&  geoms)
{
    assert(geoms.length() == levels.length());

    if (m_stn.length() <= 0)
        return;
    //
    // Flag all stations as not having a home.
    //
    for (int i = 0; i < m_stn.length(); i++)
        m_stn[i].level = -1;
    //
    // Find level and grid owning the data.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    for (int level = levels.length()-1; level >= 0; level--)
    {
        const Array<RealBox>& boxes = levels[level].gridLocations();

        MultiFab mf(levels[level].boxArray(),1,0,Fab_noallocate);

        assert(mf.boxArray().length() == boxes.length());

        for (int i = 0; i < m_stn.length(); i++)
        {
            if (m_stn[i].level < 0)
            {
                for (int j = 0; j < boxes.length(); j++)
                {
                    if (boxes[j].contains(&m_stn[i].pos[0]))
                    {
                        m_stn[i].grd   = j;
                        m_stn[i].own   = (mf.DistributionMap()[j] == MyProc);
                        m_stn[i].level = level;

                        break;
                    }
                }
            }
        }
    }
}
