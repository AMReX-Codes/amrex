//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.cpp,v 1.3 1998-11-28 21:58:43 lijewski Exp $
//
#ifdef BL_USE_NEW_HFILES
#include <cstring>
#else
#include <string.h>
#endif

#include "AmrLevel.H"
#include "ParmParse.H"
#include "StationData.H"

static int id_count = 0;

StationData::StationData () {}

StationData::~StationData () {}

inline
void
AddCoord (const Array<Real>& coord,
          StationRec&        rec)
{
    for (int n = 0; n < BL_SPACEDIM; n++)
        rec.pos[n] = coord[n];

    rec.id = id_count++;
}

void
StationData::init ()
{
    //
    // ParmParse variables:
    //
    //   StationData.vars     -- Array of name of StateData components
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //   StationData.rootname -- root name of output files.
    //
    ParmParse pp("StationData");

    if (pp.contains("vars"))
    {
        m_vars.resize(pp.countval("vars"));
        m_typ.resize(m_vars.length());
        m_ncomp.resize(m_vars.length());

        for (int i = 0; i < m_vars.length(); i++)
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
        Array<Real> data(BL_SPACEDIM);

        m_stn.resize(pp.countval("coord"));

        for (int k = 0; k < m_stn.length(); k++)
        {
            pp.getktharray("coord", k, data, 0, BL_SPACEDIM);

            AddCoord(data, m_stn[k]);
        }
    }

    m_name = "stations/stn";

    pp.query("rootname",m_name);

    if (m_name[m_name.length()-1] == '/')
        BoxLib::Error("StationData::init(): rootname must be valid filename");
    //
    // Make any directories assumed by m_name.
    //
    if (m_vars.length() > 0)
    {
        if (char* slash = strrchr(m_name.c_str(), '/'))
        {
            int idx = slash - m_name.c_str();

            assert(idx > 0);
            assert(m_name[idx] == '/');
            //
            // Some aString hanky-panky.
            //
            m_name[idx] = 0;
            aString dir = m_name.c_str();
            m_name[idx] = '/';
            //
            // Only the I/O processor makes the directory if it doesn't exist.
            //
            if (ParallelDescriptor::IOProcessor())
                if (!Utility::UtilCreateDirectory(dir, 0755))
                    Utility::CreateDirectoryFailed(dir);
        }
    }

    listStations();
}

void
StationData::openFile (int timestep)
{
   char buf[80];

   sprintf(buf, "_CPU_%04d_%06d", ParallelDescriptor::MyProc(), timestep);

   aString filename = m_name;

   filename += buf;

   m_ofile.open(filename.c_str(), ios::out|ios::app);

   m_ofile.precision(15);
}

void
StationData::listStations () const
{
    if (m_stn.length() > 0)
    {
        ofstream os("Station.List", ios::out);

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
StationData::report (Real time,
                     int  level) const
{
    if (m_stn.length() <= 0)
        return;

    static Array<Real> data;

    const int N      = m_vars.length();
    const int MyProc = ParallelDescriptor::MyProc();

    data.resize(N);

    const DescriptorList& desc_lst = AmrLevel::get_desc_lst();

    for (int i = 0; i < m_stn.length(); i++)
    {
        if (m_stn[i].own && m_stn[i].level == level)
        {
            //
            // Fill in the `data' vector.
            //
            for (int j = 0; j < N; j++)
            {
                const MultiFab& mf = desc_lst[m_typ[j]].newData();

                assert(mf.DistributionMap()[m_stn[i].grd] == MyProc);

                assert(mf.nComp() > m_ncomp[j]);

                data[j] = mf[m_stn[i].grd](m_stn[i].ix,m_ncomp[j]);
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
}

void
StationData::findGrid (const PArray<AmrLevel>& levels,
                       const Array<Geometry>&  geoms)
{
    assert(geom.length() == levels.length());

    if (m_stn.length() <= 0)
        return;
    //
    // Flag all stations as not having a home.
    //
    for (int i = 0; i < m_stn.length(); i++)
    {
        m_stn[i].level = -1;
    }
    //
    // Find level and grid owning the data.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    for (int level = levels.length(); level >= 0; level--)
    {
        const BoxArray& ba          = levels[level].boxArray();
        const Array<RealBox>& boxes = levels[level].gridLocations();

        for (int i = 0; i < m_stn.length(); i++)
        {
            if (m_stn[i].level < 0)
            {
                const Real* pos = &m_stn[i].pos[0];

                for (int j = 0; j < boxes.length(); j++)
                {
                    if (boxes[j].contains(pos))
                    {
                        m_stn[i].level = level;
                        m_stn[i].grd   = j;
                        m_stn[i].ix    = geoms[level].CellIndex(pos);
                        //
                        // Does this CPU own the data?
                        //
                        MultiFab mf(ba,1,0,Fab_noallocate);

                        m_stn[i].own = (mf.DistributionMap()[j] == MyProc);

                        break;
                    }
                }
            }
        }
    }
}
