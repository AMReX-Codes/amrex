//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.cpp,v 1.2 1998-11-25 23:58:42 lijewski Exp $
//
#ifdef BL_USE_NEW_HFILES
#include <cstring>
#else
#include <string.h>
#endif

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

    rec.id    = id_count++;
    rec.level = -1;
    rec.grd   = 0;
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

        for (int i = 0; i < m_vars.length(); i++)
        {
            pp.get("vars", m_vars[i], i);
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
StationData::report (Real time) const
{
    static Real* data      = 0;
    static int   data_size = 0;

    if (nvar > data_size)
    {
        data      = new Real[nvar];
        data_size = nvar;
    }
    for (int i = 0; i < m_stn.length(); i++)
    {
        m_stn[i].grd->getData(m_stn[i].ix,data);
        //
        // Write data to output stream.
        //
        m_ofile << m_stn[i].id << ' ' << time << ' ';

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            m_ofile << m_stn[i].pos[k] << ' ';
        }
        for (int k = 0; k < nvar; k++)
        {
            m_ofile << data[k] << ' ';
        }
        m_ofile << '\n';
    }
}

void
StationData::findGrid (const GridList* gl,
                       int             hi_lev,
                       Real            dx_lev[MAX_LEV][BL_SPACEDIM],
                       Real            prob_lo[BL_SPACEDIM])
{
    if (m_stn.length() <= 0)
        return;
    //
    // Flag all stations as not having a home.
    //
    for (int i = 0; i < m_stn.length(); i++)
    {
        m_stn[i].level = -1;
        m_stn[i].grd   = 0;
    }
    //
    // Walk list of grids from highest level down.
    //
    for (int level = hi_lev; level >= 0; level--)
    {
        for (Grid* g = gl[level].first(); g != 0; g = gl[level].next(g))
        {
            //
            // Find all station points in this grid.
            //
            const Box& gbx = g->box();

            for (int i = 0; i < m_stn.length(); i++)
            {
                if (m_stn[i].level < 0)
                {
                    for (int n = 0; n < BL_SPACEDIM; n++)
                    {
                        int ipos = int((m_stn[i].pos[n]-prob_lo[n])/dx_lev[level][n]);
                        m_stn[i].ix.setVal(n,ipos);
                    }
                    if (gbx.contains(m_stn[i].ix))
                    {
                        m_stn[i].level = level;
                        m_stn[i].grd   = g;
                    }
                }
            }
        }
    }
}
