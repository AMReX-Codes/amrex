//BL_COPYRIGHT_NOTICE

//
// $Id: StationData.cpp,v 1.1 1998-11-25 19:14:16 lijewski Exp $
//

#include "StationData.H"

//
// Remove definition if you want station data output in native binary format.
//
#define STN_ASCII 1

static int id_count = 0;

StationData::StationData () {}

StationData::~StationData () {}

void
StationData::openFile (int suffix)
{
   char buf[80];

   sprintf(buf, ".%06d", suffix);

   aString filename = m_name;

   filename += buf;

   m_ofile.open(filename.c_str(),ios::out);

#ifdef   STN_ASCII
//   Following 2 statements commented by ralph ferguson  14Feb96
//   m_ofile << "ASCII" << '\n';
//   m_ofile << m_nvar << '\n';
#else
   m_ofile << "BINARY" << '\n';
   m_ofile.write( (char*) &m_nvar, sizeof(int) );
#endif
}

void
StationData::init (int         nvar,
                   const char* root_name)
{
   m_name = root_name
   m_nvar = nvar;
}

void
StationData::addStation (Real* coord)
{
    StationRec rec;

    for (int n = 0; n < BL_SPACEDIM; n++)
        rec.pos[n] = coord[n];

    rec.id    = id_count++;
    rec.level = -1;
    rec.grd   = 0;

    m_stn.push_back(rec);
}

//
// Ascii write to checkpoint file.
//
void
StationData::dumpStation (ofstream& os,
                          int       step)
{
    os << m_nvar       << '\n';
    os << m_stn.size() << '\n';

    if (m_stn.size() > 0)
    {
        os << m_name << '\n';

        for (int i = 0; i < m_stn.size(); i++)
        {
            os << m_stn[i].id;

            for (int k = 0; k < BL_SPACEDIM; k++)
                os << ' ' << m_stn[i].pos[k];

            os << '\n';
        }
        //
        // Close current Station stream and construct a new one.
        //
        m_ofile.close();
        openFile(step);
    }
}

void
StationData::listStations ()
{
   ofstream stn_out("stn_list", ios::out);

   listStations(stn_out);
}   

//
// List all available stations.
//
void
StationData::listStations (ofstream& os)
{
   for (int i = 0; i < m_stn.size(); i++)
   {
      os << m_stn[i].id;

      for (int k = 0; k < BL_SPACEDIM; k++)
      {
         os << "    " << m_stn[i].pos[k];
      }

      os << '\n';
   }
}

//
// Ascii read from restart file.
//
void
StationData::readStation (ifstream& is,
                          int       step)
{
    char buf[80];

    int nstation = 0;

    is >> m_nvar;
    is >> nstation;

    while (is.get() != '\n')
        ;

    if (nstation > 0)
    {
        is >> buf;

        while (is.get() != '\n')
            ;

        m_name = buf;

        m_stn.resize(nstation);

        for (int i = 0; i < m_stn.size(); i++)
        {
            is >> m_stn[i].id;

            for (int k = 0; k < BL_SPACEDIM; k++)
                is >> m_stn[i].pos[k];

            while (is.get() != '\n')
                ;
        }
        openFile(step);
    }
    //
    // re-write the "stn_list" file.
    //
    listStations();
}

static Real *data      = 0;
static int   data_size = 0;

void
StationData::report (int  level,
                     Real time,
                     int  nvar)
{
    if (nvar > data_size)
    {
        data      = new Real[nvar];
        data_size = nvar;
    }
    for (int i = 0; i < m_stn.size(); i++)
    {
        if (m_stn[i].level == level)
        {
            m_stn[i].grd->getData(m_stn[i].ix,data);
            //
            // Write data to output stream.
            //
#ifdef   STN_ASCII
            m_ofile << m_stn[i].id << ' ';
//	 m_ofile << time << " (";
            m_ofile << time << ' ';
/*       the following statements commented by ralph ferguson  14Feb96
	 for (int k = 0; k < BL_SPACEDIM; k++) {
         m_ofile << m_stn[i].pos[k];
         if (k == BL_SPACEDIM-1) {
         m_ofile << ") ";
         } else {
         m_ofile << ',';
         }
	 }
         */
//	 for (k = 0; k < nvar; k++) {
            for (int k = 0; k < nvar; k++)
            {
                m_ofile << data[k] << ' ';
            }
            m_ofile << '\n';
#else
            m_ofile.write((char*) &(m_stn[i].id), sizeof(int));
            m_ofile.write((char*) &time, sizeof(time));
            m_ofile.write((char*) m_stn[i].pos, BL_SPACEDIM*sizeof(Real));
            m_ofile.write((char*) data, nvar*sizeof(Real));
#endif
        }
    }
}

void
StationData::findGrid (const GridList* gl,
                       int             hi_lev,
                       Real            dx_lev[MAX_LEV][BL_SPACEDIM],
                       Real            prob_lo[BL_SPACEDIM])
{
    if (m_stn.size() <= 0)
        return;
    //
    // Flag all stations as not having a home.
    //
    for (int i = 0; i < m_stn.size(); i++)
    {
        m_stn[i].level = -1;
        m_stn[i].grd = 0;
    }
    //
    // Walk list of grids from highest level down.
    //
    for (int level = hi_lev; level >= 0; level--)
    {
        Grid* g;
        for (g = gl[level].first(); g != 0; g = gl[level].next(g))
        {
            //
            // Find all station points in this grid.
            //
            const Box& gbx = g->box();

            for (i = 0; i < m_stn.size(); i++)
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
