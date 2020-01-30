/* Copyright 2019-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ReducedDiags.H"
#include "WarpX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Utility.H"
#include <iomanip>

using namespace amrex;

// constructor
ReducedDiags::ReducedDiags (std::string rd_name)
{

    m_rd_name = rd_name;

    ParmParse pp(m_rd_name);

    // read path
    pp.query("path", m_path);

    // read extension
    pp.query("extension", m_extension);

    // creater folder
    if (!UtilCreateDirectory(m_path, 0755))
    { CreateDirectoryFailed(m_path); }

    // check if it is a restart run
    std::string restart_chkfile = "";
    ParmParse pp_amr("amr");
    pp_amr.query("restart", restart_chkfile);
    m_IsNotRestart = restart_chkfile.empty();

    // replace / create output file
    if ( m_IsNotRestart ) // not a restart
    {
        std::ofstream ofs;
        ofs.open(m_path+m_rd_name+"."+m_extension, std::ios::trunc);
        ofs.close();
    }

    // read reduced diags frequency
    pp.query("frequency", m_freq);

    // read separator
    pp.query("separator", m_sep);

}
// end constructor

// write to file function
void ReducedDiags::WriteToFile (int step) const
{

    // open file
    std::ofstream ofs;
    ofs.open(m_path + m_rd_name + "." + m_extension,
        std::ofstream::out | std::ofstream::app);

    // write step
    ofs << step+1;

    ofs << m_sep;

    // set precision
    ofs << std::fixed << std::setprecision(14) << std::scientific;

    // write time
    ofs << WarpX::GetInstance().gett_new(0);

    // loop over data size and write
    for (int i = 0; i < m_data.size(); ++i)
    {
        ofs << m_sep;
        ofs << m_data[i];
    }
    // end loop over data size

    // end line
    ofs << std::endl;

    // close file
    ofs.close();

}
// end ReducedDiags::WriteToFile
