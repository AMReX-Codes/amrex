#include "ReducedDiags.H"
#include "AMReX_ParmParse.H"

/// constructor
ReducedDiags::ReducedDiags (std::string rd_name)
{
    /// read reduced diags frequency
    amrex::ParmParse pp(rd_name);
    pp.query("frequency", m_freq);
}
///< end constructor

/// destructor
ReducedDiags::~ReducedDiags ()
{}
///< end destructor

/// write to file function
void ReducedDiags::WriteToFile (int step, std::ofstream & ofs)
{

    /// write step
    ofs << step+1;

    /// loop over data size and write
    for (int i = 0; i < m_data.size(); ++i)
    {
        ofs << ",";
        ofs << m_data[i];
    }
    ///< end loop over data size

}
///< end ReducedDiags::WriteToFile
