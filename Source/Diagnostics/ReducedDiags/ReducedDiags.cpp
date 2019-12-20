#include "ReducedDiags.H"
#include "WarpX.H"
#include "AMReX_ParmParse.H"
#include <iomanip>

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

    /// get dt
    auto dt = WarpX::GetInstance().getdt(0);

    /// write step
    ofs << step+1;

    ofs << m_sep;

    /// set precision
    ofs << std::fixed << std::setprecision(14) << std::scientific;

    /// write time
    ofs << (step+1)*dt;

    /// loop over data size and write
    for (int i = 0; i < m_data.size(); ++i)
    {
        ofs << m_sep;
        //ofs << std::setprecision(16) << std::scientific;
        ofs << m_data[i];
    }
    ///< end loop over data size

    /// end line
    ofs << std::endl;

}
///< end ReducedDiags::WriteToFile
