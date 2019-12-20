#include "ReducedDiags.H"
#include "WarpX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"

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
    amrex::Print(ofs) << step+1;

    amrex::Print(ofs) << m_sep;

    /// write time
    amrex::Print(ofs) << (step+1)*dt;

    /// loop over data size and write
    for (int i = 0; i < m_data.size(); ++i)
    {
        amrex::Print(ofs) << m_sep;
        amrex::Print(ofs) << m_data[i];
    }
    ///< end loop over data size

    /// end line
    amrex::Print(ofs) << std::endl;

}
///< end ReducedDiags::WriteToFile
