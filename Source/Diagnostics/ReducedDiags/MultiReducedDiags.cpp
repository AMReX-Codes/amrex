#include "MultiReducedDiags.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Utility.H"
#include <fstream>

/// constructor
MultiReducedDiags::MultiReducedDiags ()
{

    /// read reduced diags flag, names, and path
    amrex::ParmParse pp("warpx");
    pp.query("plot_reduced_diags", m_plot_rd);
    if ( m_plot_rd == 0 ) { return; }
    pp.getarr("reduced_diags_names", m_rd_names);
    pp.query("reduced_diags_path", m_path);

    /// creater folder
    if (!amrex::UtilCreateDirectory(m_path, 0755))
    { amrex::CreateDirectoryFailed(m_path); }

    /// resize
    m_multi_rd.resize(m_rd_names.size());

    /// loop over all reduced diags
    for (int i_rd = 0; i_rd < m_rd_names.size(); ++i_rd)
    {

        amrex::ParmParse pp(m_rd_names[i_rd]);

        /// read reduced diags type
        std::string rd_type;
        pp.query("type", rd_type);

        /// replace old and open output file
        std::ofstream ofs;
        ofs.open(m_path+m_rd_names[i_rd]+".txt", std::ios::trunc);

        /// match diags
        if (rd_type.compare("ParticleMeanEnergy") == 0)
        {
            m_multi_rd[i_rd].reset
                ( new ParticleMeanEnergy(m_rd_names[i_rd], ofs));
        }
        ///< end if match diags

        /// close file
        ofs.close();

    }
    ///< end loop over all reduced diags

}
///< end constructor

/// destructor
MultiReducedDiags::~MultiReducedDiags ()
{}
///< end destructor

/// call functions to compute diags
void MultiReducedDiags::ComputeDiags (int step)
{
    /// loop over all reduced diags
    for (int i_rd = 0; i_rd < m_rd_names.size(); ++i_rd)
    {
        m_multi_rd[i_rd] -> ComputeDiags(step);
    }
    ///< end loop over all reduced diags
}
///< end void MultiReducedDiags::ComputeDiags

/// funciton to write data
void MultiReducedDiags::WriteToFile (int step)
{

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) { return; }

    /// loop over all reduced diags
    for (int i_rd = 0; i_rd < m_rd_names.size(); ++i_rd)
    {

        /// Judge if the diags should be done
        if ( (step+1) % m_multi_rd[i_rd]->m_freq != 0 ) { return; }

        /// open file
        std::ofstream ofs;
        ofs.open(m_path+m_rd_names[i_rd]+".txt",
            std::ofstream::out | std::ofstream::app);

        /// call the write to file function
        m_multi_rd[i_rd]->WriteToFile(step, ofs);

        /// close file
        ofs.close();

    }
    ///< end loop over all reduced diags
}
///< end void MultiReducedDiags::WriteToFile
