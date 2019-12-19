#include "MultiReducedDiags.H"
#include "AMReX_ParmParse.H"
#include <fstream>

/// constructor
MultiReducedDiags::MultiReducedDiags ()
{

    /// read reduced diags names
    amrex::ParmParse pp("warpx");
    pp.getarr("reduced_diags_names", m_rd_names);
    m_multi_rd.resize(m_rd_names.size());

    /// loop over all reduced diags
    for (int i_rd = 0; i_rd < m_rd_names.size(); ++i_rd)
    {

        amrex::ParmParse pp(m_rd_names[i_rd]);

        /// read reduced diags type
        std::string rd_type;
        pp.query("type", rd_type);

        /// match diags
        if (rd_type.compare("ParticleKineticEnergy") == 0)
        {
            /// replace old and open output file
            std::string dirname = "./";
            std::ofstream ofs;
            ofs.open(dirname+m_rd_names[i_rd]+".txt", std::ios::trunc);
            /// declare
            m_multi_rd[i_rd].reset
                ( new ParticleKineticEnergy(m_rd_names[i_rd], ofs));
            /// close file
            ofs.close();
        }
        ///< end if match diags

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
        std::string dirname = "./";
        std::ofstream ofs;
        ofs.open(dirname+m_rd_names[i_rd]+".txt",
            std::ofstream::out | std::ofstream::app);

        /// call the write to file function
        m_multi_rd[i_rd]->WriteToFile(step, ofs);

        /// close file
        ofs << std::endl;
        ofs.close();

    }
    ///< end loop over all reduced diags
}
///< end void MultiReducedDiags::WriteToFile
