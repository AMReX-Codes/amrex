#include "MultiReducedDiags.H"
#include "AMReX_ParmParse.H"

using namespace amrex;

/** constructor */
MultiReducedDiags::MultiReducedDiags ()
{

    /** read reduced diags names */
    ParmParse pp("warpx");
    pp.getarr("reduced_diags_names", m_rd_names);
    m_multi_rd.resize(m_rd_names.size());

    /** loop over all reduced diags */
    for (int i_rd = 0; i_rd < m_rd_names.size(); ++i_rd)
    {

        ParmParse pp(m_rd_names[i_rd]);

        /** read reduced diags type */
        std::string rd_type;
        pp.query("type", rd_type);

        /** initialize the diags */
        if (rd_type.compare("ParticleKineticEnergy") == 0)
        {
            m_multi_rd[i_rd].reset
                ( new ParticleKineticEnergy(m_rd_names[i_rd]) );
        }

    } // end loop over all reduced diags

} // end constructor

/** destructor */
MultiReducedDiags::~MultiReducedDiags ()
{}

void MultiReducedDiags::ComputeDiags (int step)
{
    for (int i_rd = 0; i_rd < m_rd_names.size(); ++i_rd)
    {
        m_multi_rd[i_rd] -> ComputeDiags(step);
    }
}

void MultiReducedDiags::WriteToFile (int step)
{}
