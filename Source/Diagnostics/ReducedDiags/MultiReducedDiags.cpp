#include "MultiReducedDiags.H"
#include "AMReX_ParmParse.H"

using namespace amrex;

/** constructor */
MultiReducedDiags::MultiReducedDiags ()
{

    /** read reduced diags names and types */
    ParmParse pp("warpx");
    pp.getarr("reduced_diags_names", m_reduced_diags_names);
    m_n_diags = m_reduced_diags_names.size();
    m_multi_rd.resize(m_n_diags);

    /** loop over all reduced diags */
    for (int i_rd = 0; i_rd < m_n_diags; ++i_rd)
    {

        ParmParse pp(m_reduced_diags_names[i_rd]);

        /** read reduced diags type */
        std::vector<std::string> reduced_diags_type;
        pp.getarr("type", reduced_diags_type);

        /** read reduced diags frequency */
        int reduced_diags_freq;
        pp.query("frequency", reduced_diags_freq);

        /** initialize the diags */
        if (reduced_diags_type[0].compare("ParticleKineticEnergy") == 0)
        {
            m_multi_rd[i_rd].reset( new ParticleKineticEnergy() );
        }

    } // end loop over all reduced diags

} // end constructor

/** destructor */
MultiReducedDiags::~MultiReducedDiags ()
{}

void MultiReducedDiags::ComputeDiags (int step)
{}

void MultiReducedDiags::WriteToFile (int step)
{}
