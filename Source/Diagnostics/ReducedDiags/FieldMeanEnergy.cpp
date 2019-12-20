#include "FieldMeanEnergy.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <cmath>

/// constructor
FieldMeanEnergy::FieldMeanEnergy (
std::string rd_name, std::ofstream & ofs )
: ReducedDiags{rd_name}
{
    // write header row
    ofs << "#" << "step" << m_sep << "time(s)" << m_sep;
    ofs << "total(J)" << m_sep << "E(J)" << m_sep << "B(J)";
    ofs << std::endl;
}
///< end constructor

/// destructor
FieldMeanEnergy::~FieldMeanEnergy ()
{}
///< end destructor

/// function that computes kinetic energy
void FieldMeanEnergy::ComputeDiags (int step)
{

    /// Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    /// get WarpX class object
    auto & warpx = WarpX::GetInstance();

    /// resize data array
    /// the extra one is for total energy
    m_data.resize(3,0.0);

    /// declare total energy variable
    amrex::Real Etot = 0.0;

    #pragma omp parallel reduction(+:EK, np)

    /// compute total EK
    EKtot += EK;

    /// compute total np
    nptot += np;

    /// save EK to m_data
    if ( np > 0 )
    { m_data[i_s+1] = EK / np; }
    else
    { m_data[i_s+1] = 0.0; }

    /// save total EK
    if ( nptot > 0 )
    { m_data[0] = EKtot / nptot; }
    else
    { m_data[0] = 0.0; }

}
///< end void FieldMeanEnergy::ComputeDiags
