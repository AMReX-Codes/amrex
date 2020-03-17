/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleHistogram.H"
#include "WarpX.H"
#include "Utils/WarpXUtil.H"
#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>
#include <limits>

using namespace amrex;

struct NormalizationType {
    enum {
        no_normalization = 0,
        unity_particle_weight,
        max_to_unity,
        area_to_unity
    };
};

// constructor
ParticleHistogram::ParticleHistogram (std::string rd_name)
: ReducedDiags{rd_name}
{

    ParmParse pp(rd_name);

    // read species
    std::string selected_species_name;
    pp.get("species",selected_species_name);

    // read bin parameters
    pp.get("bin_number",m_bin_num);
    pp.get("bin_max",   m_bin_max);
    pp.get("bin_min",   m_bin_min);
    m_bin_size = (m_bin_max - m_bin_min) / m_bin_num;

    // read histogram function
    std::string function_string = "";
    Store_parserString(pp,"histogram_function(t,x,y,z,ux,uy,uz)",
                       function_string);
    m_parser.reset(new ParserWrapper<m_nvars>(
        makeParser(function_string,{"t","x","y","z","ux","uy","uz"})));

    // read normalization type
    std::string norm_string = "default";
    pp.query("normalization",norm_string);

    // set normalization type
    if ( norm_string == "default" ) {
        m_norm = NormalizationType::no_normalization;
    } else if ( norm_string == "unity_particle_weight" ) {
        m_norm = NormalizationType::unity_particle_weight;
    } else if ( norm_string == "max_to_unity" ) {
        m_norm = NormalizationType::max_to_unity;
    } else if ( norm_string == "area_to_unity" ) {
        m_norm = NormalizationType::area_to_unity;
    } else {
        Abort("Unknown ParticleHistogram normalization type.");
    }

    // get MultiParticleContainer class object
    auto & mypc = WarpX::GetInstance().GetPartContainer();
    // get species names (std::vector<std::string>)
    auto const species_names = mypc.GetSpeciesNames();
    // select species
    for ( int i = 0; i < mypc.nSpecies(); ++i )
    {
        if ( selected_species_name == species_names[i] ){
            m_selected_species_id = i;
        }
    }
    // if m_selected_species_id is not modified
    if ( m_selected_species_id == -1 ){
        Abort("Unknown species for ParticleHistogram reduced diagnostic.");
    }

    // resize data array
    m_data.resize(m_bin_num,0.0_rt);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            for (int i = 0; i < m_bin_num; ++i)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(3+i) + "]";
                Real b = m_bin_min + m_bin_size*(Real(i)+0.5_rt);
                ofs << "bin" + std::to_string(1+i)
                             + "=" + std::to_string(b) + "()";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes the histogram
void ParticleHistogram::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) return;

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get time at level 0
    auto const t = warpx.gett_new(0);

    // get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    // get WarpXParticleContainer class object
    auto const & myspc = mypc.GetParticleContainer(m_selected_species_id);

    using PType = typename WarpXParticleContainer::SuperParticleType;

    // get parser
    ParserWrapper<m_nvars> *fun_partparser = m_parser.get();

    // declare local variables
    Real const bin_min  = m_bin_min;
    Real const bin_size = m_bin_size;
    const bool is_unity_particle_weight =
        (m_norm == NormalizationType::unity_particle_weight) ? true : false;

    for ( int i = 0; i < m_bin_num; ++i )
    {
        // compute the histogram
        m_data[i] = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto const w  = p.rdata(PIdx::w);
            auto const x  = p.pos(0);
            auto const y  = p.pos(1);
            auto const z  = p.pos(2);
            auto const ux = p.rdata(PIdx::ux)/PhysConst::c;
            auto const uy = p.rdata(PIdx::uy)/PhysConst::c;
            auto const uz = p.rdata(PIdx::uz)/PhysConst::c;
            auto const f = (*fun_partparser)(t,x,y,z,ux,uy,uz);
            auto const f1 = bin_min + bin_size*i;
            auto const f2 = bin_min + bin_size*(i+1);
            if ( f > f1 && f < f2 ) {
                if ( is_unity_particle_weight ) return 1.0_rt;
                else return w;
            } else return 0.0_rt;
        });
    }
    // reduced sum over mpi ranks
    ParallelDescriptor::ReduceRealSum
        (m_data.data(), m_data.size(), ParallelDescriptor::IOProcessorNumber());

    // normalize the maximum value to be one
    if ( m_norm == NormalizationType::max_to_unity )
    {
        Real f_max = 0.0_rt;
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( m_data[i] > f_max ) f_max = m_data[i];
        }
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( f_max > std::numeric_limits<Real>::min() ) m_data[i] /= f_max;
        }
        return;
    }

    // normalize the area (integral) to be one
    if ( m_norm == NormalizationType::area_to_unity )
    {
        Real f_area = 0.0_rt;
        for ( int i = 0; i < m_bin_num; ++i )
        {
            f_area += m_data[i] * m_bin_size;
        }
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( f_area > std::numeric_limits<Real>::min() ) m_data[i] /= f_area;
        }
        return;
    }

}
// end void ParticleHistogram::ComputeDiags
