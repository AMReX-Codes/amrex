/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "BeamRelevant.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"
#include <iostream>
#include <cmath>
#include <limits>

using namespace amrex;

// constructor
BeamRelevant::BeamRelevant (std::string rd_name)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
    #if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "BeamRelevant reduced diagnostics does not work for RZ coordinate.");
    #endif

    // read beam name
    ParmParse pp(rd_name);
    pp.get("species",m_beam_name);

    // resize data array
    #if (AMREX_SPACEDIM == 3)
    //  0, 1, 2: mean x,y,z
    //  3, 4, 5: mean px,py,pz
    //        6: gamma
    //  7, 8, 9: rms x,y,z
    // 10,11,12: rms px,py,pz
    //       13: rms gamma
    // 14,15,16: emittance x,y,z
    m_data.resize(17,0.0);
    #elif (AMREX_SPACEDIM == 2)
    //     0, 1: mean x,z
    //  2, 3, 4: mean px,py,pz
    //        5: gamma
    //     6, 7: rms x,z
    //  8, 9,10: rms px,py,pz
    //       11: rms gamma
    //    12,13: emittance x,z
    m_data.resize(14,0.0);
    #endif

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            #if (AMREX_SPACEDIM == 3)
            ofs << "#";
            ofs << "[1]step";             ofs << m_sep;
            ofs << "[2]time(s)";          ofs << m_sep;
            ofs << "[3]x_mean(m)";        ofs << m_sep;
            ofs << "[4]y_mean(m)";        ofs << m_sep;
            ofs << "[5]z_mean(m)";        ofs << m_sep;
            ofs << "[6]px_mean(kg m/s)";  ofs << m_sep;
            ofs << "[7]py_mean(kg m/s)";  ofs << m_sep;
            ofs << "[8]pz_mean(kg m/s)";  ofs << m_sep;
            ofs << "[9]gamma_mean";       ofs << m_sep;
            ofs << "[10]x_rms(m)";        ofs << m_sep;
            ofs << "[11]y_rms(m)";        ofs << m_sep;
            ofs << "[12]z_rms(m)";        ofs << m_sep;
            ofs << "[13]px_rms(kg m/s)";  ofs << m_sep;
            ofs << "[14]py_rms(kg m/s)";  ofs << m_sep;
            ofs << "[15]pz_rms(kg m/s)";  ofs << m_sep;
            ofs << "[16]gamma_rms";       ofs << m_sep;
            ofs << "[17]emittance_x(m)";  ofs << m_sep;
            ofs << "[18]emittance_y(m)";  ofs << m_sep;
            ofs << "[19]emittance_z(m)";  ofs << std::endl;
            #elif (AMREX_SPACEDIM == 2)
            ofs << "#";
            ofs << "[1]step";             ofs << m_sep;
            ofs << "[2]time(s)";          ofs << m_sep;
            ofs << "[3]x_mean(m)";        ofs << m_sep;
            ofs << "[4]z_mean(m)";        ofs << m_sep;
            ofs << "[5]px_mean(kg m/s)";  ofs << m_sep;
            ofs << "[6]py_mean(kg m/s)";  ofs << m_sep;
            ofs << "[7]pz_mean(kg m/s)";  ofs << m_sep;
            ofs << "[8]gamma_mean";       ofs << m_sep;
            ofs << "[9]x_rms(m)";         ofs << m_sep;
            ofs << "[10]z_rms(m)";        ofs << m_sep;
            ofs << "[11]px_rms(kg m/s)";  ofs << m_sep;
            ofs << "[12]py_rms(kg m/s)";  ofs << m_sep;
            ofs << "[13]pz_rms(kg m/s)";  ofs << m_sep;
            ofs << "[14]gamma_rms";       ofs << m_sep;
            ofs << "[15]emittance_x(m)";  ofs << m_sep;
            ofs << "[16]emittance_z(m)";  ofs << std::endl;
            #endif
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that compute beam relevant quantities
void BeamRelevant::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get MultiParticleContainer class object
    auto & mypc = WarpX::GetInstance().GetPartContainer();

    // get number of species (int)
    auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    // inverse of speed of light squared
    auto inv_c2 = 1.0 / (PhysConst::c * PhysConst::c);

    // If 2D-XZ, p.pos(1) is z, rather than p.pos(2).
    #if (AMREX_SPACEDIM == 3)
    int index_z = 2;
    #elif (AMREX_SPACEDIM == 2)
    int index_z = 1;
    #endif

    // loop over species
    for (int i_s = 0; i_s < nSpecies; ++i_s)
    {

        // only beam species does
        if (species_names[i_s] != m_beam_name) { continue; }

        // get WarpXParticleContainer class object
        auto & myspc = mypc.GetParticleContainer(i_s);

        // get mass (Real)
        auto m = myspc.getMass();

        using PType = typename WarpXParticleContainer::SuperParticleType;

        // weight sum
        auto w_sum = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::w); });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum(w_sum);

        if (w_sum < std::numeric_limits<Real>::min() )
        {
            for (int i = 0; i < m_data.size(); ++i) { m_data[i] = 0.0; }
            return;
        }

        // x mean
        auto x_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(0) * p.rdata(PIdx::w) / w_sum; });

        #if (AMREX_SPACEDIM == 3)
        // y mean
        auto y_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(1) * p.rdata(PIdx::w) / w_sum; });
        #endif

        // z mean
        auto z_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.pos(index_z) * p.rdata(PIdx::w) / w_sum; });

        // ux mean
        auto ux_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::ux) * p.rdata(PIdx::w) / w_sum; });

        // uy mean
        auto uy_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uy) * p.rdata(PIdx::w) / w_sum; });

        // uz mean
        auto uz_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        { return p.rdata(PIdx::uz) * p.rdata(PIdx::w) / w_sum; });

        // gamma mean
        auto gm_mean = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto ux = p.rdata(PIdx::ux);
            auto uy = p.rdata(PIdx::uy);
            auto uz = p.rdata(PIdx::uz);
            auto us = ux*ux + uy*uy + uz*uz;
            return std::sqrt(1.0 + us*inv_c2) * p.rdata(PIdx::w) / w_sum;
        });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum(x_mean);
        #if (AMREX_SPACEDIM == 3)
        ParallelDescriptor::ReduceRealSum(y_mean);
        #endif
        ParallelDescriptor::ReduceRealSum(z_mean);
        ParallelDescriptor::ReduceRealSum(ux_mean);
        ParallelDescriptor::ReduceRealSum(uy_mean);
        ParallelDescriptor::ReduceRealSum(uz_mean);
        ParallelDescriptor::ReduceRealSum(gm_mean);

        // x rms
        auto x_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.pos(0)-x_mean) * (p.pos(0)-x_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        #if (AMREX_SPACEDIM == 3)
        // y rms
        auto y_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.pos(1)-y_mean) * (p.pos(1)-y_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });
        #endif

        // z rms
        auto z_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.pos(index_z)-z_mean) * (p.pos(index_z)-z_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        // ux rms
        auto ux_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.rdata(PIdx::ux)-ux_mean) * (p.rdata(PIdx::ux)-ux_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        // uy rms
        auto uy_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.rdata(PIdx::uy)-uy_mean) * (p.rdata(PIdx::uy)-uy_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        // uz rms
        auto uz_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.rdata(PIdx::uz)-uz_mean) * (p.rdata(PIdx::uz)-uz_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        // gamma rms
        auto gm_rms = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto ux = p.rdata(PIdx::ux);
            auto uy = p.rdata(PIdx::uy);
            auto uz = p.rdata(PIdx::uz);
            auto us = ux*ux + uy*uy + uz*uz;
            auto gm = std::sqrt(1.0 + us*inv_c2);
            auto a  = (gm - gm_mean) * (gm - gm_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        // x times ux
        auto xux = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.pos(0)-x_mean) * (p.rdata(PIdx::ux)-ux_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        #if (AMREX_SPACEDIM == 3)
        // y times uy
        auto yuy = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.pos(1)-y_mean) * (p.rdata(PIdx::uy)-uy_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });
        #endif

        // z times uz
        auto zuz = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto a = (p.pos(index_z)-z_mean) * (p.rdata(PIdx::uz)-uz_mean);
            return a * p.rdata(PIdx::w) / w_sum;
        });

        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
            ( x_rms, ParallelDescriptor::IOProcessorNumber());
        #if (AMREX_SPACEDIM == 3)
        ParallelDescriptor::ReduceRealSum
            ( y_rms, ParallelDescriptor::IOProcessorNumber());
        #endif
        ParallelDescriptor::ReduceRealSum
            ( z_rms, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum
            (ux_rms, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum
            (uy_rms, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum
            (uz_rms, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum
            (gm_rms, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::ReduceRealSum
            (   xux, ParallelDescriptor::IOProcessorNumber());
        #if (AMREX_SPACEDIM == 3)
        ParallelDescriptor::ReduceRealSum
            (   yuy, ParallelDescriptor::IOProcessorNumber());
        #endif
        ParallelDescriptor::ReduceRealSum
            (   zuz, ParallelDescriptor::IOProcessorNumber());

        // save data
        #if (AMREX_SPACEDIM == 3)
        m_data[0]  = x_mean;
        m_data[1]  = y_mean;
        m_data[2]  = z_mean;
        m_data[3]  = ux_mean * m;
        m_data[4]  = uy_mean * m;
        m_data[5]  = uz_mean * m;
        m_data[6]  = gm_mean;
        m_data[7]  = std::sqrt(x_rms);
        m_data[8]  = std::sqrt(y_rms);
        m_data[9]  = std::sqrt(z_rms);
        m_data[10] = std::sqrt(ux_rms * m);
        m_data[11] = std::sqrt(uy_rms * m);
        m_data[12] = std::sqrt(uz_rms * m);
        m_data[13] = std::sqrt(gm_rms);
        m_data[14] = std::sqrt(std::abs(x_rms*ux_rms-xux*xux)) / PhysConst::c;
        m_data[15] = std::sqrt(std::abs(y_rms*uy_rms-yuy*yuy)) / PhysConst::c;
        m_data[16] = std::sqrt(std::abs(z_rms*uz_rms-zuz*zuz)) / PhysConst::c;
        #elif (AMREX_SPACEDIM == 2)
        m_data[0]  = x_mean;
        m_data[1]  = z_mean;
        m_data[2]  = ux_mean * m;
        m_data[3]  = uy_mean * m;
        m_data[4]  = uz_mean * m;
        m_data[5]  = gm_mean;
        m_data[6]  = std::sqrt(x_rms);
        m_data[7]  = std::sqrt(z_rms);
        m_data[8]  = std::sqrt(ux_rms * m);
        m_data[9]  = std::sqrt(uy_rms * m);
        m_data[10] = std::sqrt(uz_rms * m);
        m_data[11] = std::sqrt(gm_rms);
        m_data[12] = std::sqrt(std::abs(x_rms*ux_rms-xux*xux)) / PhysConst::c;
        m_data[13] = std::sqrt(std::abs(z_rms*uz_rms-zuz*zuz)) / PhysConst::c;
        #endif

    }
    // end loop over species

}
// end void BeamRelevant::ComputeDiags
