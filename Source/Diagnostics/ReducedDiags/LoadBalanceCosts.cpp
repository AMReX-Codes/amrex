/* Copyright 2019-2020 Michael Rowan, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "LoadBalanceCosts.H"
#include "WarpX.H"


using namespace amrex;

// constructor
LoadBalanceCosts::LoadBalanceCosts (std::string rd_name)
    : ReducedDiags{rd_name}
{

}

// function that gathers costs
void LoadBalanceCosts::ComputeDiags (int step)
{
    // get WarpX class object
    auto& warpx = WarpX::GetInstance();

    // for now, costs reduced diagnostic only works with `costs_heuristic`,
    // but this will work for timer based costs as well after changing from
    // multifab to vector
    const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cost_heuristic != nullptr,
        "LoadBalanceCosts reduced diagnostic works only with heuristic costs update.");

    // judge if the diags should be done
    // costs is initialized only if we're doing load balance
    if ( ((step+1) % m_freq != 0) || warpx.get_load_balance_int() < 1 ) { return; }

    // get number of boxes over all levels
    auto nLevels = warpx.finestLevel() + 1;
    int nBoxes = 0;
    for (int lev = 0; lev < nLevels; ++lev)
    {
        const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(lev);
        nBoxes += cost_heuristic->size();
    }

    // keep track of the max number of boxes, this is needed later on to fill
    // the jagged array (in case each step does not have the same number of boxes)
    m_nBoxesMax = std::max(m_nBoxesMax, nBoxes);

    // resize and clear data array
    const size_t dataSize = m_nDataFields*nBoxes;
    m_data.resize(dataSize, 0.0);
    m_data.assign(dataSize, 0.0);

    // read in WarpX costs_heuristic to local copy
    amrex::Vector<std::unique_ptr<amrex::Vector<amrex::Real> > > costs_heuristic;
    if (warpx.load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic)
    {
        costs_heuristic.resize(nLevels);
        for (int lev = 0; lev < nLevels; ++lev)
        {
            costs_heuristic[lev].reset(new amrex::Vector<Real>);
            const int nBoxesLev = warpx.getCostsHeuristic(lev)->size();
            costs_heuristic[lev]->resize(nBoxesLev);
        }
        warpx.ComputeCostsHeuristic(costs_heuristic);
    }

    // keeps track of correct index in array over all boxes on all levels
    int shift = 0;

    // save data
    for (int lev = 0; lev < nLevels; ++lev)
    {
        const amrex::DistributionMapping& dm = warpx.DistributionMap(lev);
        const MultiFab & Ex = warpx.getEfield(lev,0);
        for (MFIter mfi(Ex, false); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            m_data[shift + mfi.index()*m_nDataFields + 0] = (*costs_heuristic[lev])[mfi.index()];
            m_data[shift + mfi.index()*m_nDataFields + 1] = dm[mfi.index()];
            m_data[shift + mfi.index()*m_nDataFields + 2] = lev;
            m_data[shift + mfi.index()*m_nDataFields + 3] = tbx.loVect()[0];
            m_data[shift + mfi.index()*m_nDataFields + 4] = tbx.loVect()[1];
            m_data[shift + mfi.index()*m_nDataFields + 5] = tbx.loVect()[2];
        }

        // we looped through all the boxes on level lev, update the shift index
        shift += m_nDataFields*(costs_heuristic[lev]->size());
    }

    // parallel reduce to IO proc and get data over all procs
    ParallelDescriptor::ReduceRealSum(m_data.data(), m_data.size(), ParallelDescriptor::IOProcessorNumber());

    /* m_data now contains up-to-date values for:
     *  [[cost, proc, lev, i_low, j_low, k_low] of box 0 at level 0,
     *   [cost, proc, lev, i_low, j_low, k_low] of box 1 at level 0,
     *   [cost, proc, lev, i_low, j_low, k_low] of box 2 at level 0,
     *   ...
     *   [cost, proc, lev, i_low, j_low, k_low] of box 0 at level 1,
     *   [cost, proc, lev, i_low, j_low, k_low] of box 1 at level 1,
     *   [cost, proc, lev, i_low, j_low, k_low] of box 2 at level 1,
     *   ......] */

}

// write to file function for cost
void LoadBalanceCosts::WriteToFile (int step) const
{
    ReducedDiags::WriteToFile(step);

    // get WarpX class object
    auto& warpx = WarpX::GetInstance();

    if (!ParallelDescriptor::IOProcessor()) return;

    // final step is a special case, fill jagged array with NaN
    if (step == (warpx.maxStep() - (warpx.maxStep()%m_freq) - 1 ))
    {
        // open tmp file to copy data
        std::string fileTmpName = m_path + m_rd_name + ".tmp." + m_extension;
        std::ofstream ofs(fileTmpName, std::ofstream::out);
        // write header row
        // for each box on each level we saved 6 data fields: [cost, proc, lev, i_low, j_low, k_low])
        ofs << "#";
        ofs << "[1]step()";
        ofs << m_sep;
        ofs << "[2]time(s)";

        for (int boxNumber=0; boxNumber<m_nBoxesMax; ++boxNumber)
        {
            ofs << m_sep;
            ofs << "[" + std::to_string(3 + m_nDataFields*boxNumber) + "]";
            ofs << "cost_box_"+std::to_string(boxNumber)+"()";
            ofs << m_sep;
            ofs << "[" + std::to_string(4 + m_nDataFields*boxNumber) + "]";
            ofs << "proc_box_"+std::to_string(boxNumber)+"()";
            ofs << m_sep;
            ofs << "[" + std::to_string(5 + m_nDataFields*boxNumber) + "]";
            ofs << "lev_box_"+std::to_string(boxNumber)+"()";
            ofs << m_sep;
            ofs << "[" + std::to_string(6 + m_nDataFields*boxNumber) + "]";
            ofs << "i_low_box_"+std::to_string(boxNumber)+"()";
            ofs << m_sep;
            ofs << "[" + std::to_string(7 + m_nDataFields*boxNumber) + "]";
            ofs << "j_low_box_"+std::to_string(boxNumber)+"()";
            ofs << m_sep;
            ofs << "[" + std::to_string(8 + m_nDataFields*boxNumber) + "]";
            ofs << "k_low_box_"+std::to_string(boxNumber)+"()";
        }
        ofs << std::endl;

        // open the data-containing file
        std::string fileDataName = m_path + m_rd_name + "." + m_extension;
        std::ifstream ifs(fileDataName, std::ifstream::in);

        // Fill in the tmp costs file with data, padded with NaNs
        for (std::string lineIn; std::getline(ifs, lineIn);)
        {
            // count the elements in the input line
            int cnt = 0;
            std::stringstream ss(lineIn);
            std::string token;

            while (std::getline(ss, token, m_sep[0]))
            {
                cnt += 1;
                if (ss.peek() == m_sep[0]) ss.ignore();
            }

            // 2 columns for step, time; then nBoxes*nDatafields columns for data;
            // then fill the remaining columns (i.e., up to 2 + m_nBoxesMax*m_nDataFields)
            // with NaN, so the array is not jagged
            ofs << lineIn;
            for (int i=0; i<(m_nBoxesMax*m_nDataFields - (cnt - 2)); ++i)
            {
                ofs << m_sep << "NaN";
            }
            ofs << std::endl;
        }

        // close files
        ifs.close();
        ofs.close();

        // remove the original, rename tmp file
        std::remove(fileDataName.c_str());
        std::rename(fileTmpName.c_str(), fileDataName.c_str());
    }
}
