
#include <winstd.H>
#include <cstdio>
#include <cstring>

#include <AmrLevel.H>
#include <ParmParse.H>
#include <StationData.H>
#include <Utility.H>

StationRec::StationRec ()
{
    D_TERM(pos[0],=pos[1],=pos[2]) = -1;
    id = level = grd = -1;
    own = false;
}

StationData::~StationData ()
{
    m_ofile.close();
}

void
StationData::init (const PArray<AmrLevel>& levels, const int finestlevel)
{
    //
    // ParmParse variables:
    //
    //   StationData.vars     -- Names of StateData components to output
    //   StationData.coord    -- BL_SPACEDIM array of Reals
    //   StationData.coord    -- the next one
    //   StationData.coord    -- ditto ...
    //
    ParmParse pp("StationData");

    if (pp.contains("vars"))
    {
        const int N = pp.countval("vars");

        m_vars.resize(N);
        m_IsDerived.resize(N);
        m_typ.resize(N);
        m_ncomp.resize(N);

        for (int i = 0; i < N; i++)
        {
            pp.get("vars", m_vars[i], i);

            m_IsDerived[i] = false;

            if (!AmrLevel::isStateVariable(m_vars[i],m_typ[i],m_ncomp[i]))
            {
                if (AmrLevel::get_derive_lst().canDerive(m_vars[i]))
                {
                    m_IsDerived[i] = true;
                    m_ncomp[i] = 0;
                }
                else
                {
                    std::cerr << "StationData::init(): `"
                              << m_vars[i]
                              << "' is not a state or derived variable\n";
                    BoxLib::Abort();
                }
            }
        }
    }

    if (m_vars.size() > 0 && pp.contains("coord"))
    {
        static int identifier;

        Array<Real> data(BL_SPACEDIM);

        const int N  = pp.countname("coord");

        m_stn.resize(N);

        for (int k = 0; k < N; k++)
        {
            pp.getktharr("coord", k, data, 0, BL_SPACEDIM);

            if(Geometry::ProbDomain().contains(data.dataPtr())) {

              D_TERM(m_stn[k].pos[0] = data[0];,
                     m_stn[k].pos[1] = data[1];,
                     m_stn[k].pos[2] = data[2];);


	    } else {
	      RealBox rbPD(Geometry::ProbDomain());
              if(ParallelDescriptor::IOProcessor()) {
	        std::cout << "**** Error in StationData:  coord outside probdomain: "
	                  << " pd  coord = " << rbPD << "   ";
                          for(int d(0); d < BL_SPACEDIM; ++d) {
			    std::cout << data[d] << "  ";
	                  }
	        std::cout << std::endl;
	        std::cout << "     Moving point to the domain center." << std::endl;
	      }
              D_TERM(m_stn[k].pos[0] = rbPD.lo(0) + (0.5 * rbPD.length(0));,
                     m_stn[k].pos[1] = rbPD.lo(1) + (0.5 * rbPD.length(1));,
                     m_stn[k].pos[2] = rbPD.lo(2) + (0.5 * rbPD.length(2));)
	    }
            BL_ASSERT(Geometry::ProbDomain().contains(data.dataPtr()));

            m_stn[k].id = identifier++;
        }
    }


    // adjust the positions so they are not exactly on a grid line
    // cannot use [levels.size()-1] because of level capping
    const Real *fineDX = levels[finestlevel].Geom().CellSize();
    for(int i(0); i < m_stn.size(); ++i) {
      for(int d(0); d < BL_SPACEDIM; ++d) {
        Real tempPos(m_stn[i].pos[d]), dx(fineDX[d]);
        Real dxEps(0.1 * dx);
	long nCells(fabs(tempPos/dx));
	Real rR = fabs(fabs(tempPos) - (nCells * dx));
	if(rR < dxEps) {
	  tempPos += dxEps;
	  if(tempPos > Geometry::ProbDomain().hi(d)) {
	    tempPos -= 2.0 * dxEps;
	  }
          if(ParallelDescriptor::IOProcessor()) {
	    std::cout << "***** Warning:  adjusting station coord from "
	              << m_stn[i].pos[d] << "  to  " << tempPos << std::endl;
	  }
	  m_stn[i].pos[d] = tempPos;
	}
      }
      if(Geometry::ProbDomain().contains(m_stn[i].pos) == false) {
        BoxLib::Abort("Bad Stations coord.");
      }
    }

    if (m_vars.size() > 0)
    {
        //
        // Make "Stations/" directory.
        //
        // Only the I/O processor makes the directory if it doesn't exist.
        //
        static const std::string Dir("Stations");

        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(Dir, 0755))
                BoxLib::CreateDirectoryFailed(Dir);
        //
        // Everyone must wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Open the data file.
        //
        const int MyProc = ParallelDescriptor::MyProc();

        std::string datafile = BoxLib::Concatenate("Stations/stn_CPU_", MyProc, 4);

        m_ofile.open(datafile.c_str(), std::ios::out|std::ios::app);

        m_ofile.precision(30);

        BL_ASSERT(!m_ofile.bad());
        //
        // Output the list of stations.
        //
        std::ofstream os("Stations/Station.List", std::ios::out);

        for (int i = 0; i < m_stn.size(); i++)
        {
            os << m_stn[i].id;

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                os << '\t' << m_stn[i].pos[k];
            }

            os << '\n';
        }
    }
}

void
StationData::report (Real            time,
                     int             level,
                     const AmrLevel& amrlevel)
{
    if (m_stn.size() <= 0)
        return;

    static Array<Real> data;

    const int N = m_vars.size();

    data.resize(N);

    Array<MultiFab*> mfPtrs(N, 0);
    int nGhost = 0;
    for (int iVar = 0; iVar < m_vars.size(); ++iVar)
    {
        if (m_IsDerived[iVar])
        {
            mfPtrs[iVar] =
                const_cast<AmrLevel&>(amrlevel).derive(m_vars[iVar], time, nGhost);
        }
    }

    for (int i = 0; i < m_stn.size(); i++)
    {
        if (m_stn[i].own && m_stn[i].level == level)
        {
            //
            // Fill in the data vector.
            //
            for (int j = 0; j < N; j++)
            {
                if(m_IsDerived[j])
                {
                    const MultiFab* mf = mfPtrs[j];

                    BL_ASSERT(mf != NULL);
                    BL_ASSERT(mf->DistributionMap() ==
                              amrlevel.get_new_data(0).DistributionMap());
                    BL_ASSERT(mf->DistributionMap()[m_stn[i].grd] == ParallelDescriptor::MyProc());
                    //
                    // Find IntVect so we can index into FAB.
                    // We want to use Geometry::CellIndex().
                    // Must adjust the position to account for NodeCentered-ness.
                    //
                    IndexType ityp = 
                        amrlevel.get_derive_lst().get(m_vars[j])->deriveType();

                    Real pos[BL_SPACEDIM];

                    D_TERM(pos[0] = m_stn[i].pos[0] + .5 * ityp[0];,
                           pos[1] = m_stn[i].pos[1] + .5 * ityp[1];,
                           pos[2] = m_stn[i].pos[2] + .5 * ityp[2];);

                    IntVect idx = amrlevel.Geom().CellIndex(&pos[0]);

                    data[j] = (*mf)[m_stn[i].grd](idx,m_ncomp[j]);
                }
                else
                {
                    const MultiFab& mf = amrlevel.get_new_data(m_typ[j]);

                    BL_ASSERT(mf.nComp() > m_ncomp[j]);
                    BL_ASSERT(mf.DistributionMap()[m_stn[i].grd] == ParallelDescriptor::MyProc());
                    //
                    // Find IntVect so we can index into FAB.
                    // We want to use Geometry::CellIndex().
                    // Must adjust the position to account for NodeCentered-ness.
                    //
                    IndexType ityp = amrlevel.get_desc_lst()[m_typ[j]].getType();

                    Real pos[BL_SPACEDIM];

                    D_TERM(pos[0] = m_stn[i].pos[0] + .5 * ityp[0];,
                           pos[1] = m_stn[i].pos[1] + .5 * ityp[1];,
                           pos[2] = m_stn[i].pos[2] + .5 * ityp[2];);

                    IntVect idx = amrlevel.Geom().CellIndex(&pos[0]);

                    data[j] = mf[m_stn[i].grd](idx,m_ncomp[j]);
                }
            }
            //
            // Write data to output stream.
            //
            m_ofile << m_stn[i].id << ' ' << time << ' ';

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                m_ofile << m_stn[i].pos[k] << ' ';
            }
            for (int k = 0; k < N; k++)
            {
                m_ofile << data[k] << ' ';
            }
            m_ofile << '\n';
        }
    }

    m_ofile.flush();

    BL_ASSERT(!m_ofile.bad());

    for (int iVarD = 0; iVarD < m_vars.size(); ++iVarD)
    {
        if (m_IsDerived[iVarD])
        {
            delete mfPtrs[iVarD];
        }
    }
}

void
StationData::findGrid (const PArray<AmrLevel>& levels,
                       const Array<Geometry>&  geoms)
{
    BL_ASSERT(geoms.size() == levels.size());

    if (m_stn.size() <= 0)
        return;
    //
    // Flag all stations as not having a home.
    //
    for (int i = 0; i < m_stn.size(); i++)
        m_stn[i].level = -1;
    //
    // Find level and grid owning the data.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    for (int level = levels.size()-1; level >= 0; level--)
    {
        if (levels.defined(level))
        {
            Array<RealBox> boxes(levels[level].boxArray().size());

            for (int i = 0; i < boxes.size(); i++)
            {
                boxes[i] = RealBox(levels[level].boxArray()[i],
                                   levels[level].Geom().CellSize(),
                                   levels[level].Geom().ProbLo());
            }


            MultiFab mf(levels[level].boxArray(),1,0,Fab_noallocate);

            BL_ASSERT(mf.boxArray().size() == boxes.size());

            for (int i = 0; i < m_stn.size(); i++)
            {
                if (m_stn[i].level < 0)
                {
                    for (int j = 0; j < boxes.size(); j++)
                    {
                        if (boxes[j].contains(&m_stn[i].pos[0]))
                        {
                            m_stn[i].grd   = j;
                            m_stn[i].own   = (mf.DistributionMap()[j] == MyProc);
                            m_stn[i].level = level;

                            break;
                        }
                    }
                }
            }
        }
    }
}
