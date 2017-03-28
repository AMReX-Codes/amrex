
#include <AMReX_ParmParse.H>

#include <WarpX.H>
#include <WarpXConst.H>

using namespace amrex;

const int debug_lb = 0;

void
WarpX::RegridBaseLevel ()
{
    int lev = 0;

    Array<std::unique_ptr<MultiFab>> old_current;
    Array<std::unique_ptr<MultiFab>> old_Efield;
    Array<std::unique_ptr<MultiFab>> old_Bfield;

    old_current.resize(3);
    old_Efield.resize(3);
    old_Bfield.resize(3);

    const IntVect& nodalflag = IntVect::TheUnitVector();

    MFInfo info;
    info.SetNodal(nodalflag);

    // Create temp arrays and copy existing data into these temporary arrays -- srcMF and destMF
    //        have the same BoxArray and DistributionMapping here
    for (int i = 0; i < 3; ++i) {
        int ng = current[lev][i]->nGrow();
	old_current[i].reset(new MultiFab(grids[lev],dmap[lev],1,ng,info));
	old_Efield [i].reset(new MultiFab(grids[lev],dmap[lev],1,ng,info));
	old_Bfield [i].reset(new MultiFab(grids[lev],dmap[lev],1,ng,info));
        MultiFab::Copy(*old_current[i],*current[lev][i],0,0,current[lev][i]->nComp(),0);
        MultiFab::Copy( *old_Efield[i], *Efield[lev][i],0,0, Efield[lev][i]->nComp(),0);
        MultiFab::Copy( *old_Bfield[i], *Bfield[lev][i],0,0, Bfield[lev][i]->nComp(),0);
    }

    // This creates a new BoxArray and DistributionMapping and assigns the particles to that 
    // Here we need to re-define the grids for the mesh quantities as well
    bool remapped = LoadBalanceBaseLevel();

    // Copy "old" temp data into the new arrays -- here the src and dest do NOT
    //       have the same BoxArray and DistributionMapping -- only need to do this if 
    //       the grids have actually changed
    if (remapped) {
       for (int i = 0; i < 3; ++i) {
        current[lev][i]->copy(*old_current[i], 0, 0, current[lev][i]->nComp()); 
         Bfield[lev][i]->copy( *old_Bfield[i], 0, 0,  Bfield[lev][i]->nComp()); 
         Efield[lev][i]->copy( *old_Efield[i], 0, 0,  Bfield[lev][i]->nComp()); 
       }
    }
}

bool
WarpX::okToRegrid(int step)
{
    if (regrid_int < 0)
        return false;
    else
        return (step%regrid_int == 0);
}

bool
WarpX::LoadBalanceBaseLevel()
{
    int min_grid_size = 4;

    bool remapped;

    // **************************************************************************
    // Load Balance
    // **************************************************************************
    const Real eff_target = 0.8;
    const int lev = 0;

    Array<long> new_particle_cost = mypc->NumberOfParticlesInGrid(lev);
    Real neweff = getEfficiency(dmap[0], new_particle_cost);

    if (debug_lb >= 1 && ParallelDescriptor::IOProcessor()) 
    {
       long min_cost = new_particle_cost[0];
       long max_cost = new_particle_cost[0];
       for (int i = 1; i < new_particle_cost.size(); i++) 
       {
          min_cost = std::min(new_particle_cost[i],min_cost);
          max_cost = std::max(new_particle_cost[i],max_cost);
       }

         std::cout << "ORIG MIN COST / MAX COST / EFF " << min_cost << 
                                                    " " << max_cost << " " << neweff << std::endl;
    }
 
    WarpXParticleContainer* myspc = &(mypc->GetParticleContainer(0));

    // This is what we are starting with
    BoxArray new_ba= myspc->ParticleBoxArray(0);
    DistributionMapping new_dm = myspc->ParticleDistributionMap(0);
    if (debug_lb >= 1 && ParallelDescriptor::IOProcessor()) 
       std::cout << "OLD BA HAS " << new_ba.size() << " GRIDS" << std::endl;

    if (neweff < eff_target) 
    {
	    Real oldeff;
	    Array<long> old_particle_cost;
	    int heavy_grid_size = this->maxGridSize(0);

	    do {
		oldeff = neweff;
		old_particle_cost = new_particle_cost;

		if (ParallelDescriptor::IOProcessor()) 
		{
		    std::cout << "*** " << std::endl;
		    std::cout << "*** Before remapping, # of boxes: " << new_particle_cost.size()
			      << ", efficiency: " << neweff << "\n";
		}

		new_ba = myspc->ParticleBoxArray(lev);
		// This returns new_particle_cost as an *estimate* of the new cost per grid, based just
		//      on dividing the cost proportionally as the grid is divided
		splitBoxes(new_ba, new_particle_cost, old_particle_cost, heavy_grid_size);
		heavy_grid_size /= 2;

		// We use this approximate cost to get a new DistrbutionMapping so we can go ahead 
		//      and move the particles
		new_dm = getCostCountDM(new_particle_cost, new_ba);

		// We get an *estimate* of the new efficiency
		neweff = getEfficiency(new_dm, new_particle_cost);

		if (ParallelDescriptor::IOProcessor()) 
		{
		    std::cout << "*** If     remapping, # of boxes: " << new_particle_cost.size()
			      << ", approx. eff: " <<  neweff << "\n";
		}

		// Only if the new_ba and new_dm are expected to improve the efficiency, ...
		if (neweff > oldeff)  
		{
		    // Now we actually move the particles onto the new_ba with the new_dm
		    mypc->SetParticleBoxArray(lev,new_ba);
		    mypc->SetParticleDistributionMap(lev,new_dm);
		    mypc->Redistribute();
		    
		    // This counts how many particles are *actually* in each grid of the 
		    //      ParticleContainer's new ParticleBoxArray
		    new_particle_cost = mypc->NumberOfParticlesInGrid(lev);
		    
		    // Here we get the *actual* new efficiency
		    neweff = getEfficiency(new_dm, new_particle_cost);
		    
		    if (ParallelDescriptor::IOProcessor()) 
		    {
			std::cout << "*** After  remapping, # of boxes: " << new_particle_cost.size()
				  << ", actual  eff: " <<  neweff << "\n";
		    }
		}
	    } 
	    while (neweff < eff_target && neweff > oldeff && heavy_grid_size >= 2*min_grid_size);

        if (debug_lb && ParallelDescriptor::IOProcessor()) 
        {
           BoxArray new_ba = myspc->ParticleBoxArray(lev);
           std::cout << "NEW BA HAS " << new_ba.size() << " GRIDS" << std::endl;

           long min_cost = new_particle_cost[0];
           long max_cost = new_particle_cost[0];
           for (int i = 1; i < new_particle_cost.size(); i++) 
           {
              min_cost = std::min(new_particle_cost[i],min_cost);
              max_cost = std::max(new_particle_cost[i],max_cost);
           }

           std::cout << "NEW  MIN COST / MAX COST / EFF " << min_cost << 
                 " " << max_cost << " " << neweff << std::endl;
        }

        AllocLevelData(0,new_ba,new_dm);
        SetBoxArray(0, new_ba);
        SetDistributionMap(0, new_dm);

        remapped = true;

    } else {

        if (debug_lb >= 1 && ParallelDescriptor::IOProcessor()) 
	{
	    std::cout << "*** " << std::endl;
	    std::cout << "*** No remapping required: # of boxes: " << grids[0].size() 
			  << ", efficiency: " <<  neweff << "\n";
	    std::cout << "*** " << std::endl;
        }
        remapped = false;
    }
    return remapped;
}

Real
WarpX::getEfficiency(const DistributionMapping& dm, const Array<long>& cost)
{
    Array<long> cpr(ParallelDescriptor::NProcs(),0);
    Real ctot=0;
    for (int i=0, N=cost.size(); i<N; i++) {
        ctot += cost[i];
        cpr[dm[i]] += cost[i];
    }
    long cmax = *std::max_element(cpr.begin(), cpr.end());
    Real cavg = ctot / ParallelDescriptor::NProcs();
    return cavg / cmax;
}

DistributionMapping
WarpX::getCostCountDM (const Array<long>& cost, const BoxArray& ba)
{
    DistributionMapping res;
    int nprocs = ParallelDescriptor::NProcs();
    const Real factor = 1.5; // A process can get up to 'factor' times of the average number of boxes.
    int nmax = (cost.size()+nprocs-1) / nprocs * factor;
    Real eff;
    res.KnapSackProcessorMap(cost, nprocs, &eff, true, nmax);
    return res;
}

void
WarpX::splitBoxes (BoxArray& ba, Array<long>& newcost, const Array<long>& cost_in, int heavy_grid_size)
{
    long totcost = 0;
    for (int i = 0; i < cost_in.size(); i++)
       totcost += cost_in[i];

    long avgcost    = totcost / ParallelDescriptor::NProcs();
    long cost_split = avgcost / 2;
    long cost_merge = avgcost / 3;

    newcost = cost_in;

    int base_box_size = heavy_grid_size;
    int half_box_size = base_box_size/2;
    if (half_box_size*2 != base_box_size) return;

    Array<long> cost = newcost;

    BoxList newbl;
    newcost.clear();
	
    for (int i=0, N=ba.size(); i<N; i++) 
    {
	    const Box& bx = ba[i];
	    const Real ct = cost[i];
	    if (ct > cost_split &&
		D_TERM(bx.length(0) == base_box_size,
		    && bx.length(1) == base_box_size,
		    && bx.length(2) == base_box_size))
	    { // split it
		BoxList bltmp(bx);
		bltmp.maxSize(half_box_size);
		BL_ASSERT(bltmp.size() == D_TERM(2,*2,*2));
		Real hct = ct / bltmp.size();
		for (int i=0; i<bltmp.size(); i++)
		    newcost.push_back(hct);
		newbl.catenate(bltmp);
	    } 
	    else if (D_TERM(bx.length(0) == half_box_size,
			 && bx.length(1) == half_box_size,
		         && bx.length(2) == half_box_size)) 
	    { // maybe we can merge
		int nm = D_TERM(2,*2,*2);
		if (i+nm-1 >= N) {
		    // not enough boxes to merge
		    newbl.push_back(bx);
		    newcost.push_back(ct);
		    continue;
		} else {
		    // they must have same size and they must merge into a single box of base_box_size
		    bool samesize = true;
		    Box mergedbox(bx);
		    Real mergedcost = ct;
		    for (int j=i+1; j<i+nm; j++) {
			if (! bx.sameSize(ba[j])) {
			    samesize = false;
			    break;
			} else {
			    mergedbox.minBox(ba[j]);
			    mergedcost += cost[j];
			}
		    }
		
		    if (samesize && 
			D_TERM(mergedbox.length(0) == base_box_size,
			    && mergedbox.length(1) == base_box_size,
			    && mergedbox.length(2) == base_box_size)) 
		    {
			if (mergedcost < cost_merge) {
			    newbl.push_back(mergedbox);
			    newcost.push_back(mergedcost);
			} else {
			    for (int j=i; j<i+nm; j++) {
				newbl.push_back(ba[j]);
				newcost.push_back(cost[j]);
			    }
			}
			i += nm-1;  // skip some boxes becuase they have been processed
		    } else {
			newbl.push_back(bx);
			newcost.push_back(ct);
		    }
		}
	    }
	    else 
	    {
		newbl.push_back(bx);
		newcost.push_back(ct);
	    }
    }

    ba = BoxArray(newbl);

    BL_ASSERT(ba.size() == newcost.size());
}
