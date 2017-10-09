#include <iostream>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BoxArray.H>

using namespace amrex;

void
splitBoxes (BoxArray& ba, Vector<long>& newcost, const Vector<long>& cost_in, int heavy_grid_size)
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

    Vector<long> cost = newcost;

    BoxList newbl;
    newcost.clear();
	
    for (int i=0, N=ba.size(); i<N; i++) 
    {
	    const Box& bx = ba[i];
	    const Real ct = cost[i];
	    if (ct > cost_split &&
		AMREX_D_TERM(bx.length(0) == base_box_size,
		    && bx.length(1) == base_box_size,
		    && bx.length(2) == base_box_size))
	    { // split it
		BoxList bltmp(bx);
		bltmp.maxSize(half_box_size);
		BL_ASSERT(bltmp.size() == AMREX_D_TERM(2,*2,*2));
		Real hct = ct / bltmp.size();
		for (int i=0; i<bltmp.size(); i++)
		    newcost.push_back(hct);
		newbl.catenate(bltmp);
	    } 
	    else if (AMREX_D_TERM(bx.length(0) == half_box_size,
			 && bx.length(1) == half_box_size,
		         && bx.length(2) == half_box_size)) 
	    { // maybe we can merge
		int nm = AMREX_D_TERM(2,*2,*2);
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
			AMREX_D_TERM(mergedbox.length(0) == base_box_size,
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
