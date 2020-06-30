
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "CheckPair.H"

#include "MDParticleContainer.H"

using namespace amrex;

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int nsteps;
    int num_rebuild;
    int num_ppc;
    bool print_min_dist;
    bool print_neighbor_list;
    bool print_num_particles;
    bool write_particles;
    Real cfl;
};

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}

void get_test_params(TestParams& params)
{
    ParmParse pp;
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("nsteps", params.nsteps);
    pp.get("print_minimum_distance", params.print_min_dist);
    pp.get("print_neighbor_list", params.print_neighbor_list);
    pp.get("write_particles", params.write_particles);
    pp.get("num_rebuild", params.num_rebuild);
    pp.get("num_ppc", params.num_ppc);
    pp.get("cfl", params.cfl);
    pp.get("print_num_particles", params.print_num_particles);
}

void main_main ()
{

    amrex::Print() << "Running MD benchmark \n";

    TestParams params;
    get_test_params(params);

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 0;
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    const int ncells = 1;
    MDParticleContainer pc(geom, dm, ba, ncells);

    int npc = params.num_ppc;
    IntVect nppc = IntVect(AMREX_D_DECL(npc, npc, npc));

    pc.InitParticles(nppc, 1.0, 0.0);

    if (params.print_num_particles) 
      amrex::Print() << "Num particles after init is " << pc.TotalNumberOfParticles() << "\n";

    int num_rebuild = params.num_rebuild;

    Real cfl = params.cfl;
    
    Real min_d = std::numeric_limits<Real>::max();

    for (int step = 0; step < params.nsteps; ++step) {

	Real dt = pc.computeStepSize(cfl);

	if (step % num_rebuild == 0)
	{
	  if (step > 0) pc.RedistributeLocal();

	  pc.fillNeighbors();

	  pc.buildNeighborList(CheckPair());
	} 
	else
	{
	  pc.updateNeighbors();
	}

        if (params.print_min_dist) 
	   min_d = std::min(min_d, pc.minDistance());

        if (params.print_neighbor_list) 
           pc.printNeighborList();

	pc.computeForces();

	pc.moveParticles(dt);
    }

    pc.RedistributeLocal();

    if (params.print_min_dist     ) amrex::Print() << "Min distance  is " << min_d << "\n";
    if (params.print_num_particles) amrex::Print() << "Num particles is " << pc.TotalNumberOfParticles() << "\n";

    if (params.write_particles) pc.writeParticles(params.nsteps);
}
