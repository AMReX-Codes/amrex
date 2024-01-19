#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>
#include <AMReX_BoxArray.H>
#include <MappedPC.H>
#include <TerrainPC.H>

using namespace amrex;

enum struct GridType {
    Terrain, Mapped
};

enum struct ProbType {
    Annulus, Stretched, Hill
};

enum struct AdvectionType {
    mac, cc, nd
};

struct TestParams
{
    // Number of grid cells in each direction
    IntVect size;

    // Periodicity of each coordinate direction
    IntVect is_periodic;

    // Maximum size of any grid in each coordinate direction
    int max_grid_size;

    // Number of particles to be initialized in each cell
    int num_ppc;

    // How many time steps to run
    int nsteps;

    // How many levels of refinement to run (default = 1, aka "single level")
    int nlevs;

    // Types of meshes currently include
    //    "terrain" -- the grid is regularly spaced in the first (AMREX_SPACEDIM-1) directions,
    //                 but the last direction is irregularly spaced and the height at nodes is
    //                 stored in the a_z_loc array
    //    "mapped"  -- the grid is irregularly spaced in all directions and the locations at nodes
    //                 are stored in the a_xyz_loc array
    GridType grid_type = GridType::Terrain;

    // Types of advection currently include
    //    "cc"  -- all components of the velocity are stored at cell centers
    //    "nd"  -- all components of the velocity are stored at cell corners (nodes)
    //    "mac" -- normal components of velocity are stored on the cell faces
    std::string advection_type = "mac";

    // Combinations that are currently allowed and can be tested here
    // "terrain" + {cc, nd, mac}
    // "mapped"  + {nd}

    // Sample mesh initializations currently include
    // "annulus"   -- this is of grid_type GridType::Mapped and is a periodic rectangle mapped into an annulus
    // "stretched" -- this is of grid_type GridType::Terrain in which the grid spacing dz only depends on z
    // "hill"      -- this is of grid_type GridType::Terrain in which the grid spacing dz varies with x and y as well as z
    ProbType prob_type = ProbType::Hill;

    // We have separate calls for cc/nd vs mac.  The location of the velocity in the cc_nd call
    //    is determined at run-time by checking the boxArray type of the velocity MultiFab

    // Currently the domain size is hard-wired to [0:1, 0.5] if 2D and [0:1, 0:0.5, 0:0.25] if 3D
};

void Test ();

void InitAnnulus (MultiFab& a_xyz_loc, Geometry& geom);
void InitStretched (MultiFab& a_z_loc  , Geometry& geom);
void InitHill (MultiFab& a_z_loc  , Geometry& geom);

void InitUmac  (MultiFab* umac, const MultiFab& a_xyz_loc, Geometry& geom, ProbType prob_type);
void InitUCC (MultiFab& u   , const MultiFab& a_xyz_loc, Geometry& geom, ProbType& prob_type);
void InitUND (MultiFab& u   , const MultiFab& a_xyz_loc, Geometry& geom, ProbType& prob_type);

void get_test_params(TestParams& params)
{
    ParmParse pp("");
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("num_ppc", params.num_ppc);
    pp.get("is_periodic", params.is_periodic);
    pp.get("nsteps", params.nsteps);
    pp.get("nlevs", params.nlevs);

    pp.get("advection_type", params.advection_type);
    AMREX_ALWAYS_ASSERT(params.advection_type == "mac" || params.advection_type == "cc" || params.advection_type == "nd");

    std::string grid_type_string;
    pp.get("grid_type"     , grid_type_string);
    AMREX_ALWAYS_ASSERT(grid_type_string == "mapped" || grid_type_string == "terrain");
    params.grid_type = (grid_type_string == "mapped") ? GridType::Mapped : GridType::Terrain;

    std::string prob_type_string;
    pp.get("prob_type", prob_type_string);
    AMREX_ALWAYS_ASSERT(prob_type_string == "annulus" || prob_type_string == "stretched" || prob_type_string == "hill");
    if (prob_type_string == "annulus") params.prob_type = ProbType::Annulus;
    if (prob_type_string == "stretched") params.prob_type = ProbType::Stretched;
    if (prob_type_string == "hill") params.prob_type = ProbType::Hill;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    Test();

    amrex::Finalize();
}

void Test()
{
    BL_PROFILE("Test");
    TestParams params;
    get_test_params(params);

    int is_per[] = {AMREX_D_DECL(params.is_periodic[0],
                                 params.is_periodic[1],
                                 params.is_periodic[2])};

    RealBox real_box;
    real_box.setLo(0, 0.0);
    real_box.setHi(0, 1.0);

    real_box.setLo(1, 0.0);
    real_box.setHi(1, 0.5);

#if (AMREX_SPACEDIM == 3)
    real_box.setLo(2, 0.0);
    real_box.setHi(2, 0.5);
#endif

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box base_domain(domain_lo, domain_hi);

    Vector<Geometry> geom(params.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);

    Vector<BoxArray> ba(params.nlevs);
    Vector<DistributionMapping> dm(params.nlevs);
    IntVect lo(0);
    IntVect size = params.size;
    for (int lev = 0; lev < params.nlevs; ++lev)
    {
        ba[lev].define(Box(domain_lo, domain_hi));
        ba[lev].maxSize(params.max_grid_size);
        dm[lev].define(ba[lev]);
        lo += size/2;
        size *= 2;
    }

    // We currently assume a single-level problem
    int lev = 0;

    auto dx = geom[lev].CellSize();

    // We define both types of particles here, but in separate particle containers
    MappedPC   mapped_pc(geom[lev], dm[lev], ba[lev]);
    TerrainPC terrain_pc(geom[lev], dm[lev], ba[lev]);

    IntVect nppc(params.num_ppc);

    // **************************************************************************************
    // Define xyz_phys on nodes
    // **************************************************************************************
    BoxArray ba_nd(ba[lev]); ba_nd.surroundingNodes();

    // This has one component, the "height"
    MultiFab   a_z_loc(ba_nd,dm[lev],1,1);

    // This has AMREX_SPACEDIM components to define (x,y,z) as a function of (i,j,k)
    MultiFab a_xyz_loc(ba_nd,dm[lev],AMREX_SPACEDIM,1);

    if (params.grid_type == GridType::Mapped and params.prob_type == ProbType::Annulus) {
        InitAnnulus(a_xyz_loc, geom[lev]);
    } else if (params.grid_type == GridType::Terrain and params.prob_type == ProbType::Stretched) {
        InitStretched(a_z_loc, geom[lev]);
    } else if (params.grid_type == GridType::Terrain and params.prob_type == ProbType::Hill) {
        InitHill(a_z_loc, geom[lev]);
    } else {
        amrex::Abort("This grid_type and prob_type combination not allowed.");
    }

    // **************************************************************************************
    // Define velocity components (for now it is constant velocity in "i" direction)
    // **************************************************************************************
    MultiFab ucc;                  // Cell-centered
    MultiFab und;                  // Node-centered
    MultiFab umac[AMREX_SPACEDIM]; // Face-centered

    if (params.grid_type == GridType::Terrain && params.advection_type == "mac")
    {
        BoxArray ba_x(ba[lev]); ba_x.convert(IntVect(AMREX_D_DECL(1,0,0)));
        umac[0].define(ba_x,dm[lev],1,1); umac[0].setVal(0.0);
        umac[0].FillBoundary(geom[lev].periodicity());

        BoxArray ba_y(ba[lev]); ba_y.convert(IntVect(AMREX_D_DECL(0,1,0)));
        umac[1].define(ba_y,dm[lev],1,1); umac[1].setVal(0.0);
        umac[1].FillBoundary(geom[lev].periodicity());

#if (AMREX_SPACEDIM == 3)
        BoxArray ba_z(ba[lev]); ba_z.convert(IntVect(AMREX_D_DECL(0,0,1)));
        umac[2].define(ba_z,dm[lev],1,1); umac[2].setVal(0.0);
        umac[2].FillBoundary(geom[lev].periodicity());
#endif

        InitUmac(&umac[0], a_z_loc, geom[lev], params.prob_type);

    } else if (params.grid_type == GridType::Terrain && params.advection_type == "cc") {
        ucc.define(ba[lev],dm[lev],AMREX_SPACEDIM,1);
        InitUCC (ucc, a_z_loc, geom[lev], params.prob_type);
        ucc.FillBoundary(geom[lev].periodicity());

    } else if (params.grid_type == GridType::Terrain && params.advection_type == "nd") {
        und.define(ba_nd,dm[lev],AMREX_SPACEDIM,1);
        InitUND (und, a_z_loc, geom[lev], params.prob_type);
        und.FillBoundary(geom[lev].periodicity());

    } else if (params.grid_type == GridType::Mapped && params.advection_type == "nd") {
        und.define(ba_nd,dm[lev],AMREX_SPACEDIM,1);
        InitUND (und, a_xyz_loc, geom[lev], params.prob_type);
        und.FillBoundary(geom[lev].periodicity());

    } else {
        amrex::Abort("Unknown grid_type and advection_type combination.");
    }

    // **************************************************************************************
    // Initialize and write out particle locations
    // **************************************************************************************
    if (params.grid_type == GridType::Mapped)
    {
        mapped_pc.InitParticles(a_xyz_loc);
        mapped_pc.WritePlotFile("plot", "particles");
    }
    else if (params.grid_type == GridType::Terrain)
    {
        terrain_pc.InitParticles(a_z_loc);
        terrain_pc.WritePlotFile("plot", "particles");
    }

    // **************************************************************************************
    // Advance the particle positions in time based on the face-based velocity
    // **************************************************************************************
    std::string plotfilename;

    amrex::Real max_vel;
    if (params.advection_type == "mac") {
        max_vel = umac[0].max(0,0,false);
    } else if (params.advection_type == "ucc") {
        max_vel = ucc.max(0,0,false);
    } else if (params.advection_type == "und") {
        max_vel = und.max(0,0,false);
    }

    // This is assuming velocity only in x-direction
    amrex::Real dt = 0.9 * dx[0] / max_vel;
    amrex::Print() << "COMPUTING DT TO BE " << dt << " BASED ON MAX VEL " << max_vel << std::endl;

    for (int nt = 0; nt < params.nsteps; nt++)
    {
        if (params.grid_type == GridType::Terrain && params.advection_type == "mac") {
            amrex::Print() << "Advecting at time " << nt << " using MAC velocities" << std::endl;
            terrain_pc.AdvectWithUmac(&umac[0], 0, dt, a_z_loc);

        } else if (params.grid_type == GridType::Terrain && params.advection_type == "cc") {
            amrex::Print() << "Advecting at time " << nt << " using cell-centered velocities" << std::endl;
            terrain_pc.AdvectWithUCC(ucc, 0, dt, a_z_loc);

        } else if (params.grid_type == GridType::Terrain && params.advection_type == "nd") {
            amrex::Print() << "Advecting at time " << nt << " using node-centered velocities" << std::endl;
            terrain_pc.AdvectWithUND(und, 0, dt, a_z_loc);

        } else if (params.grid_type == GridType::Mapped && params.advection_type == "nd") {
            ucc.FillBoundary(geom[0].periodicity());
            amrex::Print() << "Advecting at time " << nt << " using node-centered velocities" << std::endl;
            mapped_pc.AdvectWithUND(und, 0, dt, a_xyz_loc);
        }

        plotfilename = Concatenate("plot", nt, 5);
        if (params.grid_type == GridType::Terrain) {
            terrain_pc.WritePlotFile(plotfilename, "particles");
        } else if (params.grid_type == GridType::Mapped) {
            mapped_pc.WritePlotFile(plotfilename, "particles");
        }
    }
}
