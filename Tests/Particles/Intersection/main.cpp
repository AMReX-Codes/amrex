#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleLocator.H>

using namespace amrex;

struct TestParams {
    IntVect size;
    int max_grid_size;
    int is_periodic;
    int nlevs;
};

void get_test_params(TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("is_periodic", params.is_periodic);
    pp.get("nlevs", params.nlevs);
}

void testIntersection()
{
    TestParams params;
    get_test_params(params, "intersect");
    
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = params.is_periodic;
    
    Vector<IntVect> rr(params.nlevs-1);
    for (int lev = 1; lev < params.nlevs; lev++)
        rr[lev-1] = IntVect(D_DECL(2,2,2));
    
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    IntVect domain_lo(D_DECL(0 , 0, 0));
    IntVect domain_hi(D_DECL(params.size[0]-1, params.size[1]-1, params.size[2]-1));
    const Box base_domain(domain_lo, domain_hi);
    
    Vector<Geometry> geom(params.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < params.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }
    
    Vector<BoxArray> ba(params.nlevs);
    IntVect lo = IntVect(D_DECL(0, 0, 0));
    IntVect size = params.size;
    for (int lev = 0; lev < params.nlevs; ++lev)
    {
        ba[lev].define(Box(lo, lo+params.size-1));
        ba[lev].maxSize(params.max_grid_size);
        lo += size/2;
        size *= 2;
    }

    Vector<ParticleLocator> ploc(params.nlevs);

    for (int lev = 0; lev < params.nlevs; ++lev)
    {
        ploc[lev].build(ba[lev]);

        auto assign_grid = ploc[lev].getGridAssignor();

        for (int i = 0; i < ba[lev].size(); ++i) 
        {
            const Box& box = ba[lev][i];

            Gpu::HostVector<IntVect> host_cells;
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) host_cells.push_back(iv);
            //host_cells.push_back(box.smallEnd());
        
            int num_cells = host_cells.size();
      
            Gpu::DeviceVector<IntVect> device_cells(num_cells);
            Gpu::thrust_copy(host_cells.begin(), host_cells.end(), device_cells.begin());
            
            Gpu::DeviceVector<int> device_grids(num_cells);
            
            const auto cells_ptr = device_cells.dataPtr();
            const auto grids_ptr = device_grids.dataPtr();
            amrex::ParallelFor(num_cells, [=] AMREX_GPU_DEVICE (int j) noexcept
            {
                grids_ptr[j] = assign_grid(cells_ptr[j]);
            });

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            
            reduce_op.eval(num_cells, reduce_data,
            [=] AMREX_GPU_DEVICE (int j) -> ReduceTuple
            {
                return {grids_ptr[j] != i};
            }); 
      
            ReduceTuple hv = reduce_data.value();
            
            int num_wrong = amrex::get<0>(hv);
            AMREX_ALWAYS_ASSERT(num_wrong == 0);
        }
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    testIntersection();
    
    amrex::Finalize();
}
