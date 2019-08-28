#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleLocator.H>

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
};

void testParticleMesh(TestParams& parms)
{

  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
  IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
  const Box domain(domain_lo, domain_hi);

  // This sets the boundary conditions to be doubly or triply periodic
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) 
    is_per[i] = 1; 
  Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);
  
  BoxArray ba(domain);
  ba.maxSize(parms.max_grid_size);

  ParticleLocator ploc;
  ploc.build(ba);

  auto assign_grid = ploc.getGridAssignor();

  for (int i = 0; i < ba.size(); ++i) 
  {
      const Box& box = ba[i];

      Gpu::HostVector<IntVect> host_cells;
      for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) host_cells.push_back(iv);

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

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms;
  
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  
  testParticleMesh(parms);
  
  amrex::Finalize();
}
