#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <unistd.h>
#include <cstdio>

using namespace amrex;

struct InputParams {
  int ncells;
  int max_grid_size;
  int ncomp;
  int nlevs;

  //int restart_check = 0;
  int nplotfile = 1;
  //int sleeptime = 0;
};


void loadParameters(InputParams& input)
{
    ParmParse pp;
    pp.get("ncells", input.ncells);
    pp.get("max_grid_size", input.max_grid_size);
    pp.get("ncomp", input.ncomp);
    pp.get("nlevs", input.nlevs);

    pp.query("nplotfile", input.nplotfile);
    //pp.query("sleeptime", input.sleeptime);
    //pp.query("restart_check", input.restart_check);
    //pp.query("grids_from_file", input.grids_from_file);
}

void set_grids_nested (const InputParams& input,
                       Vector<Box>& domains,
                       Vector<BoxArray>& grids,
                       Vector<IntVect>& ref_ratio);

void test ();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    test();
    amrex::Finalize();
    std::cout<<"Finalized "<<std::endl;
}

struct TestField
{
  TestField(const InputParams& inputs)
  {
    const int nghost = 0;

    Vector<Box> domains;
    Vector<BoxArray> ba;
    set_grids_nested(inputs, domains, ba, m_Ref_ratio);

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM] = {AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object for each level
    m_Geom.resize(inputs.nlevs);

    m_Geom[0].define(domains[0], &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < inputs.nlevs; lev++) {
        m_Geom[lev].define(domains[lev], &real_box, CoordSys::cartesian, is_per);
    }

    Vector<DistributionMapping> dmap(inputs.nlevs);

    m_mf.resize(inputs.nlevs);
    for (int lev = 0; lev < inputs.nlevs; lev++) {
        std::cout<<ba[lev]<<std::endl;
        std::cout<<dmap[lev]<<std::endl;
        dmap[lev] = DistributionMapping{ba[lev]};
        m_mf[lev] = std::make_unique<MultiFab>(ba[lev], dmap[lev], inputs.ncomp, nghost);
        m_mf[lev]->setVal(lev);
    }

    for (int i = 0; i < inputs.ncomp; ++i) {
      m_Varnames.push_back("component_" + std::to_string(i));
    }

    //Vector<int> level_steps(inputs.nlevs, 0);
    m_Level_steps.assign(inputs.nlevs, 0);
  }

  Vector<int> m_Level_steps;

  Vector<std::string> m_Varnames;
  Real m_Time = Real(0.0);
  Vector<IntVect> m_Ref_ratio;
  Vector<Geometry> m_Geom;
  Vector<std::unique_ptr<MultiFab> > m_mf;
};


void saveFile(char const* fname, const InputParams& inputs, const TestField& testField)
{
#ifdef AMREX_USE_OPENPMD_API
  openpmd_api::InitHandler(fname);

  for (int ts = 0; ts < inputs.nplotfile; ts++)
    {
      openpmd_api::SetStep(ts);
      openpmd_api::WriteMultiLevel(//fname,
                                   amrex::GetVecOfConstPtrs(testField.m_mf), testField.m_Varnames,
                                   testField.m_Geom, testField.m_Time, /*testField.m_Level_steps,*/ testField.m_Ref_ratio);
      openpmd_api::CloseStep(ts);

      //saveFile(fname, ts, inputs, testField);
      amrex::Print()<<"Timestep: "<<ts<<" done \n";
    }

  openpmd_api::CloseHandler();

#else
  /*WriteMultiLevelPlotfile(fname, inputs.nlevs, amrex::GetVecOfConstPtrs(testField.m_mf),
                                testField.m_Varnames, testField.m_Geom,
                                testField.m_Time, testField.m_Level_steps, testField.m_Ref_ratio);
  */

#endif
}

void test ()
{
    InputParams inputs;
    loadParameters(inputs);
    TestField testField(inputs);


    //char fname[512]="adios_diag";
    char fname[512]="";
    saveFile(fname, inputs, testField);

    /* ParallelDescriptor::Barrier(); */
}


void set_grids_nested (const InputParams& input,
                       Vector<Box>& domains,
                       Vector<BoxArray>& grids,
                       Vector<IntVect>& ref_ratio)
{
  //int ncells, max_grid_size, nlevs;

    AMREX_ALWAYS_ASSERT(input.nlevs < 2); // relax this later

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(input.ncells-1, input.ncells-1, input.ncells-1));

    domains.resize(input.nlevs);
    domains[0].setSmall(domain_lo);
    domains[0].setBig(domain_hi);

    ref_ratio.resize(input.nlevs-1);
    for (int lev = 1; lev < input.nlevs; lev++) {
        ref_ratio[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    grids.resize(input.nlevs);
    grids[0].define(domains[0]);

    // Now we make the refined level be the center eighth of the domain
    if (input.nlevs > 1) {
        int n_fine = input.ncells * ref_ratio[0][0];
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        grids[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < input.nlevs; lev++) {
        grids[lev].maxSize(input.max_grid_size);
    }

    for (int lev = 1; lev < input.nlevs; lev++) {
        domains[lev] = amrex::refine(domains[lev-1], ref_ratio[lev-1]);
    }
}
