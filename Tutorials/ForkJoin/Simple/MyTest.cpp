#include "MyTest.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ForkJoin.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();
    initData();
}

static void
myFunction(ForkJoin& fj)
{
    const MultiFab& myData_split = fj.get_mf("data_split");
    const MultiFab& myData_single = fj.get_mf("data_single");
    const MultiFab& myData_all = fj.get_mf("data_all");

    int taskID = fj.MyTask();

    const ForkJoin::ComponentSet bounds_split = fj.ComponentBounds("data_split");
    const ForkJoin::ComponentSet bounds_single = fj.ComponentBounds("data_single");
    const ForkJoin::ComponentSet bounds_all = fj.ComponentBounds("data_all");

    Print() << "Task ID: " << taskID << std::endl;
    Print() << "  Number of boxes:  " << myData_split.boxArray().size() << std::endl;
    Print() << "  Number of ranks:  " << ParallelDescriptor::NProcs() << std::endl;
    if (myData_split.nComp() > 0) {
        Print() << "  Number of split comps:      " << myData_split.nComp() << std::endl;
        Print() << "  Split component bounds:     " << bounds_split.lo << ":" << bounds_split.hi -1 << std::endl;
        Print() << "  Number of split grow cells: " << myData_split.nGrow() << std::endl;
    } else {
        Print() << "  No split data on this task " << std::endl;
    }

    if (myData_single.nComp() > 0) {
        Print() << "  Number of single comps:      " << myData_single.nComp() << std::endl;
        Print() << "  Single component bounds:     " << bounds_single.lo << ":" << bounds_single.hi -1 << std::endl;
        Print() << "  Number of single grow cells: " << myData_single.nGrow() << std::endl;
    } else {
        Print() << "  No single data on this task " << std::endl;
    }

    if (myData_all.nComp() > 0) {
        Print() << "  Number of all comps:      " << myData_all.nComp() << std::endl;
        Print() << "  All component bounds:     " << bounds_all.lo << ":" << bounds_all.hi -1 << std::endl;
        Print() << "  Number of all grow cells: " << myData_all.nGrow() << std::endl;
    } else {
        Print() << "  No all data on this task " << std::endl;
    }
}

void
MyTest::runTest ()
{
    BL_ASSERT(n_tasks > 0);
    ForkJoin fj(n_tasks);
    fj.SetVerbose(fj_verbose);
    fj.set_task_output_dir(outdir);

    fj.reg_mf(data_split,"data_split",ForkJoin::Strategy::split,ForkJoin::Intent::in);
    fj.reg_mf(data_single,"data_single",ForkJoin::Strategy::single,ForkJoin::Intent::in,n_tasks-1);
    fj.reg_mf(data_all,"data_all",ForkJoin::Strategy::duplicate,ForkJoin::Intent::in);

    if (data_single.nGrow() > 0) {
        fj.modify_ngrow("data_single", IntVect(data_single.nGrow() - 1)); // reduce grow cells by one
    }

    fj.fork_join(myFunction);
}

void
MyTest::readParameters ()
{
    ParmParse pp;

    int n_cell;
    pp.get("n_cell",n_cell);

    int max_grid_size;
    pp.get("max_grid_size",max_grid_size);

    int n_comp_split;
    pp.get("n_comp_split",n_comp_split);

    int n_comp_single;
    pp.get("n_comp_single",n_comp_single);

    int n_comp_all;
    pp.get("n_comp_all",n_comp_all);

    Box domain(IntVect(D_DECL(0,0,0)),IntVect(n_cell-1,n_cell-1,n_cell-1));
    BoxArray grids(domain);
    grids.maxSize(max_grid_size);

    pp.get("n_tasks",n_tasks);

    int n_grow = 1;
    DistributionMapping dmap(grids);
    data_split.define(grids,dmap,n_comp_split,n_grow);
    data_single.define(grids,dmap,n_comp_single,n_grow+1);
    data_all.define(grids,dmap,n_comp_all,n_grow+1);

    pp.query("outdir",outdir);
    fj_verbose = false;
    pp.query("v",fj_verbose);
}

void
MyTest::initData ()
{
    for (int n=0; n<data_split.nComp(); ++n) {
        data_split.setVal(n+1,n,1);
    }
    for (int n=0; n<data_single.nComp(); ++n) {
        data_single.setVal(10*(n+1),n,1);
    }
    for (int n=0; n<data_all.nComp(); ++n) {
        data_all.setVal(100*(n+1),n,1);
    }
}

