#include "hgparallel.h"

task_copy::task_copy(MultiFab& mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
: m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx)
#ifdef BL_USE_MPI
, s_tmp(0), d_tmp(0)
#endif
{
#ifdef BL_USE_MPI
    if ( m_mf.DistributionMap()[m_sgrid] == ParallelDescriptor::MyProc() )
    {
	d_tmp = new FArrayBox(bx, m_mf.nComp());
	MPI_Irecv(d_tmp->dataPtr());
    }
    if ( m_mf.DistributionMap()[m_dgrid] == ParallelDescriptor::MyProc() ) 
    {
	s_tmp = new FArrayBox(bx, m_mf.nComp());
	s_tmp->copy(m_smf[m_sgrid]);
	MPI_Isend(s_tmp->dataPtr(), s_tmp->box().numPts() * s_tmp->nComp(), MPI_DOUBLE );
    }
#else
    m_ready = true;
    m_mf[m_dgrid].copy(m_smf[m_sgrid], m_bx);
#endif
}

task_copy::task_copy(MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
: m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), s_bx(sb), m_sgrid(sgrid)
{
    m_ready = true;
    m_mf[m_dgrid].copy(m_smf[m_sgrid], s_bx, 0, m_bx, 0, mf.nComp());
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());

task_copy::~task_copy()
{
#ifdef BL_USE_MPI
    delete s_tmp;
    delete d_tmp;
#endif
}

bool task_copy::ready()
{
#ifdef BL_USE_MPI
    int flag;
    MPI_Status status;
    MPI_Test(&m_request, &flag, &status);
    return flag == 1;
#endif
    return m_ready;
}

task_copy_local::task_copy_local(FArrayBox& fab_, const MultiFab& smf_, int grid, const Box& bx)
: m_fab(fab_), m_smf(smf_), m_sgrid(grid), m_bx(bx)
{
    m_ready = true;
    m_fab.copy(m_smf[m_sgrid]);
}

task_copy_local::~task_copy_local()
{
}

bool task_copy_local::ready()
{
    return m_ready;
}

void task_list::execute()
{
    while ( !tasks.empty() )
    {
	task* t = tasks.front();
	tasks.pop_front();
	if ( t->ready() )
	{
	    delete t;
	}
	else
	{
	    tasks.push_back(t);
	}
    }
}

void task_list::add_task(task* t)
{
    tasks.push_back(t);
}

task_fab_get::task_fab_get(const MultiFab& r_, int grid_) 
: r(r_), grid(grid_), bx(r_.box(grid_)) {}

task_fab_get::task_fab_get(const MultiFab& r_, int grid_, const Box& bx_) 
: r(r_), grid(grid_), bx(bx_) {}

const FArrayBox& task_fab_get::fab()
{
    return r[grid];
}

task_fab_get::~task_fab_get()
{
}

bool task_fab_get::ready()
{
    return true;
}

