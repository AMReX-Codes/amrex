#include "hgparallel.h"

// TASK_COPY
task_copy::task_copy(MultiFab& mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
: m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx)
#ifdef BL_USE_MPI
, s_tmp(0), d_tmp(0)
#endif
{
}

task_copy::~task_copy()
{
#ifdef BL_USE_MPI
    delete s_tmp;
    delete d_tmp;
#endif
}

void task_copy::init(sequence_number sno, MPI_Comm comm)
{
#ifdef BL_USE_MPI
    assert( is_local(m_mf, m_dgrid) || is_local(m_smf, m_sgrid) );
    if ( is_local(m_mf, m_dgrid) && is_local(m_smf, m_sgrid) )
    {
	m_ready = true;
    }
    else if ( is_local(m_mf, m_dgrid) )
    {
	d_tmp = new FArrayBox(m_bx, m_mf.nComp());
	MPI_Irecv(d_tmp->dataPtr(), m_bx.numPts()*s_tmp->nComp(), MPI_DOUBLE, sno, processor_number(m_smf, m_sgrid), comm, &m_request);
    }
    else if ( is_local(m_smf, m_sgrid) ) 
    {
	s_tmp = new FArrayBox(m_bx, m_mf.nComp());
	s_tmp->copy(m_smf[m_sgrid]);
	MPI_Isend(s_tmp->dataPtr(), m_bx.numPts()*s_tmp->nComp(), MPI_DOUBLE, sno, processor_number(m_mf,  m_dgrid), comm, &m_request);
    }
#else
    m_ready = true;
#endif
}

bool task_copy::is_off_processor() const
{
    return is_remote(m_mf, m_dgrid) && is_remote(m_smf, m_sgrid);
}

task_copy::task_copy(MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
: m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), s_bx(sb), m_sgrid(sgrid)
{
    m_ready = true;
    m_mf[m_dgrid].copy(m_smf[m_sgrid], s_bx, 0, m_bx, 0, mf.nComp());
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());

bool task_copy::ready()
{
#ifdef BL_USE_MPI
    int flag;
    MPI_Status status;
    MPI_Test(&m_request, &flag, &status);
    if ( flag )
    {
	return true;
    }
    return false;
#else
    m_mf[m_dgrid].copy(m_smf[m_sgrid], m_bx);
    return true;
#endif
}

// TASK_COPY_LOCAL

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

// TASK_LIST

task_list::task_list(MPI_Comm comm_) 
{
    seq_no = 0;
    MPI_Comm_dup(comm_, &comm);
}

task_list::~task_list()
{
    seq_no = 0;
    MPI_Comm_free(&comm);
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
    seq_no++;
    if ( t->is_off_processor() )
    {
	delete t;
    }
    else
    {
	t->init(seq_no, comm);
	tasks.push_back(t);
    }
}

// TASK_FAB_GET

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

