
#include "hgparallel.h"

// TASK_COPY
task_copy::task_copy(MultiFab& mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
: m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx), s_bx(bx), m_local(false)
#ifdef BL_USE_MPI
, tmp(0)
#endif
{
}

task_copy::task_copy(MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
: m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), s_bx(sb), m_sgrid(sgrid), m_local(false)
#ifdef BL_USE_MPI
, tmp(0)
#endif
{
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());

task_copy::~task_copy()
{
#ifdef BL_USE_MPI
    delete tmp;
#endif
}

bool task_copy::init(sequence_number sno, MPI_Comm comm)
{
#ifdef BL_USE_MPI
    if ( is_local(m_mf, m_dgrid) && is_local(m_smf, m_sgrid) )
    {
	m_local = true;
    }
    else if ( is_local(m_mf, m_dgrid) )
    {
	tmp = new FArrayBox(m_bx, m_mf.nComp());
	int res = MPI_Irecv(tmp->dataPtr(), m_bx.numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), sno, comm, &m_request);
	if ( res != 0 )
	    BoxLib::Error("Failed MPI_Irecv");
    }
    else if ( is_local(m_smf, m_sgrid) ) 
    {
	tmp = new FArrayBox(s_bx, m_mf.nComp());
	tmp->copy(m_smf[m_sgrid]);
	int res = MPI_Isend(tmp->dataPtr(), s_bx.numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_mf,  m_dgrid), sno, comm, &m_request);
	if ( res != 0 )
	    BoxLib::Error("Failed MPI_Isend");
    }
    else
    {
	// neither fab lives on local processor
	return false;
    }
#else
    m_local = true;
#endif
    return true;
}


bool task_copy::ready()
{
    if ( m_local )
    {
	m_mf[m_dgrid].copy(m_smf[m_sgrid], s_bx, 0, m_bx, 0, m_mf.nComp());
	return true;
    }
#ifdef BL_USE_MPI
    int flag;
    MPI_Status status;
    int res = MPI_Test(&m_request, &flag, &status);
    if ( res != 0 )
	BoxLib::Error("Failed MPI_Test");
    if ( flag )
    {
	if ( is_local(m_mf, m_dgrid) )
	    m_mf[m_dgrid].copy(*tmp, s_bx, 0, m_bx, 0, m_mf.nComp());
	return true;
    }
    return false;
#endif
}

// TASK_COPY_LOCAL

task_copy_local::task_copy_local(FArrayBox& fab_, const MultiFab& smf_, int grid, const Box& bx)
: m_fab(fab_), m_smf(smf_), m_sgrid(grid), m_bx(bx)
{
    m_local = true;
    m_fab.copy(m_smf[m_sgrid]);
}

task_copy_local::~task_copy_local()
{
}

bool task_copy_local::ready()
{
    abort(); return m_local;
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
    if ( ! t->init(seq_no, comm) )
    {
	delete t;
    }
    else
    {
	tasks.push_back(t);
    }
}

// TASK_FAB_GET

task_fab_get::task_fab_get(const MultiFab& d_, int dgrid_, const MultiFab& s_, int sgrid_, const Box& bx_) 
: d(d_), dgrid(dgrid_), s(s_), sgrid(sgrid_), bx(bx_) {}

const FArrayBox& task_fab_get::fab()
{
    return s[sgrid];
}

task_fab_get::~task_fab_get()
{
}

bool task_fab_get::ready()
{
    abort(); return true;
}

