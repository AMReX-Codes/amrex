#include "hgparallel.h"

task_copy::task_copy(MultiFab* mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
: m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx)
{
    m_ready = true;
    m_mf->operator[](m_dgrid).copy(m_smf[m_sgrid], m_bx);
}

task_copy::task_copy(FArrayBox* fab_, const MultiFab& smf_, int grid, const Box& bx)
: m_mf(0), m_fab(fab_), m_smf(smf_), m_sgrid(grid), m_bx(bx)
{
    m_ready = true;
    m_fab->copy(m_smf[m_sgrid]);
}

task_copy::task_copy(MultiFab* mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
: m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), s_bx(sb), m_sgrid(sgrid)
{
    m_ready = true;
    m_mf->operator[](m_dgrid).copy(m_smf[m_sgrid], s_bx, 0, m_bx, 0, mf->nComp());
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());

task_copy::~task_copy()
{
}

bool
task_copy::ready()
{
    return m_ready;
}

void
task_list::execute()
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

void
task_list::add_task(task* t)
{
    tasks.push_back(t);
}