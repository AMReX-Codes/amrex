#include "amr_defs.H"
#include "hgparallel.h"

bool task::depend_ready()
{
    list< task** >::iterator lit = dependencies.begin();
    while ( lit != dependencies.end() )
    {
	if ( **lit )
	    return false;
    }
    return true;
}

// TASK_COPY
task_copy::task_copy(MultiFab& mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
: m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx), m_sbx(bx), m_local(false)
#ifdef BL_USE_MPI
, tmp(0), m_request(MPI_REQUEST_NULL)
#endif
{
}

bool task_copy::depends_on_q(const task* t1) const
{
    const task_copy* t1tc = dynamic_cast<const task_copy*>(t1);
    if ( t1tc )
    {
    }
    else
    {
	BoxLib::Abort("Nightmare");
    }
    return false;
}

task_copy::task_copy(MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
: m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), m_sbx(sb), m_sgrid(sgrid), m_local(false)
#ifdef BL_USE_MPI
, tmp(0), m_request(MPI_REQUEST_NULL)
#endif
{
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());

task_copy::~task_copy()
{
#ifdef BL_USE_MPI
    delete tmp;
    assert( m_request == MPI_REQUEST_NULL);
#endif
}

#if 0
bool task_copy::init(sequence_number sno, MPI_Comm comm)
{
    task::init(sno, comm);
    assert( m_sbx.numPts() == m_bx.numPts() );
    if ( is_local(m_mf, m_dgrid) || is_local(m_smf, m_sgrid) )
    {
	return true;
    }
    else
    {
	return false;
    }
}

bool task_copy::ready()
{
    if ( ! depend_ready() ) return false;
    if ( ! m_started )
    {
	if ( is_local(m_mf, m_dgrid) && is_local(m_smf, m_sgrid) )
	{
	    m_local = true;
	}
	else if ( is_local(m_mf, m_dgrid) )
	{
	    tmp = new FArrayBox(m_sbx, m_smf.nComp());
	    int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, m_comm, &m_request);
	    if ( res != 0 )
		ParallelDescriptor::Abort(res);
	    assert( m_request != MPI_REQUEST_NULL );
	}
	else if ( is_local(m_smf, m_sgrid) ) 
	{
	    tmp = new FArrayBox(m_sbx, m_smf.nComp());
	    // before I can post the receive, I have to ensure that there are no dependent zones in the
	    // grid
	    tmp->copy(m_smf[m_sgrid], m_sbx);
	    HG_DEBUG_OUT( "<< Norm(S) of tmp " << m_sno << " " << tmp->norm(m_sbx, 2) << endl );
	    HG_DEBUG_OUT( "<<<Box(S) of tmp "   << m_sno << " " << tmp->box() << endl );
	    // printRange(debug_out, *tmp, m_sbx, 0, tmp->nComp());
	    int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_mf,  m_dgrid), m_sno, m_comm, &m_request);
	    if ( res != 0 )
		ParallelDescriptor::Abort(res);
	    assert( m_request != MPI_REQUEST_NULL );
	}
	else
	{
	    BoxLib::Error("task_copy::ready: Can't be here");
	    // neither fab lives on local processor
	}
	m_started = true;
    }
    if ( m_local )
    {
	HG_DEBUG_OUT( "Norm(L) " << m_sno << " " << m_smf[m_sgrid].norm(m_sbx, 2) << endl );
	// printRange(debug_out, m_smf[m_sgrid], m_sbx, 0, m_smf.nComp());
	m_mf[m_dgrid].copy(m_smf[m_sgrid], m_sbx, 0, m_bx, 0, m_mf.nComp());
	return true;
    }
    int flag;
    MPI_Status status;
    assert( m_request != MPI_REQUEST_NULL );
    int res = MPI_Test(&m_request, &flag, &status);
    if ( res != 0 )
	ParallelDescriptor::Abort( res );
    if ( flag )
    {
	assert ( m_request == MPI_REQUEST_NULL );
	if ( is_local(m_mf, m_dgrid) )
	{
	    int count;
	    assert( status.MPI_SOURCE == processor_number(m_smf, m_sgrid) );
	    assert( status.MPI_TAG    == m_sno );
	    int res = MPI_Get_count(&status, MPI_DOUBLE, &count);
	    if ( res != 0 )
		ParallelDescriptor::Abort( res );
	    assert(count == tmp->box().numPts()*tmp->nComp());
	    HG_DEBUG_OUT( ">> Norm(R) of tmp " << m_sno << " " << tmp->norm(m_sbx, 2) << endl );
	    HG_DEBUG_OUT( ">>>Box(R) of tmp "   << m_sno << " " << tmp->box() << endl );
	    // printRange(debug_out, *tmp, m_sbx, 0, tmp->nComp());
	    m_mf[m_dgrid].copy(*tmp, m_sbx, 0, m_bx, 0, m_smf.nComp());
	}
	return true;
    }
    return false;
}

#else

bool task_copy::init(sequence_number sno, MPI_Comm comm)
{
    m_sno = sno;
    assert( m_sbx.numPts() == m_bx.numPts() );
    // debug_out << "task_copy::init " << m_mf.nComp() << ' ' << m_smf.nComp() << endl;
#ifdef BL_USE_MPI
    if ( is_local(m_mf, m_dgrid) && is_local(m_smf, m_sgrid) )
    {
	m_local = true;
    }
    else if ( is_local(m_mf, m_dgrid) )
    {
	tmp = new FArrayBox(m_sbx, m_smf.nComp());
	int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else if ( is_local(m_smf, m_sgrid) ) 
    {
	tmp = new FArrayBox(m_sbx, m_smf.nComp());
	// before I can post the receive, I have to ensure that there are no dependent zones in the
	// grid
	tmp->copy(m_smf[m_sgrid], m_sbx);
	HG_DEBUG_OUT( "<< Norm(S) of tmp " << m_sno << " " << tmp->norm(m_sbx, 2) << endl );
	HG_DEBUG_OUT( "<<<Box(S) of tmp "   << m_sno << " " << tmp->box() << endl );
	// printRange(debug_out, *tmp, m_sbx, 0, tmp->nComp());
	int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_mf,  m_dgrid), m_sno, comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
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
    if ( !depend_ready() ) return false;
    if ( m_local )
    {
	HG_DEBUG_OUT( "Norm(L) " << m_sno << " " << m_smf[m_sgrid].norm(m_sbx, 2) << endl );
	// printRange(debug_out, m_smf[m_sgrid], m_sbx, 0, m_smf.nComp());
	m_mf[m_dgrid].copy(m_smf[m_sgrid], m_sbx, 0, m_bx, 0, m_mf.nComp());
	return true;
    }
#ifdef BL_USE_MPI
    int flag;
    MPI_Status status;
    assert( m_request != MPI_REQUEST_NULL );
    int res = MPI_Test(&m_request, &flag, &status);
    if ( res != 0 )
	ParallelDescriptor::Abort( res );
    if ( flag )
    {
	assert ( m_request == MPI_REQUEST_NULL );
	if ( is_local(m_mf, m_dgrid) )
	{
	    int count;
	    assert( status.MPI_SOURCE == processor_number(m_smf, m_sgrid) );
	    assert( status.MPI_TAG    == m_sno );
	    int res = MPI_Get_count(&status, MPI_DOUBLE, &count);
	    if ( res != 0 )
		ParallelDescriptor::Abort( res );
	    assert(count == tmp->box().numPts()*tmp->nComp());
	    HG_DEBUG_OUT( ">> Norm(R) of tmp " << m_sno << " " << tmp->norm(m_sbx, 2) << endl );
	    HG_DEBUG_OUT( ">>>Box(R) of tmp "   << m_sno << " " << tmp->box() << endl );
	    // printRange(debug_out, *tmp, m_sbx, 0, tmp->nComp());
	    m_mf[m_dgrid].copy(*tmp, m_sbx, 0, m_bx, 0, m_smf.nComp());
	}
	return true;
    }
    return false;
#endif
}

#endif

void task_copy::hint() const
{
    HG_DEBUG_OUT( "task_copy : ");
    if ( m_local )
    {
	HG_DEBUG_OUT( "L" );
    }
    else if ( is_local(m_smf, m_sgrid) )
    {
	HG_DEBUG_OUT( "S" );
    }
    else
    {
    	HG_DEBUG_OUT( "R" );
    }
    HG_DEBUG_OUT( 
	' ' <<
	m_sno << ' ' <<
	m_bx  << ' ' << m_dgrid << ' ' <<
	m_sbx  << ' ' << m_sgrid << ' ' <<
	endl );	// to flush
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

bool task_list::def_verbose = true;

task_list::task_list(MPI_Comm comm_)
    : seq_no(0), verbose(def_verbose)
{
    int res = MPI_Comm_dup(comm_, &comm);
    if ( res != 0 )
	ParallelDescriptor::Abort( res );
}

task_list::~task_list()
{
    seq_no = 0;
    int res = MPI_Comm_free(&comm);
    if ( res != 0 )
	ParallelDescriptor::Abort( res );
}

void task_list::execute()
{
    list< task** > dead_tasks;
    while ( !tasks.empty() )
    {
	task** t = tasks.front();
	tasks.pop_front();
	if ( verbose )
	    (*t)->hint();
	if ( (*t)->ready() )
	{
	    delete *t;
	    *t = 0;
	    dead_tasks.push_back(t);
	}
	else
	{
	    tasks.push_back(t);
	}
    }
    while ( !dead_tasks.empty() )
    {
	delete dead_tasks.front();
	dead_tasks.pop_front();
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
	task** tp = new task*( t );
	// loop here over existing tasks, see if the depend on this one,
	// if so, add them to the dependency
	list< task**>::const_iterator cit = tasks.begin();
	while ( cit != tasks.end() )
	{
	    if ( t->depends_on_q( **cit ) )
	    {
		// do something
	    }
	    cit++;
	}
	tasks.push_back( tp );
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

