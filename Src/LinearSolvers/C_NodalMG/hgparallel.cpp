#include "amr_defs.H"
#include "hgparallel.h"

#include <typeinfo>

bool HG_is_debugging = false;

bool task::depend_ready()
{
    list< task** >::iterator lit = dependencies.begin();
    while ( lit != dependencies.end() )
    {
	task** t = *lit;
	if ( *t == 0 )
	{
	    list< task** >::iterator tmp = lit++;
	    dependencies.erase(tmp);
	}
	else
	{
	    lit++;
	}
    }
    return dependencies.empty();
}

bool task::recommit(list<task*>*)
{
    return false;
}

void task::_hint() const
{
    HG_DEBUG_OUT( 
	"(" << typeid(*this).name() << ' ' << m_sno << ' ' << m_started << ' ' << m_comm << ' '
	);
}

void task::hint() const
{
    task::_hint();
    HG_DEBUG_OUT(")" << endl);
}

// TASK_LIST

bool task_list::def_verbose = true;

task_list::task_list(MPI_Comm comm_)
    : seq_no(0), seq_delta(10), verbose(def_verbose)
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

void task_list::add_task(task* t)
{
    seq_no += seq_delta;
    add_task(t, seq_no);
}

void task_list::add_task(task* t, task::sequence_number seq_no_)
{
    if ( ! t->init(seq_no_, comm) )
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
		t->depend_on(**cit);
	    }
	    cit++;
	}
	tasks.push_back( tp );
    }
}

void task_list::execute()
{
    if ( HG_is_debugging ) MPI_Barrier(comm);
    int l_progress = tasks.size() * 3;
    HG_DEBUG_OUT("Processing List " << comm << " with " << tasks.size() << " elements " << endl);
    list< task** > dead_tasks;
    // The dead_task list is used, because the tasks being processed also appear
    // in tasks dependecy lists.
    while ( !tasks.empty() )
    {
	task** t = tasks.front();
	tasks.pop_front();
	if ( verbose ) 
	    (*t)->hint();
	HG_DEBUG_OUT("*** Trying " << t << endl;)
	if ( (*t)->ready() )
	{
	    HG_DEBUG_OUT("*** Finished " << t << endl);
	    list<task*> tl;
	    if ( (*t)->recommit(&tl) )
	    {
		int i = 0;
		for(list<task*>::iterator tli = tl.begin(); tli != tl.end(); ++tli)
		{
		    add_task(*tli, (*tli)->get_sequence_number() + ++i);
		}
	    }
	    delete *t;
	    *t = 0;
	    dead_tasks.push_back(t);
	}
	else
	{
	    HG_DEBUG_OUT("*** Retry " << t << endl);
	    tasks.push_back(t);
	}
	if ( l_progress-- < 0 )
	{
	    BoxLib::Error("task_list::execute(): No Progress");
	}
    }
    while ( !dead_tasks.empty() )
    {
	delete dead_tasks.front();
	dead_tasks.pop_front();
    }
}

bool task_list::execute_no_block()
{
    HG_DEBUG_OUT("No Block Processing List " << comm << " with " << tasks.size() << " elements " << endl);
    list< task** > dead_tasks;
    list< task** >::iterator tli = tasks.begin();
    while ( tli != tasks.end() )
    {
	task** t = *tli;
	if ( verbose )
	    (*t)->hint();
	HG_DEBUG_OUT("*** No Block Trying " << t << endl;)
	if ( (*t)->ready() )
	{
	    HG_DEBUG_OUT("*** No Block Finished " << t << endl);
	    list<task*> tl;
	    if ( (*t)->recommit(&tl) )
	    {
		int i = 0;
		for(list<task*>::iterator tli = tl.begin(); tli != tl.end(); ++tli)
		{
		    add_task(*tli, (*tli)->get_sequence_number() + ++i);
		}
	    }
	    delete *t;
	    *t = 0;
	    dead_tasks.push_back(t);
	    list< task** >::iterator tmp = tli++;
	    tasks.erase(tmp);
	}
	else
	{
	    tli++;
	}
    }
    while ( !dead_tasks.empty() )
    {
	delete dead_tasks.front();
	dead_tasks.pop_front();
    }
    return tasks.empty();
}

// TASK_COPY
task_copy::task_copy(MultiFab& mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
    : m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx), m_sbx(bx), m_local(false), tmp(0), m_request(MPI_REQUEST_NULL)
{
}

task_copy::task_copy(MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
    : m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), m_sbx(sb), m_sgrid(sgrid), m_local(false), tmp(0), m_request(MPI_REQUEST_NULL)
{
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());
static bool eq(const MultiFab& a, const MultiFab& b)
{
    return &a == &b;
}


bool task_copy::depends_on_q(const task* t1) const
{
    if ( !eq(m_mf, m_smf) ) return false;
    if ( const task_copy* t1tc = dynamic_cast<const task_copy*>(t1) )
    {
	const Box& t1_bx = t1tc->m_bx;
	if ( m_bx.intersects(t1_bx) || m_sbx.intersects(t1_bx) ) return true;
    }
    return false;
}


task_copy::~task_copy()
{
    delete tmp;
    assert( m_request == MPI_REQUEST_NULL);
}

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

void task_copy::startup()
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
	tmp->copy(m_smf[m_sgrid], m_sbx);
	int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_mf,  m_dgrid), m_sno, m_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else
    {
	BoxLib::Abort( "task_copy::ready(): Can't be here" );
    }
    m_started = true;
}

bool task_copy::ready()
{
    if ( ! depend_ready() ) return false;
    if ( ! m_started ) startup();
    if ( m_local )
    {
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
	    m_mf[m_dgrid].copy(*tmp, m_sbx, 0, m_bx, 0, m_smf.nComp());
	}
	return true;
    }
    return false;
}

void task_copy::hint() const
{
    task::_hint();
    if ( is_local(m_smf, m_sgrid) && is_local(m_mf, m_dgrid))
    {
	HG_DEBUG_OUT( "L" );
    }
    else if ( is_local(m_smf, m_sgrid) )
    {
	HG_DEBUG_OUT( "S" );
    }
    else if ( is_local(m_mf, m_dgrid) )
    {
    	HG_DEBUG_OUT( "R" );
    }
    else
    {
	HG_DEBUG_OUT( "?" );
    }
    HG_DEBUG_OUT(
	' ' <<
	m_bx  << ' ' << m_dgrid << ' ' <<
	m_sbx  << ' ' << m_sgrid << ' '
	);
    HG_DEBUG_OUT( ")" << endl );
}

// TASK_COPY_LOCAL

task_copy_local::task_copy_local(FArrayBox* fab_, const Box& bx, const MultiFab& smf_, int grid)
    : m_fab(fab_), m_smf(smf_), m_sgrid(grid), m_bx(bx), tmp(0), m_local(false)
{
}

task_copy_local::~task_copy_local()
{
    if ( tmp ) HG_DEBUG_OUT("task_copy_local::~task_copy_local(): delete tmp" << endl);
    delete tmp;
}

bool task_copy_local::init(sequence_number sno, MPI_Comm comm)
{
    task::init( sno, comm);
    assert ( m_fab ==0 || m_fab->nComp() == m_smf.nComp() );
    if ( m_fab != 0 || is_local(m_smf, m_sgrid) ) return true;
    return false;
}

void task_copy_local::hint() const
{
    task::_hint();
    if ( m_fab !=0 && is_local(m_smf, m_sgrid) ) HG_DEBUG_OUT( "L" );
    else if ( m_fab != 0) HG_DEBUG_OUT("R");
    else if ( is_local( m_smf, m_sgrid ) ) HG_DEBUG_OUT("S");
    else HG_DEBUG_OUT("?");
    HG_DEBUG_OUT(
	m_bx <<  ' ' <<  m_sgrid << ' '
	);
    HG_DEBUG_OUT( ")" << endl );
}

void task_copy_local::startup()
{
    if ( m_fab !=0 && is_local(m_smf, m_sgrid) )
    {
	m_local = true;
    }
    else if ( m_fab != 0 )
    {
	tmp = new FArrayBox(m_bx, m_smf.nComp());
	int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, m_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else if ( is_local(m_smf, m_sgrid) ) 
    {
	tmp = new FArrayBox(m_bx, m_smf.nComp());
	tmp->copy(m_smf[m_sgrid], m_bx);
	int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, m_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else
    {
	BoxLib::Abort( "task_copy_local::ready(): Can't Happen" );
    }
    m_started = true;
}

bool task_copy_local::ready()
{
    if ( ! depend_ready() ) return false;
    if ( ! m_started ) startup();
    if ( m_local )
    {
	m_fab->copy(m_smf[m_sgrid], m_bx);
	return true;
    }
    int flag;
    MPI_Status status;
    assert ( m_request != MPI_REQUEST_NULL );
    int res = MPI_Test(&m_request, &flag, &status);
    if ( res != 0 )
	ParallelDescriptor::Abort( res );
    HG_DEBUG_OUT("task_copy_local::ready(): " << flag << endl);
    if ( flag )
    {
	assert ( m_request == MPI_REQUEST_NULL );
	if ( m_fab )
	{
	    int count;
	    assert( status.MPI_SOURCE == processor_number(m_smf, m_sgrid) );
	    assert( status.MPI_TAG    == m_sno );
	    int res = MPI_Get_count(&status, MPI_DOUBLE, &count);
	    if ( res != 0 )
		ParallelDescriptor::Abort( res );
	    assert(count == tmp->box().numPts()*tmp->nComp());
	    m_fab->copy(*tmp, m_bx);
	}
	return true;
    }
    return false;
}

// TASK_FAB

const FArrayBox& task_fab::fab()
{
    assert(target != 0);
    return *target;
}

bool task_fab::init(sequence_number sno, MPI_Comm comm)
{
    task::init(sno, comm);
    if ( m_local_target )
    {
	target = new FArrayBox(region, ncomp);
	return true;
    }
    return false;
}


// task_fec_base

task_fec_base::task_fec_base(const list<int>& tll_, const Box& freg_, MultiFab& s_, int igrid_)
    : tll(tll_), freg(freg_), s(s_), igrid(igrid_)
{
}
task_fec_base::task_fec_base( MultiFab& s_, int igrid_)
: s(s_), igrid(igrid_)
{
}
task_fec_base::~task_fec_base()
{
    HG_DEBUG_OUT("task_fec_base::~task_fec_base()" << endl);
    for( vector<task_fab*>::iterator tfi = tfvect.begin(); tfi != tfvect.end(); ++tfi)
    {
	delete *tfi;
    }
}

bool task_fec_base::init(sequence_number sno, MPI_Comm comm)
{
    task::init(sno, comm);
    bool result = is_local(s, igrid);
    for(vector<task_fab*>::iterator tli = tfvect.begin(); tli != tfvect.end(); ++tli)
    {
	bool tresult = (*tli)->init(sno, comm);
	result = tresult ||  result;
    }
    for ( list<int>::const_iterator tli = tll.begin(); tli != tll.end(); ++tli )
    {
	bool tresult = is_local(s, *tli);
	result = result || tresult;
    }
    HG_DEBUG_OUT("task_fec_base::init(): result = " << result << endl);
    return result;
}

bool task_fec_base::ready()
{
    bool result = true;
    for(vector<task_fab*>::iterator tfi = tfvect.begin(); tfi != tfvect.end(); ++tfi)
    {
	bool tresult = (*tfi)->ready();
	result = tresult && result;
    }
    return result;
}

void task_fec_base::push_back(task_fab* tf)
{
    tfvect.push_back(tf);
}

bool task_fec_base::is_local_target() const
{
    return is_local(s, igrid);
}

FArrayBox& task_fec_base::target_fab()
{
    assert ( is_local_target() );
    return s[igrid];
}

int task_fec_base::grid_number() const
{
    return igrid;
}

const FArrayBox& task_fec_base::task_fab_result(int n)
{
    assert ( n >= 0 && n < tfvect.size() );
    return tfvect[n]->fab();
}

bool task_fec_base::recommit(list<task*>* tl)
{
    for(list<int>::const_iterator tlli = tll.begin(); tlli != tll.end(); ++tlli)
    {
	tl->push_back(
	    new task_copy(s, *tlli, s, igrid, freg)
	    );
    }
    return tl->size() > 0;
}

void task_fec_base::clean_up()
{
}
