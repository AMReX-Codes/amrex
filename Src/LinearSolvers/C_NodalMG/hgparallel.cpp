#include "amr_defs.H"
#include "hgparallel.h"

#ifdef HG_DEBUG
#ifdef BL_USE_NEW_HFILES
#include <typeinfo>
#else
#include <typeinfo.h>
#endif
#endif

bool HG_is_debugging = false;
MPI_Comm HG::mpi_comm = MPI_COMM_WORLD;
int HG::mpi_tag_ub;

void HG::MPI_init()
{
    static int first = 0;
    if ( first == 0 )
    {
	first = 1;
	int res;
	// int flag;
	// res = MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &mpi_tag_ub, &flag);
	// if ( res != 0  )
	//    ParallelDescriptor::Abort( res );
	res = MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
	if ( res != 0 )
	    ParallelDescriptor::Abort( res );
    }
}

void HG::MPI_finish()
{
}

// TASK

task::task(task_list& tl_) 
    : m_sno(tl_.get_then_advance()), m_started(false), m_task_list(tl_) 
{
    assert(m_sno != 0);
}

task::~task() 
{
    assert(m_sno != 0);
}

bool task::is_started() const
{
    return m_started;
}

bool task::startup()
{
    m_started = true;
    return true;
}

task::sequence_number task::get_sequence_number() const
{
    return m_sno;
}

void task::print_dependencies(ostream& os) const
{
    os << "Task " << get_sequence_number() << " depends on ( ";
    for ( list<task_proxy>::const_iterator lit = dependencies.begin(); lit != dependencies.end(); ++lit)
    {
	if ( ! (*lit).is_finished() )
	{
	    os << (*lit)->get_sequence_number() << " ";
	}
    }
    os << ") ";
}

bool task::depend_ready()
{
    for ( list<task_proxy>::iterator lit = dependencies.begin(); lit != dependencies.end(); )
    {
	task_proxy t = *lit;
	if ( t.is_finished() )
	{
	    list< task_proxy >::iterator tmp = lit++;
	    dependencies.erase(tmp);
	}
	else
	{
	    lit++;
	}
    }
    return dependencies.empty();
}

void task::_hint() const
{
#ifdef HG_DEBUG
    assert(m_sno != 0);
    HG_DEBUG_OUT( 
	"(" << typeid(*this).name() << ' ' << m_sno << ' ' << m_started << ' '
	);
    print_dependencies(debug_out);
#endif
}

void task::hint() const
{
    task::_hint();
    HG_DEBUG_OUT(")" << endl);
}

bool task::depends_on_q(const task* t1) const
{
    return false;
}

void task::depend_on(const task_proxy& t1)
{
    dependencies.push_back( t1 );
}

bool task::ready()
{
    assert(is_started());
    return true;
}

// The list...
// TASK_LIST

bool task_list::def_verbose = false;

task_list::task_list()
    : seq_no(1), seq_delta(1), verbose(def_verbose)
{
}

task_list::~task_list()
{
    seq_no = 1;
}

task::sequence_number task_list::get_then_advance()
{
    task::sequence_number tmp = seq_no;
    seq_no += seq_delta;
    return tmp;
}

task::task_proxy task_list::add_task(task* t)
{
    task::task_proxy result(t);
    tasks.push_back( result );
    return result;
}

void task_list::print_dependencies(ostream& os) const
{
    os << "Task list ( " << endl;
    for( list< task::task_proxy >::const_iterator tli = tasks.begin(); tli != tasks.end(); ++tli )
    {
	(*tli)->print_dependencies(os);
	os << endl;
    }
    os << ")" << endl;
}

void task_list::execute()
{
    if ( HG_is_debugging ) MPI_Barrier(HG::mpi_comm);
    int l_progress = tasks.size() * 3;
    if ( verbose )
    {
#ifdef HG_DEBUG
	HG_DEBUG_OUT("Processing List " << HG::mpi_comm << " with " << tasks.size() << " elements " << endl);
	print_dependencies(debug_out);
#endif
    }
    while ( !tasks.empty() )
    {
	task::task_proxy t = tasks.front();
	tasks.pop_front();
	if ( verbose ) 
	{
	    t->hint();
	}
	if ( t->depend_ready() )
	{
	    if ( ! t->is_started() )
	    {
		if ( ! t->startup() )
		{
		    t.set_finished();
		    continue;
		}
	    }
	    assert(t->is_started());
	    if ( t->depend_ready() && t->ready() )
	    {
		t.set_finished();
		continue;
	    }
	}
        tasks.push_back(t);
	if ( l_progress-- < 0 )
	{
	    // BoxLib::Error("task_list::execute(): No Progress");
	}
    }
    seq_no = 1;
}

bool task_list::empty() const 
{ 
    return tasks.empty(); 
}
int task_list::size() const 
{
    return tasks.size(); 
}

list<task::task_proxy>::const_iterator task_list::begin() const
{
    return tasks.begin();
}

list<task::task_proxy>::const_iterator task_list::end() const
{
    return tasks.end();
}

list<task::task_proxy>::iterator task_list::begin()
{
    return tasks.begin();
}

list<task::task_proxy>::iterator task_list::end()
{
    return tasks.end();
}


// TASK_COPY
task_copy::task_copy(task_list& tl_, MultiFab& mf, int dgrid, const MultiFab& smf, int sgrid, const Box& bx)
    : task(tl_), m_mf(mf), m_dgrid(dgrid), m_smf(smf), m_sgrid(sgrid), m_bx(bx), m_sbx(bx), m_local(false), tmp(0), m_request(MPI_REQUEST_NULL)
{
    _do_depend();
}

task_copy::task_copy(task_list& tl_, MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb)
    : task(tl_), m_mf(mf), m_bx(db), m_dgrid(dgrid), m_smf(smf), m_sbx(sb), m_sgrid(sgrid), m_local(false), tmp(0), m_request(MPI_REQUEST_NULL)
{
    _do_depend();
}
			// r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());
static inline bool eq(const MultiFab& a, const MultiFab& b)
{
    return &a == &b;
}


bool task_copy::depends_on_q(const task* t1) const
{
    if ( !eq(m_mf, m_smf) ) return false;
    if ( const task_copy* t1tc = dynamic_cast<const task_copy*>(t1) )
    {
	if ( m_sgrid == t1tc->m_dgrid && m_sbx.intersects(t1tc->m_bx) ) return true;
	if ( m_dgrid == t1tc->m_dgrid && m_bx.intersects(t1tc->m_bx) ) return true;
	if ( m_sgrid == t1tc->m_sgrid && m_sbx.intersects(t1tc->m_sbx) ) return true;
	if ( m_dgrid == t1tc->m_sgrid && m_bx.intersects(t1tc->m_sbx) ) return true;
    }
    return false;
}

void task_copy::_do_depend()
{
    
    for ( list< task::task_proxy >::const_iterator cit = m_task_list.begin(); cit != m_task_list.end(); ++cit)
    {
	if ( depends_on_q( **cit ) )
	{
	    depend_on(*cit);
	}
    }
}
  
task_copy::~task_copy()
{
    delete tmp;
    assert( m_request == MPI_REQUEST_NULL);
}

bool task_copy::startup()
{
    bool result = true;
    m_started = true;
    if ( is_local(m_mf, m_dgrid) && is_local(m_smf, m_sgrid) )
    {
	m_local = true;
    }
    else if ( is_local(m_mf, m_dgrid) )
    {
	tmp = new FArrayBox(m_sbx, m_smf.nComp());
	int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, HG::mpi_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else if ( is_local(m_smf, m_sgrid) ) 
    {
	tmp = new FArrayBox(m_sbx, m_smf.nComp());
	tmp->copy(m_smf[m_sgrid], m_sbx);
	int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_mf,  m_dgrid), m_sno, HG::mpi_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else
    {
	result = false;
    }
    return result;
}

bool task_copy::ready()
{
    assert(is_started());
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
	    m_mf[m_dgrid].copy(*tmp, tmp->box(), 0, m_bx, 0, m_smf.nComp());
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
	'(' << m_dgrid << "," << m_sgrid << ')' <<
	m_sbx << ' ' <<
	m_bx  << ' '
	);
    HG_DEBUG_OUT( ")" << endl );
}

// TASK_COPY_LOCAL

task_copy_local::task_copy_local(task_list& tl_, FArrayBox* fab_, int target_proc_id, const Box& bx, const MultiFab& smf_, int grid)
    : task(tl_), m_fab(fab_), m_smf(smf_), m_sgrid(grid), m_bx(bx), tmp(0), m_local(false), m_target_proc_id(target_proc_id)
{
}

task_copy_local::~task_copy_local()
{
    if ( tmp ) HG_DEBUG_OUT("task_copy_local::~task_copy_local(): delete tmp" << endl);
    delete tmp;
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

bool task_copy_local::startup()
{
    bool result = true;
    m_started = true;
    if ( m_fab !=0 && is_local(m_smf, m_sgrid) )
    {
	m_local = true;
    }
    else if ( m_fab != 0 )
    {
	tmp = new FArrayBox(m_bx, m_smf.nComp());
	int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, HG::mpi_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else if ( is_local(m_smf, m_sgrid) ) 
    {
	tmp = new FArrayBox(m_bx, m_smf.nComp());
	tmp->copy(m_smf[m_sgrid], m_bx);
	int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, m_target_proc_id, m_sno, HG::mpi_comm, &m_request);
	if ( res != 0 )
	    ParallelDescriptor::Abort(res);
	assert( m_request != MPI_REQUEST_NULL );
    }
    else
    {
	result = false;
    }
    return result;
}

bool task_copy_local::ready()
{
    assert(is_started());
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

task_fab::task_fab(task_list& tl_, const MultiFab&t_, int tt_, const Box& region_, int ncomp_)
    : task(tl_), m_target_proc_id(processor_number(t_,tt_)), region(region_), ncomp(ncomp_), target(0) 
{
    assert(m_sno != 0);
    if ( is_local(t_, tt_) )
    {
	target = new FArrayBox(region, ncomp);
    }
}

task_fab::~task_fab()
{
    delete target;
}

int task_fab::target_proc_id() const
{
    return m_target_proc_id;
}

const FArrayBox& task_fab::fab()
{
    assert(target != 0);
    return *target;
}

// task_fec_base

task_fec_base::task_fec_base( task_list& tl_, MultiFab& s_, int igrid_)
    : task(tl_), s(s_), igrid(igrid_)
{
}

task_fec_base::~task_fec_base()
{
    HG_DEBUG_OUT("task_fec_base::~task_fec_base()" << endl);
}

void task_fec_base::push_back(task_fab* tf)
{
    task_proxy tp(m_task_list.add_task(tf));
    tfvect.push_back(tp);
    depend_on( tp );
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
    task_fab* tf = dynamic_cast<task_fab*>(tfvect[n].get());
    return tf->fab();
}
