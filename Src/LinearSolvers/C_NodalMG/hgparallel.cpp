//BL_COPYRIGHT_NOTICE

#include "amr_defs.H"
#include "hgparallel.H"

#include <ParmParse.H>

#ifdef HG_DEBUG
#ifdef BL_USE_NEW_HFILES
#include <typeinfo>
#else
#include <typeinfo.h>
#endif
#endif

bool HG_is_debugging       = false;
#ifdef BL_USE_MPI
MPI_Comm HG::mpi_comm      = MPI_COMM_WORLD;
#endif
int HG::max_live_tasks     = 50;
int HG::multigrid_maxiter  = 100;
int HG::cgsolve_maxiter    = 250;
int HG::pverbose           = 0;
bool HG::initialized       = false;
double HG::cgsolve_tolfact = 1.0e-3;

void HG::MPI_init ()
{
    if (!initialized)
    {
        ParmParse pp("HG");
        pp.query("cgsolve_maxiter", cgsolve_maxiter);
        pp.query("multigrid_maxiter", multigrid_maxiter);
        pp.query("cgsolve_tolfact", cgsolve_tolfact);
        pp.query("max_live_tasks", max_live_tasks);
        pp.query("pverbose", pverbose);
        if (ParallelDescriptor::IOProcessor())
        {
            cout << "HG.cgsolve_maxiter = " << cgsolve_maxiter << endl;
            cout << "HG.multigrid_maxiter = " << multigrid_maxiter << endl;
            cout << "HG.cgsolve_tolfact = " << cgsolve_tolfact << endl;
            cout << "HG.max_live_tasks = " << max_live_tasks << endl;
            cout << "HG.pverbose = " << pverbose << endl;
        }
        initialized = true;
#ifdef BL_USE_MPI
        int res;
        // int flag;
        // res = MPI_Attr_get(MPI_COMM_WORLD, MPI_TAG_UB, &mpi_tag_ub, &flag);
        // if (res != 0 )
        //    ParallelDescriptor::Abort(res);
        res = MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
        if (res != 0)
            ParallelDescriptor::Abort(res);
#endif
    }
}

void HG::MPI_finish () {}

task::task (task_list& tl_) 
    :
    m_task_list(tl_),
    m_sno(tl_.get_then_advance()),
    m_cnt(1),
    m_finished(false),
    m_started(false)
{
    BLassert(m_sno != 0);
}

task::~task () 
{
    BLassert(m_sno != 0);
}

bool
task::startup ()
{
    return m_started = true;
}

bool
task::work_to_do () const
{
    return true;
}

void
task::print_dependencies (ostream& os) const
{
    os << "Task " << get_sequence_number() << " depends on ( ";
    for (list<task_proxy>::const_iterator lit = dependencies.begin(); lit != dependencies.end(); ++lit)
    {
        if (!(*lit).is_finished())
        {
            os << (*lit)->get_sequence_number() << " ";
        }
    }
    os << ") ";
}

bool
task::depend_ready ()
{
    for (list<task_proxy>::iterator lit = dependencies.begin(); lit != dependencies.end();)
    {
        if ((*lit).is_finished())
        {
            dependencies.erase(lit++);
        }
        else
        {
            ++lit;
        }
    }

    return dependencies.empty();
}

void
task::_hint () const
{
#ifdef HG_DEBUG
    BLassert(m_sno != 0);
    HG_DEBUG_OUT("(" << typeid(*this).name() << ' ' << m_sno << ' ' << m_started << ' ');
    print_dependencies(debug_out);
#endif
}

void
task::hint () const
{
    task::_hint();
    HG_DEBUG_OUT(")" << endl);
}

bool
task::depends_on_q (const task* t1) const
{
    return false;
}

void
task::depend_on (const task_proxy& t1)
{
    if (!t1.null())
        dependencies.push_back(t1);
}

bool
task::ready ()
{
    BLassert(is_started());
    return true;
}

task_list::task_list ()
    :
    seq_no(1),
    seq_delta(1),
    verbose(HG::pverbose)
{}

task_list::~task_list ()
{
    seq_no = 1;
}

task::task_proxy task_list::add_task (task* t)
{
    BLassert(t != 0);

    if (t->work_to_do())
    {
        tasks.push_back(task::task_proxy(t));
        return tasks.back();
    }
    else
    {
        delete t;
        return task::task_proxy(0);
    }
}

void
task_list::print_dependencies (ostream& os) const
{
    os << "Task list ( " << endl;
    for (list<task::task_proxy>::const_iterator tli = tasks.begin(); tli != tasks.end(); ++tli)
    {
        (*tli)->print_dependencies(os);
        os << endl;
    }
    os << ")" << endl;
}

void
task_list::execute ()
{
#ifdef BL_USE_MPI
    if (HG_is_debugging)
        MPI_Barrier(HG::mpi_comm);
#endif

    if (verbose)
    {
#ifdef HG_DEBUG
        HG_DEBUG_OUT("Processing List ");
#ifdef BL_USE_MPI
        HG_DEBUG_OUT(HG::mpi_comm);
#endif
        HG_DEBUG_OUT(" with " << tasks.size() << " elements " << endl);
        print_dependencies(debug_out);
#endif
    }

    int live_tasks = 0;
restart:
    while (!tasks.empty())
    {
        list<task::task_proxy>::iterator tli = tasks.begin();

        while (tli != tasks.end())
        {
            task::task_proxy t = *tli;

            BLassert(!t.null());

            if (verbose)
                t->hint();
            if (t->depend_ready())
            {
                if (!t->is_started())
                {
                    if (live_tasks > HG::max_live_tasks)
                        goto restart;
                    if (!t->startup())
                    {
                        t.set_finished();
                        tasks.erase(tli++);
                        continue;
                    }
                    live_tasks++;
                }
                if (t->ready())
                {
                    t.set_finished();
                    live_tasks--;
                    tasks.erase(tli++);
                    continue;
                }
            }

            ++tli;
        }
    }
    BLassert(live_tasks == 0);

    seq_no = 1;
}

            // r[jgrid].copy(r[igrid], bb, 0, b, 0, r.nComp());

void
task_copy::init ()
{
    if (work_to_do())
    {
        _do_depend();

        if (is_local(m_mf,m_dgrid) && is_local(m_smf,m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_mf[m_dgrid].copy(m_smf[m_sgrid],m_sbx,0,m_bx,0,m_mf.nComp());
                //
                // Flip the work_to_do() bit.
                //
                m_done = true;
            }
        }
    }
}

task_copy::task_copy (task_list&      tl_,
                      MultiFab&       mf,
                      int             dgrid,
                      const MultiFab& smf,
                      int             sgrid,
                      const Box&      bx)
    :
    task(tl_),
#ifdef BL_USE_MPI
    m_request(MPI_REQUEST_NULL),
#endif
    tmp(0),
    m_mf(mf),
    m_smf(smf),
    m_dgrid(dgrid),
    m_sgrid(sgrid),
    m_bx(bx),
    m_sbx(bx),
    m_local(false),
    m_done(false)
{
    init();
}

task_copy::task_copy (task_list&      tl_,
                      MultiFab&       mf,
                      int             dgrid,
                      const Box&      db,
                      const MultiFab& smf,
                      int             sgrid,
                      const Box&      sb)
    :
    task(tl_),
#ifdef BL_USE_MPI
    m_request(MPI_REQUEST_NULL),
#endif
    tmp(0),
    m_mf(mf),
    m_smf(smf),
    m_dgrid(dgrid),
    m_sgrid(sgrid),
    m_bx(db),
    m_sbx(sb),
    m_local(false),
    m_done(false)
{
    init();
}

task_copy::task_copy (task_list&        tl_,
                      MultiFab&         mf,
                      int               dgrid,
                      const MultiFab&   smf,
                      int               sgrid,
                      const Box&        bx,
                      const task_proxy& tp)
    :
    task(tl_),
#ifdef BL_USE_MPI
    m_request(MPI_REQUEST_NULL),
#endif
    tmp(0),
    m_mf(mf),
    m_smf(smf),
    m_dgrid(dgrid),
    m_sgrid(sgrid),
    m_bx(bx),
    m_sbx(bx),
    m_local(false),
    m_done(false)
{
    depend_on(tp);

    init();
}

inline
bool
mfeq (const MultiFab& a, const MultiFab& b)
{
    return &a == &b;
}

bool
task_copy::depends_on_q (const task* t1) const
{
    if (!work_to_do()) return false;

    if (!mfeq(m_mf,m_smf)) return false;

    if (const task_copy* t1tc = dynamic_cast<const task_copy*>(t1))
    {
        if (m_sgrid == t1tc->m_dgrid && m_sbx.intersects(t1tc->m_bx)) return true;
        if (m_dgrid == t1tc->m_dgrid && m_bx.intersects(t1tc->m_bx)) return true;
        if (m_sgrid == t1tc->m_sgrid && m_sbx.intersects(t1tc->m_sbx)) return true;
        if (m_dgrid == t1tc->m_sgrid && m_bx.intersects(t1tc->m_sbx)) return true;
    }

    return false;
}

void
task::_do_depend ()
{
    if (ParallelDescriptor::NProcs() == 1)
        //
        // In this case just doing tasks in order will be OK.
        //
        return;

    for (list<task::task_proxy>::const_iterator cit = m_task_list.begin(); cit != m_task_list.end(); ++cit)
    {
        BLassert(!(*cit).null());

        if (depends_on_q(**cit))
            dependencies.push_back(*cit);
    }
}
  
task_copy::~task_copy ()
{
    delete tmp;
#ifdef BL_USE_MPI
    BLassert(m_request == MPI_REQUEST_NULL);
#endif
}

bool
task_copy::work_to_do () const
{
    return (is_local(m_mf,m_dgrid) || is_local(m_smf,m_sgrid)) && !m_done;
}

bool
task_copy::startup ()
{
    m_started = true;

    bool result = true;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (is_local(m_mf,m_dgrid))
        {
            tmp = new FArrayBox(m_sbx, m_smf.nComp());
            int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, HG::mpi_comm, &m_request);
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BLassert(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_smf,m_sgrid)) 
        {
            tmp = new FArrayBox(m_sbx,m_smf.nComp());
            tmp->copy(m_smf[m_sgrid],m_sbx);
            int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_mf,  m_dgrid), m_sno, HG::mpi_comm, &m_request);
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BLassert(m_request != MPI_REQUEST_NULL);
        }
        else
        {
            result = false;
        }
    }
#endif

    return result;
}

bool
task_copy::ready ()
{
    BLassert(is_started());

    if (m_local)
    {
        BLassert(!m_done);
        m_mf[m_dgrid].copy(m_smf[m_sgrid],m_sbx,0,m_bx,0,m_mf.nComp());
        return true;
    }

#ifdef BL_USE_MPI
    int flag;
    MPI_Status status;
    BLassert(m_request != MPI_REQUEST_NULL);
    int res = MPI_Test(&m_request, &flag, &status);
    if (res != 0)
        ParallelDescriptor::Abort(res);
    if (flag)
    {
        BLassert(m_request == MPI_REQUEST_NULL);
        if (is_local(m_mf, m_dgrid))
        {
            int count;
            BLassert(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
            BLassert(status.MPI_TAG    == m_sno);
            if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
                ParallelDescriptor::Abort(res);
            BLassert(count == tmp->box().numPts()*tmp->nComp());
            m_mf[m_dgrid].copy(*tmp, tmp->box(), 0, m_bx, 0, m_smf.nComp());
        }
        return true;
    }
#endif

    return false;
}

void
task_copy::hint () const
{
    task::_hint();
    if (is_local(m_smf, m_sgrid) && is_local(m_mf, m_dgrid))
    {
        HG_DEBUG_OUT("L");
    }
    else if (is_local(m_smf, m_sgrid))
    {
        HG_DEBUG_OUT("S");
    }
    else if (is_local(m_mf, m_dgrid))
    {
        HG_DEBUG_OUT("R");
    }
    else
    {
        HG_DEBUG_OUT("?");
    }
    HG_DEBUG_OUT('(' << m_dgrid << "," << m_sgrid << ')' << m_sbx << ' ' << m_bx  << ' ' );
    HG_DEBUG_OUT(")" << endl);
}

task_copy_local::task_copy_local (task_list&      tl_,
                                  FArrayBox*      fab_,
                                  int             target_proc_id,
                                  const Box&      bx,
                                  const MultiFab& smf_,
                                  int             grid)
    :
    task(tl_),
    tmp(0),
    m_fab(fab_),
    m_smf(smf_),
    m_bx(bx),
    m_sgrid(grid),
    m_target_proc_id(target_proc_id),
#ifdef BL_USE_MPI
    m_request(MPI_REQUEST_NULL),
#endif
    m_local(false),
    m_done(false)
{
    if (work_to_do())
    {
        _do_depend();

        if (m_fab != 0 && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_fab->copy(m_smf[m_sgrid],m_bx);
                //
                // Flip the work_to_do() bit.
                //
                m_done = true;
            }
        }
    }
}

task_copy_local::~task_copy_local ()
{
    HG_DEBUG_OUT("task_copy_local::~task_copy_local(): delete tmp" << endl);
    delete tmp;
}

bool
task_copy_local::work_to_do () const
{
    return (m_fab != 0 || is_local(m_smf,m_sgrid)) && !m_done;
}

bool
task_copy_local::startup ()
{
    m_started = true;

    bool result = true;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (m_fab != 0)
        {
            tmp = new FArrayBox(m_bx, m_smf.nComp());
            int res = MPI_Irecv(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, processor_number(m_smf, m_sgrid), m_sno, HG::mpi_comm, &m_request);
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BLassert(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_smf, m_sgrid)) 
        {
            tmp = new FArrayBox(m_bx, m_smf.nComp());
            tmp->copy(m_smf[m_sgrid], m_bx);
            int res = MPI_Isend(tmp->dataPtr(), tmp->box().numPts()*tmp->nComp(), MPI_DOUBLE, m_target_proc_id, m_sno, HG::mpi_comm, &m_request);
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BLassert(m_request != MPI_REQUEST_NULL);
        }
        else
        {
            result = false;
        }
    }
#endif

    return result;
}

bool
task_copy_local::ready ()
{
    BLassert(is_started());

    if (m_local)
    {
        BLassert(!m_done);
        m_fab->copy(m_smf[m_sgrid], m_bx);
        return true;
    }
        
#ifdef BL_USE_MPI
    int flag;
    MPI_Status status;
    BLassert (m_request != MPI_REQUEST_NULL);
    int res = MPI_Test(&m_request, &flag, &status);
    if (res != 0)
        ParallelDescriptor::Abort(res);
    if (flag)
    {
        BLassert (m_request == MPI_REQUEST_NULL);
        if (m_fab)
        {
            int count;
            BLassert(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
            BLassert(status.MPI_TAG    == m_sno);
            int res = MPI_Get_count(&status, MPI_DOUBLE, &count);
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BLassert(count == tmp->box().numPts()*tmp->nComp());
            m_fab->copy(*tmp, m_bx);
        }
        return true;
    }
#endif

    return false;
}

void
task_copy_local::hint () const
{
    task::_hint();
    if (m_fab !=0 && is_local(m_smf, m_sgrid))
        HG_DEBUG_OUT("L");
    else if (m_fab != 0)
        HG_DEBUG_OUT("R");
    else if (is_local(m_smf, m_sgrid))
        HG_DEBUG_OUT("S");
    else
        HG_DEBUG_OUT("?");
    HG_DEBUG_OUT(m_bx <<  ' ' <<  m_sgrid << ' ');
    HG_DEBUG_OUT(")" << endl);
}

bool
task_copy_local::depends_on_q (const task* t1) const
{
    if (!work_to_do()) return false;

    if (const task_copy_local* t1tc = dynamic_cast<const task_copy_local*>(t1))
    {
        if (!mfeq(m_smf, t1tc->m_smf)) return false;
        if (m_sgrid == t1tc->m_sgrid && m_bx.intersects(t1tc->m_bx)) return true;
    }

    return false;
}

task_fab::task_fab (task_list&      tl_,
                    const MultiFab& t_,
                    int             tt_,
                    const Box&      region_,
                    int             ncomp_)
    :
    task(tl_),
    target(0),
    region(region_),
    ncomp(ncomp_),
    m_target_proc_id(processor_number(t_,tt_))
{
    BLassert(m_sno != 0);
    if (is_local(t_, tt_))
        target = new FArrayBox(region, ncomp);
}

task_fab::~task_fab ()
{
    delete target;
}

bool
task_fab::work_to_do () const
{
    return ParallelDescriptor::MyProc() == m_target_proc_id;
}

task_fec_base::task_fec_base (task_list& tl_,
                              MultiFab&  s_,
                              int        igrid_)
    :
    task(tl_),
    s(s_),
    igrid(igrid_),
    done(false)
{}

task_fec_base::~task_fec_base ()
{
    HG_DEBUG_OUT("task_fec_base::~task_fec_base()" << endl);
}

bool
task_fec_base::work_to_do () const
{
    return (is_local_target() || !tfvect.empty()) && !done;
}

void
task_fec_base::push_back (task_fab* tf)
{
    task_proxy tp = m_task_list.add_task(tf);

    if (!tp.null())
    {
        tfvect.push_back(tp);
        dependencies.push_back(tp);
    }
}
