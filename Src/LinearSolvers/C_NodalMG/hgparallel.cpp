//BL_COPYRIGHT_NOTICE

#include <ParmParse.H>

#include <amr_defs.H>
#include <hgparallel.H>
#include <boundary.H>


#ifdef HG_DEBUG
#ifdef BL_USE_NEW_HFILES
#include <typeinfo>
#include <cstdio>
#else
#include <typeinfo.h>
#include <stdio.h>
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

void
HG::MPI_init ()
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
            cout << "HG.cgsolve_maxiter = "   << cgsolve_maxiter   << '\n';
            cout << "HG.multigrid_maxiter = " << multigrid_maxiter << '\n';
            cout << "HG.cgsolve_tolfact = "   << cgsolve_tolfact   << '\n';
            cout << "HG.max_live_tasks = "    << max_live_tasks    << '\n';
            cout << "HG.pverbose = "          << pverbose          << '\n';
        }

        initialized = true;

#ifdef BL_USE_MPI
        int res = MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
        if (res != 0)
            ParallelDescriptor::Abort(res);
#endif
    }
}

void HG::MPI_finish () {}

task::task (task_list& tl_)
    :
    m_task_list(tl_),
    m_sno(1),
    m_cnt(1),
    m_finished(false),
    m_started(false)
{
    BL_ASSERT(m_sno > 0);
}

task::~task ()
{
    BL_ASSERT(m_sno > 0);
}

bool
task::startup (long&, long&)
{
    return m_started = true;
}

bool
task::work_to_do () const
{
    return true;
}

bool
task::need_to_communicate (int& /*with*/) const
{
    return false;
}

void
task::print_dependencies (ostream& os) const
{
    os << "Task " << get_sequence_number() << " depends on ( ";
    for (list<task_proxy>::const_iterator lit = dependencies.begin();
	 lit != dependencies.end();
	 ++lit)
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
    for (list<task_proxy>::iterator lit = dependencies.begin();
	 lit != dependencies.end();)
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
    BL_ASSERT(m_sno > 0);
    HG_DEBUG_OUT("(" << typeid(*this).name() << ' '
		 << m_sno << ' ' << m_started << ' ');
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
    BL_ASSERT(is_started());
    return true;
}

task_list::task_list ()
    :
    seq_no(ParallelDescriptor::NProcs(), 1),
    verbose(HG::pverbose)
{}

task_list::~task_list ()
{
}

task::task_proxy
task_list::add_task (task* t)
{
    BL_ASSERT(t != 0);

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
    os << "Task list ( " << '\n';
    for (list<task::task_proxy>::const_iterator tli = tasks.begin();
	 tli != tasks.end(); ++tli)
    {
        (*tli)->print_dependencies(os);
        os << '\n';
    }
    os << ")" << endl;
}

void
task_list::execute (const char* msg)
{
#ifdef BL_USE_MPI
    if (HG_is_debugging)
        MPI_Barrier(HG::mpi_comm);
#endif
    //
    // Assign message tag IDs ...
    //
    list<task::task_proxy>::iterator tli = tasks.begin();

    for ( ; tli != tasks.end(); ++tli)
    {
        int with = -1;

        if ((*tli)->need_to_communicate(with))
        {
            (*tli)->set_sequence_number(get_then_advance(with));
        }
    }

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

    int live_tasks    = 0;
    long total_sndcnt = 0, total_rcvcnt = 0, maxpacketsize = 0;

restart:
    while (!tasks.empty())
    {
        list<task::task_proxy>::iterator tli = tasks.begin();

        while (tli != tasks.end())
        {
            task::task_proxy t = *tli;

            BL_ASSERT(!t.null());

            if (verbose) t->hint();

            if (t->depend_ready())
            {
                if (!t->is_started())
                {
                    if (live_tasks > HG::max_live_tasks) goto restart;

                    long sndcnt = 0, rcvcnt = 0;

                    if (!t->startup(sndcnt, rcvcnt))
                    {
                        t.set_finished();
                        tasks.erase(tli++);
                        continue;
                    }

                    total_sndcnt += sndcnt;
                    total_rcvcnt += rcvcnt;
                    maxpacketsize = Max(maxpacketsize, Max(sndcnt, rcvcnt));
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
    BL_ASSERT(live_tasks == 0);

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        seq_no[i] = 1;

#ifdef HG_DEBUG
    if (verbose && maxpacketsize > 0)
    {
        char buf[512];

        sprintf(buf,
                "CPU(%d) %s: Sent: %ld Rcvd: %ld MaxPkt: %ld",
                ParallelDescriptor::MyProc(),
                msg,
                total_sndcnt,
                total_rcvcnt,
                maxpacketsize);

        cout << buf << endl;
    }
#endif
}

void
task_copy::init ()
{
    if (work_to_do())
    {
        _do_depend();

        if (is_local(m_mf, m_dgrid) && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_mf[m_dgrid].copy(m_smf[m_sgrid], m_sbx, 0, m_bx, 0, m_mf.nComp());
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

    if (!mfeq(m_mf, m_smf)) return false;

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

    for (list<task::task_proxy>::const_iterator cit = m_task_list.begin();
	 cit != m_task_list.end(); ++cit)
    {
        BL_ASSERT(!(*cit).null());

        if (depends_on_q(**cit))
            dependencies.push_back(*cit);
    }
}

task_copy::~task_copy ()
{
    delete tmp;
#ifdef BL_USE_MPI
    BL_ASSERT(m_request == MPI_REQUEST_NULL);
#endif
}

bool
task_copy::work_to_do () const
{
    return (is_local(m_mf, m_dgrid) || is_local(m_smf, m_sgrid)) && !m_done;
}

bool
task_copy::need_to_communicate (int& with) const
{
    bool result = false;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (is_local(m_mf, m_dgrid))
        {
            with   = processor_number(m_smf, m_sgrid);
            result = true;
        }
        else if (is_local(m_smf, m_sgrid))
        {
            with   = processor_number(m_mf, m_dgrid);
            result = true;
        }
    }
#endif

    return result;
}

bool
task_copy::startup (long& sndcnt, long& rcvcnt)
{
    m_started = true;

    bool result = true;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (is_local(m_mf, m_dgrid))
        {
            static RunStats irecv_stats("hg_irecv");
            tmp = new FArrayBox(m_sbx, m_smf.nComp());
            irecv_stats.start();
            rcvcnt = tmp->box().numPts()*tmp->nComp();
            int res = MPI_Irecv(tmp->dataPtr(),
                                rcvcnt,
                                MPI_DOUBLE,
                                processor_number(m_smf, m_sgrid),
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            rcvcnt *= sizeof(double);
            irecv_stats.end();

            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_smf, m_sgrid))
        {
            static RunStats isend_stats("hg_isend");
            tmp = new FArrayBox(m_sbx, m_smf.nComp());
            tmp->copy(m_smf[m_sgrid], m_sbx);
            isend_stats.start();
            sndcnt = tmp->box().numPts()*tmp->nComp();
            int res = MPI_Isend(tmp->dataPtr(),
                                sndcnt,
                                MPI_DOUBLE,
                                processor_number(m_mf, m_dgrid),
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            sndcnt *= sizeof(double);
            isend_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
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
    BL_ASSERT(is_started());

    if (m_local)
    {
        BL_ASSERT(!m_done);
        m_mf[m_dgrid].copy(m_smf[m_sgrid], m_sbx, 0, m_bx, 0, m_mf.nComp());
        return true;
    }

#ifdef BL_USE_MPI
    int flag, res;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);
    static RunStats test_stats("hg_test");
    test_stats.start();
    if ((res = MPI_Test(&m_request, &flag, &status)) != 0)
        ParallelDescriptor::Abort(res);
    test_stats.end();
    if (flag)
    {
        BL_ASSERT(m_request == MPI_REQUEST_NULL);
        if (is_local(m_mf, m_dgrid))
        {
#ifndef NDEBUG
            int count;
            BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
            BL_ASSERT(status.MPI_TAG    == m_sno);
            if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(count == tmp->box().numPts()*tmp->nComp());
#endif
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
    HG_DEBUG_OUT('(' << m_dgrid << "," << m_sgrid << ')'
		 << m_sbx << ' ' << m_bx  << ' ' );
    HG_DEBUG_OUT(")" << endl);
}

task_local_base::task_local_base (task_list&      tl_,
                                  FArrayBox*      fab_,
                                  int             target_proc_id,
                                  const Box&      bx,
                                  const Box&      region,
                                  const MultiFab& smf_,
                                  int             grid)
    :
    task(tl_),
#ifdef BL_USE_MPI
    m_request(MPI_REQUEST_NULL),
#endif
    tmp(0),
    m_fab(fab_),
    m_smf(smf_),
    m_bx(bx),
    m_region(region),
    m_sgrid(grid),
    m_target_proc_id(target_proc_id),
    m_local(false),
    m_done(false)
{
}

task_local_base::~task_local_base ()
{
    delete tmp;
}

bool
task_local_base::work_to_do () const
{
    return (m_fab != 0 || is_local(m_smf, m_sgrid)) && !m_done;
}

bool
task_local_base::need_to_communicate (int& with) const
{
    bool result = false;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (m_fab != 0)
        {
            with   = processor_number(m_smf, m_sgrid);
            result = true;
        }
        else if (is_local(m_smf,  m_sgrid))
        {
            with   = m_target_proc_id;
            result = true;
        }
    }
#endif

    return result;
}

bool
task_local_base::startup (long& sndcnt, long& rcvcnt)
{
    m_started = true;

    bool result = true;

#ifdef BL_USE_MPI
    if (!m_local)
    {
        if (m_fab != 0)
        {
            static RunStats irecv_stats("hg_irecv");
            tmp = new FArrayBox(m_bx, m_smf.nComp());
            irecv_stats.start();
            rcvcnt = tmp->box().numPts()*tmp->nComp();
            int res = MPI_Irecv(tmp->dataPtr(),
                                rcvcnt,
                                MPI_DOUBLE,
                                processor_number(m_smf, m_sgrid),
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            rcvcnt *= sizeof(double);
            irecv_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_smf, m_sgrid))
        {
            static RunStats isend_stats("hg_isend");
            tmp = new FArrayBox(m_bx, m_smf.nComp());
            tmp->copy(m_smf[m_sgrid], m_bx);
            isend_stats.start();
            sndcnt = tmp->box().numPts()*tmp->nComp();
            int res = MPI_Isend(tmp->dataPtr(),
                                sndcnt,
                                MPI_DOUBLE,
                                m_target_proc_id,
                                m_sno,
                                HG::mpi_comm,
                                &m_request);
            sndcnt *= sizeof(double);
            isend_stats.end();
            if (res != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(m_request != MPI_REQUEST_NULL);
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
task_local_base::depends_on_q (const task* t1) const
{
    if (!work_to_do()) return false;

    if (const task_local_base* t1tc = dynamic_cast<const task_local_base*>(t1))
    {
        if (!(m_fab == t1tc->m_fab)) return false;
        if (!mfeq(m_smf, t1tc->m_smf)) return false;
        if (m_region.intersects(t1tc->m_region)) return true;
    }

    return false;
}

void
task_local_base::hint () const
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

task_copy_local::task_copy_local (task_list&      tl_,
                                  FArrayBox*      fab_,
                                  int             target_proc_id,
                                  const Box&      bx,
                                  const MultiFab& smf_,
                                  int             grid)
    :
    task_local_base(tl_,fab_,target_proc_id,bx,bx,smf_,grid)
{
    if (work_to_do())
    {
        _do_depend();

        if (m_fab != 0 && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_fab->copy(m_smf[m_sgrid], m_bx);
                //
                // Flip the work_to_do() bit.
                //
                m_done = true;
            }
        }
    }
}

bool
task_copy_local::ready ()
{
    BL_ASSERT(is_started());

    if (m_local)
    {
        BL_ASSERT(!m_done);
        m_fab->copy(m_smf[m_sgrid], m_bx);
        return true;
    }

#ifdef BL_USE_MPI
    int flag, res;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);
    static RunStats test_stats("hg_test");
    test_stats.start();
    if ((res = MPI_Test(&m_request, &flag, &status)) != 0)
        ParallelDescriptor::Abort(res);
    test_stats.end();
    if (flag)
    {
        BL_ASSERT(m_request == MPI_REQUEST_NULL);
        if (m_fab)
        {
#ifndef NDEBUG
            int count;
            BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
            BL_ASSERT(status.MPI_TAG    == m_sno);
            if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(count == tmp->box().numPts()*tmp->nComp());
#endif
            m_fab->copy(*tmp, m_bx);
        }
        return true;
    }
#endif

    return false;
}

task_bdy_fill::task_bdy_fill (task_list&          tl_,
                              const amr_boundary* bdy_,
                              FArrayBox*          fab_,
                              int                 target_proc_id,
                              const Box&          region_,
                              const MultiFab&     src_,
                              int                 grid_,
                              const Box&          domain_)
    :
    task_local_base(tl_,fab_,target_proc_id,Box(),region_,src_,grid_),
    m_domain(domain_),
    m_bdy(bdy_)
{
    BL_ASSERT(m_bdy != 0);

    Box tmpb = m_bdy->anImage(m_region,m_smf.boxArray()[m_sgrid],m_domain);

    m_bx = tmpb;
    //
    // This is a GROSS hack FIXME:  the growth should probably be set
    // by the refinement between the coarse/fine domains.
    //
    m_bx = ::grow(tmpb, 4);

    if (work_to_do())
    {
        _do_depend();

        if (m_fab != 0 && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_bdy->fill(*m_fab, m_region, m_smf[m_sgrid], m_domain);
                //
                // Flip the work_to_do() bit.
                //
                m_done = true;
            }
        }
    }
}

bool
task_bdy_fill::ready ()
{
    BL_ASSERT(is_started());

    if (m_local)
    {
        BL_ASSERT(!m_done);
	m_bdy->fill(*m_fab, m_region, m_smf[m_sgrid], m_domain);
	return true;
    }

#ifdef BL_USE_MPI
    int flag, res;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);
    static RunStats test_stats("hg_test");
    test_stats.start();
    if ((res = MPI_Test(&m_request, &flag, &status)) != 0)
	ParallelDescriptor::Abort(res);
    test_stats.end();
    if (flag)
    {
	BL_ASSERT(m_request == MPI_REQUEST_NULL);
	if (m_fab)
	{
#ifndef NDEBUG
	    int count;
	    BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
	    BL_ASSERT(status.MPI_TAG    == m_sno);
	    if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
		ParallelDescriptor::Abort(res);
	    BL_ASSERT(count == tmp->box().numPts()*tmp->nComp());
#endif
	    m_bdy->fill(*m_fab, m_region, *tmp, m_domain);
	}
	return true;
    }
#endif

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
    m_target_proc_id(processor_number(t_, tt_))
{
    BL_ASSERT(m_sno > 0);
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
    // HG_DEBUG_OUT("task_fec_base::~task_fec_base()" << endl);
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
