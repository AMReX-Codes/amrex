#include <winstd.H>

#include <ParmParse.H>
#include <amr_defs.H>
#include <hgparallel.H>
#include <boundary.H>
#include <Profiler.H>


#ifdef HG_DEBUG
#include <algorithm>
#include <iomanip>
#include <typeinfo>
#include <cstdio>
#include <fstream>
#include <sstream>
#endif

bool HG_is_debugging       = false;
MPI_Comm HG::mpi_comm      = MPI_COMM_WORLD;
int HG::max_live_tasks     = 100;
int HG::multigrid_maxiter  = 100;
int HG::cgsolve_maxiter    = 250;
int HG::pverbose           = 0;
bool HG::initialized       = false;
double HG::cgsolve_tolfact = 1.0e-3;

#ifdef HG_DEBUG
std::ofstream debug_out;
#endif

void
HG::MPI_init ()
{
    if (!initialized)
    {
        ParmParse pp("HG");

        pp.query("cgsolve_maxiter", cgsolve_maxiter);
	pp.query("debugging", HG_is_debugging);
        pp.query("multigrid_maxiter", multigrid_maxiter);
        pp.query("cgsolve_tolfact", cgsolve_tolfact);
        pp.query("max_live_tasks", max_live_tasks);
        pp.query("pverbose", pverbose);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "HG.cgsolve_maxiter = "   << cgsolve_maxiter   << '\n';
            std::cout << "HG.multigrid_maxiter = " << multigrid_maxiter << '\n';
            std::cout << "HG.cgsolve_tolfact = "   << cgsolve_tolfact   << '\n';
            std::cout << "HG.max_live_tasks = "    << max_live_tasks    << '\n';
            std::cout << "HG.pverbose = "          << pverbose          << '\n';
        }

        initialized = true;

        ParallelDescriptor::Comm_dup(MPI_COMM_WORLD, mpi_comm);

#ifdef HG_DEBUG
        std::ostringstream fname;
        fname << "gu" << ParallelDescriptor::NProcs()
              << "_" << ParallelDescriptor::MyProc() << std::ends;
        debug_out.open(fname.str().c_str(), std::ios::trunc);
        if ( debug_out.fail() ) BoxLib::Error( "Failed to open debug file" );
        debug_out << std::setprecision(15);
#endif
    }
}

void
HG::MPI_finish ()
{
#ifdef HG_DEBUG
    debug_out.close();
#endif
}


// task
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
task::need_to_communicate (int& /*with*/) const
{
    return false;
}

void
task::print_dependencies (std::ostream& os) const
{
    os << "Task " << get_sequence_number() << " depends on ( ";

    for (std::list<task_proxy>::const_iterator lit = dependencies.begin();
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
    std::list<task_proxy>::iterator lit = dependencies.begin();

    while (lit != dependencies.end())
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
    HG_DEBUG_OUT(")" << std::endl);
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

void
task::_do_depend ()
{
    if (ParallelDescriptor::NProcs() == 1)
        //
        // In this case just doing tasks in order will be OK.
        //
        return;

    for (std::list<task::task_proxy>::const_iterator cit = m_task_list.begin();
	 cit != m_task_list.end();
         ++cit)
    {
        BL_ASSERT(!(*cit).null());

        if (depends_on_q(**cit))
            dependencies.push_back(*cit);
    }
}


// task_list
task_list::task_list ()
    :
    seq_no(ParallelDescriptor::NProcs(), 1),
    verbose(HG::pverbose != 0)
{}

task_list::~task_list () {}

task::task_proxy
task_list::add_task (task* t)
{
    BL_ASSERT(t != 0);

    if (t->is_finished())
    {
        delete t;
        return task::task_proxy(0);
    }
    else
    {
        tasks.push_back(task::task_proxy(t));
        return tasks.back();
    }
}

void
task_list::print_dependencies (std::ostream& os) const
{
    os << "Task list ( " << '\n';

    for (std::list<task::task_proxy>::const_iterator tli = tasks.begin();
	 tli != tasks.end();
         ++tli)
    {
        (*tli)->print_dependencies(os);
        os << '\n';
    }

    os << ')' << std::endl;
}

void
task_list::execute (const char* msg)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::execute()");

    if (HG_is_debugging)
        ParallelDescriptor::Barrier(HG::mpi_comm);
    //
    // Assign message tag IDs ...
    //
    std::list<task::task_proxy>::iterator tli = tasks.begin();

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
        HG_DEBUG_OUT(HG::mpi_comm);
        HG_DEBUG_OUT(" with " << tasks.size() << " elements " << std::endl);
        print_dependencies(debug_out);
#endif
    }

    int live_tasks = 0;
    long total_snd = 0, total_rcv = 0, maxpacket = 0;

restart:
    while (!tasks.empty())
    {
        std::list<task::task_proxy>::iterator tli = tasks.begin();

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

                    total_snd += sndcnt;
                    total_rcv += rcvcnt;
                    maxpacket  = std::max(maxpacket,std::max(sndcnt,rcvcnt));
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
    if (verbose && maxpacket > 0)
    {
        char buf[512];

        sprintf(buf,
                "CPU(%d) %s: Sent: %ld Rcvd: %ld MaxPkt: %ld",
                ParallelDescriptor::MyProc(),
                msg,
                total_snd,
                total_rcv,
                maxpacket);

        std::cout << buf << std::endl;
    }
#endif
}

inline
bool
mfeq (const MultiFab& a, const MultiFab& b)
{
    return &a == &b;
}

// task_copy_base:
task_copy_base::task_copy_base(task_list&      tl_,
			       MultiFab&       dmf_,
			       int             dgrid_,
			       const Box&      db_,
			       const MultiFab& smf_,
			       int             sgrid_,
			       const Box&      sb_)
  :
  task(tl_),
  m_request(MPI_REQUEST_NULL),
  m_tmp(0),
  m_dmf(dmf_),
  m_dgrid(dgrid_),
  m_dbx(db_),
  m_smf(smf_),
  m_sgrid(sgrid_),
  m_sbx(sb_),
  m_local(false)
{}

task_copy_base::~task_copy_base ()
{
    delete m_tmp;
}

bool
task_copy_base::depends_on_q (const task* t1) const
{
    if (!mfeq(m_dmf, m_smf)) return false;

    if (const task_copy_base* t1tc = dynamic_cast<const task_copy_base*>(t1))
    {
        if (m_sgrid == t1tc->m_dgrid && m_sbx.intersects(t1tc->m_dbx)) return true;
        if (m_dgrid == t1tc->m_dgrid && m_dbx.intersects(t1tc->m_dbx)) return true;
        if (m_sgrid == t1tc->m_sgrid && m_sbx.intersects(t1tc->m_sbx)) return true;
        if (m_dgrid == t1tc->m_sgrid && m_dbx.intersects(t1tc->m_sbx)) return true;
    }

    return false;
}

bool
task_copy_base::need_to_communicate (int& with) const
{
    bool result = false;

    if (!m_local)
    {
        if (is_local(m_dmf, m_dgrid))
        {
            with   = processor_number(m_smf, m_sgrid);
            result = true;
        }
        else if (is_local(m_smf, m_sgrid))
        {
            with   = processor_number(m_dmf, m_dgrid);
            result = true;
        }
    }

    return result;
}

bool
task_copy_base::startup (long& sndcnt, long& rcvcnt)
{
    m_started = true;

    bool result = true;

    if (!m_local)
    {
        if (is_local(m_dmf, m_dgrid))
        {
            m_tmp = new FArrayBox(m_sbx, m_smf.nComp());
            rcvcnt = m_tmp->box().numPts()*m_tmp->nComp();

            m_request = ParallelDescriptor::Arecv(m_tmp->dataPtr(),
                                                  rcvcnt,
                                                  processor_number(m_smf,m_sgrid),
                                                  m_sno,
                                                  HG::mpi_comm).req();
            rcvcnt *= sizeof(double);

            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_smf, m_sgrid))
        {
            m_tmp = new FArrayBox(m_sbx, m_smf.nComp());
            m_tmp->copy(m_smf[m_sgrid], m_sbx);
            sndcnt = m_tmp->box().numPts()*m_tmp->nComp();

            m_request = ParallelDescriptor::Asend(m_tmp->dataPtr(),
                                                  sndcnt,
                                                  processor_number(m_dmf,m_dgrid),
                                                  m_sno,
                                                  HG::mpi_comm).req();
            sndcnt *= sizeof(double);

            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else
        {
            result = false;
        }
    }

    return result;
}

void
task_copy_base::hint () const
{
    task::_hint();
    if (is_local(m_smf, m_sgrid) && is_local(m_dmf, m_dgrid))
    {
        HG_DEBUG_OUT("L");
    }
    else if (is_local(m_smf, m_sgrid))
    {
        HG_DEBUG_OUT("S");
    }
    else if (is_local(m_dmf, m_dgrid))
    {
        HG_DEBUG_OUT("R");
    }
    else
    {
        HG_DEBUG_OUT("?");
    }
    HG_DEBUG_OUT('(' << m_dgrid << "," << m_sgrid << ')'
		 << m_dbx << ' ' << m_sbx  << ' ' );
    HG_DEBUG_OUT(")" << std::endl);
}

// task_copy
task_copy::task_copy (task_list&      tl_,
                      MultiFab&       mf,
                      int             dgrid,
                      const MultiFab& smf,
                      int             sgrid,
                      const Box&      bx)
    :
    task_copy_base(tl_, mf, dgrid, bx, smf, sgrid, bx)
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
    task_copy_base(tl_, mf, dgrid, db, smf, sgrid, sb)
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
    task_copy_base(tl_, mf, dgrid, bx, smf, sgrid, bx)
{
    depend_on(tp);
    init();
}

void
task_copy::init ()
{
    if (is_local(m_dmf, m_dgrid) || is_local(m_smf, m_sgrid))
    {
        _do_depend();

        if (is_local(m_dmf, m_dgrid) && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_dmf[m_dgrid].copy(m_smf[m_sgrid], m_sbx, 0, m_dbx, 0, m_dmf.nComp());

                m_finished = true;
            }
        }
    }
    else
    {
        m_finished = true;
    }
}

bool
task_copy::ready ()
{
    BL_ASSERT(is_started());

    if (m_local)
    {
        BL_ASSERT(!m_finished);
        m_dmf[m_dgrid].copy(m_smf[m_sgrid],m_sbx,0,m_dbx,0,m_dmf.nComp());
        return true;
    }

    int flag;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);

    ParallelDescriptor::Test(m_request, flag, status);

    if (flag)
    {
        BL_ASSERT(m_request == MPI_REQUEST_NULL);
        if (is_local(m_dmf, m_dgrid))
        {
#if defined(BL_USE_MPI) && !defined(NDEBUG)
            int count,res;
            BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
            BL_ASSERT(status.MPI_TAG    == m_sno);
            if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(count == m_tmp->box().numPts()*m_tmp->nComp());
#endif
            m_dmf[m_dgrid].copy(*m_tmp,m_tmp->box(),0,m_dbx,0,m_smf.nComp());
        }
        return true;
    }

    return false;
}


// task_local_base
task_local_base::task_local_base (task_list&      tl_,
                                  FArrayBox*      fab_,
                                  int             target_proc_id,
                                  const Box&      bx,
                                  const Box&      region,
                                  const MultiFab& smf_,
                                  int             grid)
    :
    task(tl_),
    m_request(MPI_REQUEST_NULL),
    m_tmp(0),
    m_fab(fab_),
    m_smf(smf_),
    m_bx(bx),
    m_region(region),
    m_sgrid(grid),
    m_target_proc_id(target_proc_id),
    m_local(false)
{}

task_local_base::~task_local_base ()
{
    delete m_tmp;
}

bool
task_local_base::need_to_communicate (int& with) const
{
    bool result = false;

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

    return result;
}

bool
task_local_base::startup (long& sndcnt, long& rcvcnt)
{
    m_started = true;

    bool result = true;

    if (!m_local)
    {
        if (m_fab != 0)
        {
            m_tmp = new FArrayBox(m_bx, m_smf.nComp());
            rcvcnt = m_tmp->box().numPts()*m_tmp->nComp();

            m_request = ParallelDescriptor::Arecv(m_tmp->dataPtr(),
                                                  rcvcnt,
                                                  processor_number(m_smf,m_sgrid),
                                                  m_sno,
                                                  HG::mpi_comm).req();
            rcvcnt *= sizeof(double);

            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else if (is_local(m_smf, m_sgrid))
        {
            m_tmp = new FArrayBox(m_bx, m_smf.nComp());
            m_tmp->copy(m_smf[m_sgrid], m_bx);
            sndcnt = m_tmp->box().numPts()*m_tmp->nComp();

            m_request = ParallelDescriptor::Asend(m_tmp->dataPtr(),
                                                  sndcnt,
                                                  m_target_proc_id,
                                                  m_sno,
                                                  HG::mpi_comm).req();
            sndcnt *= sizeof(double);

            BL_ASSERT(m_request != MPI_REQUEST_NULL);
        }
        else
        {
            result = false;
        }
    }

    return result;
}

bool
task_local_base::depends_on_q (const task* t1) const
{
    if (const task_local_base* t1tc = dynamic_cast<const task_local_base*>(t1))
    {
        if (!(m_fab == t1tc->m_fab))     return false;
        if (!mfeq(m_smf, t1tc->m_smf))   return false;
        if (m_region.intersects(t1tc->m_region)) return true;
        if (m_sgrid != t1tc->m_sgrid && m_bx.intersects(t1tc->m_bx))
            return true;
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
    HG_DEBUG_OUT(")" << std::endl);
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
    if (m_fab != 0 || is_local(m_smf, m_sgrid))
    {
        _do_depend();

        if (m_fab != 0 && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_fab->copy(m_smf[m_sgrid], m_bx);

                m_finished = true;
            }
        }
    }
    else
    {
        m_finished = true;
    }
}

bool
task_copy_local::ready ()
{
    BL_ASSERT(is_started());

    if (m_local)
    {
        BL_ASSERT(!m_finished);
        m_fab->copy(m_smf[m_sgrid], m_bx);
        return true;
    }

    int flag;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);

    ParallelDescriptor::Test(m_request, flag, status);

    if (flag)
    {
        BL_ASSERT(m_request == MPI_REQUEST_NULL);
        if (m_fab)
        {
#if defined(BL_USE_MPI) && !defined(NDEBUG)
            int count,res;
            BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
            BL_ASSERT(status.MPI_TAG    == m_sno);
            if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
                ParallelDescriptor::Abort(res);
            BL_ASSERT(count == m_tmp->box().numPts()*m_tmp->nComp());
#endif
            m_fab->copy(*m_tmp, m_bx);
        }
        return true;
    }

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
    m_bx = BoxLib::grow(tmpb, 4);

    if (m_fab != 0 || is_local(m_smf, m_sgrid))
    {
        _do_depend();

        if (m_fab != 0 && is_local(m_smf, m_sgrid))
        {
            m_local = true;

            if (dependencies.empty())
            {
                m_bdy->fill(*m_fab, m_region, m_smf[m_sgrid], m_domain);

                m_finished = true;
            }
        }
    }
    else
    {
        m_finished = true;
    }
}

bool
task_bdy_fill::ready ()
{
    BL_ASSERT(is_started());

    if (m_local)
    {
        BL_ASSERT(!m_finished);
	m_bdy->fill(*m_fab, m_region, m_smf[m_sgrid], m_domain);
	return true;
    }

    int flag;
    MPI_Status status;
    BL_ASSERT(m_request != MPI_REQUEST_NULL);

    ParallelDescriptor::Test(m_request, flag, status);

    if (flag)
    {
	BL_ASSERT(m_request == MPI_REQUEST_NULL);
	if (m_fab)
	{
#if defined(BL_USE_MPI) && !defined(NDEBUG)
	    int count,res;
	    BL_ASSERT(status.MPI_SOURCE == processor_number(m_smf, m_sgrid));
	    BL_ASSERT(status.MPI_TAG    == m_sno);
	    if ((res = MPI_Get_count(&status, MPI_DOUBLE, &count)) != 0)
		ParallelDescriptor::Abort(res);
	    BL_ASSERT(count == m_tmp->box().numPts()*m_tmp->nComp());
#endif
	    m_bdy->fill(*m_fab, m_region, *m_tmp, m_domain);
	}
	return true;
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
    m_target_proc_id(processor_number(t_, tt_))
{
    BL_ASSERT(m_sno > 0);

    if (is_local(t_, tt_))
    {
        target = new FArrayBox(region, ncomp);
    }
    else
    {
        m_finished = true;
    }
}

task_fab::~task_fab ()
{
    delete target;
}

task_fec_base::task_fec_base (task_list& tl_,
                              MultiFab&  s_,
                              int        igrid_)
    :
    task(tl_),
    s(s_),
    igrid(igrid_)
{
    if (!is_local_target()) m_finished = true;
}

task_fec_base::~task_fec_base () {}

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
