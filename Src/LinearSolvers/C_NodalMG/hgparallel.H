
#ifndef _HGPARALLEL_H_
#define _HGPARALLEL_H_

#ifdef BL_USE_NEW_HFILES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <MultiFab.H>

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

inline
int
processor_number (const MultiFab& r,
                  int             igrid)
{
    return r.DistributionMap()[igrid];
}

inline
bool
is_remote (const MultiFab& r,
           int             igrid)
{
    return ParallelDescriptor::MyProc() != processor_number(r, igrid);
}

inline
bool
is_local (const MultiFab& r,
          int             igrid)
{
    return ParallelDescriptor::MyProc() == processor_number(r, igrid);
}

struct HG
{
    static void MPI_init ();
    static void MPI_finish ();
#ifdef BL_USE_MPI
    static MPI_Comm mpi_comm;
#endif
    static int    max_live_tasks;
    static int    multigrid_maxiter;
    static int    cgsolve_maxiter;
    static double cgsolve_tolfact;
    static int    pverbose;
private:
    static bool   initialized;
};

#ifdef BL_USE_NEW_HFILES
#include <list>
using namespace std;
#else
#include <list.h>
#endif

class task_list;

class task
{
public:

    class task_proxy
    {
    public:
	explicit task_proxy (task* t_ = 0);
	task_proxy (const task_proxy& r);
	~task_proxy ();
	task_proxy& operator= (const task_proxy& r);
	task* operator-> () const;
	task* operator* () const;
	task* get () const;
	void set_finished ();
	bool is_finished () const;
	bool null () const;
    private:
        //
        // The data.
        //
        task* m_t;
    };

    typedef unsigned int sequence_number;
    explicit task (task_list& task_list_);
    virtual ~task ();
    virtual bool ready ();
    virtual bool startup (long& sdncnt, long& rcvcnt);
    bool is_started () const ;
    virtual bool depends_on_q (const task* t1) const;
    bool is_finished () const;
    virtual bool need_to_communicate (int& with) const;
    void depend_on (const task_proxy& t1);
    bool depend_ready ();
    virtual void hint () const;
    void print_dependencies (ostream& os) const;
    sequence_number get_sequence_number () const;
    void set_sequence_number (sequence_number sno);

protected:
    friend class task_proxy;
    void _do_depend ();
    void _hint () const;
    //
    // The data.
    //
    task_list&       m_task_list;
    list<task_proxy> dependencies;
    sequence_number  m_sno;
    unsigned int     m_cnt;
    bool             m_finished;
    bool             m_started;

private:
    //
    // Not defined.
    //
    task (const task&);
    task& operator= (const task&);
};

inline
task::task_proxy::task_proxy (task* t_)
{
    m_t = t_;
}

inline
task::task_proxy::task_proxy (const task_proxy& r)
{
    if ((m_t = r.m_t)) m_t->m_cnt++;
}

inline
task::task_proxy::~task_proxy ()
{
    if (m_t && --m_t->m_cnt == 0) delete m_t;
}

inline
task::task_proxy&
task::task_proxy::operator= (const task_proxy& r)
{
    if (this != &r)
    {
	if (m_t && --m_t->m_cnt == 0) delete m_t;

	if ((m_t = r.m_t)) m_t->m_cnt++;
    }
    return *this;
}

inline
task*
task::task_proxy::operator-> () const
{
    BL_ASSERT(!null());
    return m_t;
}

inline
task*
task::task_proxy::operator* () const
{
    BL_ASSERT(!null());
    return m_t;
}

inline
task*
task::task_proxy::get () const
{
    BL_ASSERT(!null());
    return m_t;
}

inline
void
task::task_proxy::set_finished ()
{
    BL_ASSERT(!null());
    m_t->m_finished = true;
}

inline
bool
task::task_proxy::is_finished () const
{
    BL_ASSERT(!null());
    return m_t->m_finished;
}

inline
bool
task::task_proxy::null () const
{
    return m_t == 0;
}

inline
bool
task::is_started () const
{
    return m_started;
}

inline
bool
task::is_finished () const
{
    return m_finished;
}

inline
task::sequence_number
task::get_sequence_number () const
{
    return m_sno;
}

inline
void
task::set_sequence_number (task::sequence_number sno)
{
    m_sno = sno;
}

class task_list
{
public:
    task_list ();
    ~task_list ();
    task::task_proxy add_task (task* t);
    //
    // Executes once through the task list, return true if any elements left.
    //
    void execute (const char* msg);

    list<task::task_proxy>::const_iterator begin () const;
    list<task::task_proxy>::iterator begin ();
    list<task::task_proxy>::const_iterator end () const;
    list<task::task_proxy>::iterator end ();

    bool empty () const;
    int size () const;

    task::sequence_number get_then_advance (int proc);
    void print_dependencies (ostream& os) const;

private:
    //
    // The data.
    //
    list<task::task_proxy>       tasks;
    Array<task::sequence_number> seq_no;
    bool                         verbose;
    //
    // Some class-wide static data.
    //
    static bool def_verbose;
};

inline
list<task::task_proxy>::const_iterator
task_list::begin () const
{
    return tasks.begin();
}

inline
list<task::task_proxy>::iterator
task_list::begin ()
{
    return tasks.begin();
}

inline
list<task::task_proxy>::const_iterator
task_list::end () const
{
    return tasks.end();
}

inline
list<task::task_proxy>::iterator
task_list::end ()
{
    return tasks.end();
}

inline
bool
task_list::empty () const
{
    return tasks.empty();
}

inline
int
task_list::size () const
{
    return tasks.size();
}

inline
task::sequence_number
task_list::get_then_advance (int proc)
{
    task::sequence_number sno = seq_no[proc];
    seq_no[proc] = sno + 1;
    return sno;
}

class task_copy_base
  :
    public task
{
public:

    task_copy_base(task_list&      tl_,
                   MultiFab&       dmf_,
                   int             dgrid_,
                   const Box&      db_,
                   const MultiFab& smf_,
                   int             sgrid_,
                   const Box&      sb_);

    virtual ~task_copy_base();
    virtual bool depends_on_q (const task* t) const;
    virtual bool startup (long&,long&);
    virtual bool need_to_communicate (int& with) const;
    virtual void hint () const;

protected:
#ifdef BL_USE_MPI
  MPI_Request     m_request;
#endif
  FArrayBox*      m_tmp;
  MultiFab&       m_dmf;
  const int       m_dgrid;
  const Box       m_dbx;
  const MultiFab& m_smf;
  const int       m_sgrid;
  const Box       m_sbx;
  bool            m_local;
};

class task_copy
    :
    public task_copy_base
{
public:

    task_copy (task_list&      tl_,
               MultiFab&       mf,
               int             dgrid,
               const MultiFab& smf,
               int             sgrid,
               const Box&      bx);

    task_copy (task_list&      tl_,
               MultiFab&       mf,
               int             dgrid,
               const Box&      db,
               const MultiFab& smf,
               int             sgrid,
               const Box&      sb);

    task_copy (task_list&        tl_,
               MultiFab&         mf,
               int               dgrid,
               const MultiFab&   smf,
               int               sgrid,
               const Box&        bx,
               const task_proxy& tp);

    virtual bool ready ();

private:
    //
    // Common function called by constructors.
    //
    void init ();
};

class task_local_base
    :
    public task
{
public:

    task_local_base (task_list&      tl_,
                     FArrayBox*      fab_,
                     int             target_proc_id,
                     const Box&      bx,
                     const Box&      region,
                     const MultiFab& mf_,
                     int             grid_);

    virtual ~task_local_base ();
    virtual bool startup (long&, long&);
    virtual bool depends_on_q (const task* t1) const;
    virtual bool need_to_communicate (int& with) const;
    virtual void hint () const;

protected:

#ifdef BL_USE_MPI
    MPI_Request     m_request;
#endif
    FArrayBox*      m_tmp;
    FArrayBox*      m_fab;
    const MultiFab& m_smf;
    Box             m_bx;
    Box             m_region;
    const int       m_sgrid;
    const int       m_target_proc_id;
    bool            m_local;
};

class task_copy_local
    :
    public task_local_base
{
public:

    task_copy_local (task_list&      tl_,
                     FArrayBox*      fab_,
                     int             target_proc_id,
                     const Box&      bx,
                     const MultiFab& mf_,
                     int             grid_);

    virtual bool ready ();
};

//
// A forward declaration.
//
class amr_boundary;

class task_bdy_fill
    :
    public task_local_base
{
public:

    task_bdy_fill (task_list&                tl_,
                   const amr_boundary* bdy_,
                   FArrayBox*                fab_,
                   int                       target_proc_id,
                   const Box&                region_,
                   const MultiFab&           smf_,
                   int                       grid_,
                   const Box&                domain_);

    virtual bool ready ();

protected:

    const Box&          m_domain;
    const amr_boundary* m_bdy;
};

class task_fab
    :
    public task
{
public:

    task_fab (task_list&      tl_,
              const MultiFab& t_,
              int             tt_,
              const Box&      region_,
              int             ncomp_);

    virtual ~task_fab ();

    const FArrayBox& fab ();

protected:

    int target_proc_id () const;
    //
    // The data.
    //
    FArrayBox* target;
    const Box  region;
    const int  ncomp;
    int        m_target_proc_id;
};

inline
const FArrayBox&
task_fab::fab ()
{
    BL_ASSERT(target != 0);
    return *target;
}

inline
int
task_fab::target_proc_id () const
{
    return m_target_proc_id;
}

class level_interface;
class amr_boundary;

class task_fill_patch
    : public task_fab
{
public:

    task_fill_patch (task_list&                tl_,
                     const MultiFab&           t_,
                     int                       tt_,
                     const Box&                region_,
                     const MultiFab&           r_,
                     const level_interface&    lev_interface_,
                     const amr_boundary* bdy_,
                     int                       idim_ /* = -1*/,
                     int                       index_ /*= -1*/);

    virtual ~task_fill_patch ();

private:

    bool fill_patch_blindly ();
    bool fill_exterior_patch_blindly ();
    void fill_patch ();

    const MultiFab&        r;
    const level_interface& lev_interface;
    const amr_boundary*    bdy;
    const int              idim;
    const int              index;
};

class task_fec_base
    :
    public task
{
public:

    task_fec_base (task_list& tl_,
                   MultiFab&  s_,
                   int        igrid_);

    virtual ~task_fec_base ();

protected:

    void push_back (task_fab* tf);
    bool is_local_target () const;
    FArrayBox& target_fab ();
    int grid_number () const;
    const FArrayBox& task_fab_result (int n);
    //
    // The data.
    //
    vector<task::task_proxy> tfvect;
    MultiFab&                s;
    const int                igrid;
};

inline
bool
task_fec_base::is_local_target () const
{
    return is_local(s, igrid);
}

inline
FArrayBox&
task_fec_base::target_fab ()
{
    BL_ASSERT(is_local_target());
    return s[igrid];
}

inline
int
task_fec_base::grid_number () const
{
    return igrid;
}

inline
const FArrayBox&
task_fec_base::task_fab_result (int n)
{
    BL_ASSERT(is_local_target());
    BL_ASSERT(n >= 0 && n < tfvect.size());
    task_fab* tf = dynamic_cast<task_fab*>(tfvect[n].get());
    BL_ASSERT(tf != 0);
    return tf->fab();
}

#endif /*_HGPARALLEL_H_*/
