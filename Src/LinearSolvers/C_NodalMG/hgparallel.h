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
#else
//#error Must define BL_USE_MPI in this file
#endif

inline int processor_number (const MultiFab& r,
                             int             igrid)
{
    return r.DistributionMap()[igrid];
}

inline bool is_remote (const MultiFab& r,
                       int             igrid)
{
    return ParallelDescriptor::MyProc() != processor_number(r, igrid);
}

inline bool is_local (const MultiFab& r,
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
    static bool initialized;
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
        explicit task_proxy (task* t_ = 0)
        {
            if (t_)
            {
                m_t        = new task*(t_);
                m_cnt      = new unsigned int(1);
                m_finished = new bool(false);
            }
            else
            {
                m_t        = 0;
                m_cnt      = 0;
                m_finished = 0;
            }
        }
        task_proxy (const task_proxy& r)
        {
            link(r);
        }
        ~task_proxy ()
        {
            unlink();
        }
        task_proxy& operator= (const task_proxy& r)
        {
            if (this != &r)
            {
                unlink();
                link(r);
            }
            return *this;
        }
        task* operator-> () const
        {
            assert(!null());
            return *m_t;
        }
        task* operator* () const
        {
            assert(!null());
            return *m_t;
        }
        task* get () const
        {
            assert(!null());
            return *m_t;
        }
        void set_finished ()
        {
            *m_finished = true;
        }
        bool is_finished () const
        {
            return *m_finished;
        }
        bool null () const
        {
            return m_t == 0;
        }
    private:
        void unlink ()
        {
            if (m_cnt)
            {
                if (--*m_cnt == 0)
                {
                    delete *m_t;
                    delete m_t;
                    delete m_cnt;
                    delete m_finished;
                }
                m_t        = 0;
                m_cnt      = 0;
                m_finished = 0;
            }
        }
        void link (const task_proxy& r)
        {
            m_t        = r.m_t;
            m_cnt      = r.m_cnt;
            m_finished = r.m_finished;
            if (m_cnt)
                ++*m_cnt;
        }
        //
        // The data.
        //
        task**        m_t;
        bool*         m_finished;
        unsigned int* m_cnt;
    };

    typedef unsigned int sequence_number;
    explicit task (task_list& task_list_);
    virtual ~task ();
    virtual bool ready ();
    virtual bool startup ();
    bool is_started () const { return m_started; }
    virtual bool depends_on_q (const task* t1) const;
    virtual bool work_to_do () const;
    void depend_on (const task_proxy& t1);
    bool depend_ready ();
    virtual void hint () const;
    void print_dependencies (ostream& os) const;
    sequence_number get_sequence_number () const { return m_sno; }
protected:
    void _do_depend ();
    void _hint () const;
    //
    // The data.
    //
    const sequence_number m_sno;
    list<task_proxy>      dependencies;
    bool                  m_started;
    task_list&            m_task_list;
private:
    //
    // Not defined.
    //
    task (const task&);
    task& operator= (const task&);
};

class task_list
{
public:
    explicit task_list ();
    ~task_list ();
    task::task_proxy add_task (task* t);
    //
    // Executes once through the task list, return true if any elements left.
    //
    void execute ();

    list<task::task_proxy>::const_iterator begin () const { return tasks.begin(); }
    list<task::task_proxy>::iterator begin () { return tasks.begin(); }
    list<task::task_proxy>::const_iterator end () const { return tasks.end(); }
    list<task::task_proxy>::iterator end () { return tasks.end(); }

    bool empty () const { return tasks.empty(); }
    int size () const { return tasks.size(); }

    task::sequence_number get_then_advance ()
    {
        task::sequence_number tmp = seq_no;
        seq_no += seq_delta;
        return tmp;
    }
    void print_dependencies (ostream& os) const;
private:
    //
    // The data.
    //
    list< task::task_proxy > tasks;
    task::sequence_number    seq_no;
    int                      seq_delta;
    bool                     verbose;
    static bool              def_verbose;
};

class task_copy : public task
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

    virtual ~task_copy ();
    virtual bool ready ();
    virtual bool depends_on_q (const task* t) const;
    virtual void hint () const;
    virtual bool startup ();
    virtual bool work_to_do () const;
protected:
    //
    // Common function called by constructors.
    //
    void init ();
    //
    // The data.
    //
    FArrayBox*      tmp;
    MultiFab&       m_mf;
    const MultiFab& m_smf;
    int             m_dgrid;
    const int       m_sgrid;
    const Box       m_bx;
    const Box       m_sbx;
#ifdef BL_USE_MPI
    MPI_Request     m_request;
#endif
    bool            m_local;
    bool            m_done;
};

class task_copy_local : public task
{
public:
    task_copy_local (task_list&      tl_,
                     FArrayBox*      fab_,
                     int             target_proc_id,
                     const Box&      bx,
                     const MultiFab& mf_,
                     int             grid_);

    virtual ~task_copy_local ();
    virtual bool ready ();
    virtual void hint () const;
    virtual bool startup ();
    virtual bool depends_on_q (const task* t1) const;
    virtual bool work_to_do () const;
private:
    //
    // The data.
    //
    FArrayBox*      tmp;
    FArrayBox*      m_fab;
    const MultiFab& m_smf;
    const Box       m_bx;
    const int       m_sgrid;
    const int       m_target_proc_id;
#ifdef BL_USE_MPI
    MPI_Request     m_request;
#endif
    bool            m_local;
    bool            m_done;
};

class task_fab : public task
{
public:
    task_fab (task_list&      tl_,
              const MultiFab& t_,
              int             tt_,
              const Box&      region_,
              int             ncomp_);

    virtual ~task_fab ();
    virtual const FArrayBox& fab ();
    virtual bool work_to_do () const;
protected:
    int target_proc_id () const { return m_target_proc_id; }
    //
    // The data.
    //
    const Box  region;
    const int  ncomp;
    FArrayBox* target;
    int        m_target_proc_id;
};

class level_interface;
class amr_boundary_class;

class task_fill_patch : public task_fab
{
public:
    task_fill_patch (task_list&                tl_,
                     const MultiFab&           t_,
                     int                       tt_,
                     const Box&                region_,
                     const MultiFab&           r_,
                     const level_interface&    lev_interface_,
                     const amr_boundary_class* bdy_,
                     int                       idim_ /* = -1*/,
                     int                       index_ /*= -1*/);

    virtual ~task_fill_patch ();
    virtual bool work_to_do () const;
private:
    bool fill_patch_blindly ();
    bool fill_exterior_patch_blindly ();
    void fill_patch ();
private:
    const MultiFab&           r;
    const level_interface&    lev_interface;
    const amr_boundary_class* bdy;
    const int                 idim;
    const int                 index;
};

class task_fec_base : public task
{
public:
    task_fec_base (task_list& tl_,
                   MultiFab&  s_,
                   int        igrid_);

    virtual ~task_fec_base ();
protected:
    void push_back (task_fab* tf);
    bool is_local_target () const { return is_local(s, igrid); }
    FArrayBox& target_fab ()
    {
        assert(is_local_target());
        return s[igrid];
    }
    int grid_number () const { return igrid; }
    const FArrayBox& task_fab_result (int n)
    {
        assert(is_local_target());
        assert(n >= 0 && n < tfvect.size());
        task_fab* tf = dynamic_cast<task_fab*>(tfvect[n].get());
        assert(tf != 0);
        return tf->fab();
    }
    virtual bool work_to_do () const;
private:
    //
    // The data.
    //
    MultiFab&                s;
    const int                igrid;
    vector<task::task_proxy> tfvect;
};

#endif /*_HGPARALLEL_H_*/

