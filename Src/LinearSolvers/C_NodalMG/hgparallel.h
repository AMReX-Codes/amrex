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

inline int processor_number(const MultiFab&r, int igrid)
{
    return r.DistributionMap()[igrid];
}

inline bool is_remote(const MultiFab& r, int igrid)
{
    return ParallelDescriptor::MyProc() != processor_number(r, igrid);
}

inline bool is_local(const MultiFab& r, int igrid)
{
    return ParallelDescriptor::MyProc() == processor_number(r, igrid);
}

struct HG
{
    static void MPI_init();
    static void MPI_finish();
#ifdef BL_USE_MPI
    static MPI_Comm mpi_comm;
    static int mpi_tag_ub;
#endif
    static int max_live_tasks;
    static int multigrid_maxiter;
    static int cgsolve_maxiter;
    static double cgsolve_tolfact;
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
	explicit task_proxy(task* t_ = 0)
	    : m_t(new task*(t_)), m_cnt( t_ ? new unsigned int(1) : 0), m_finished( t_ ? new bool (false) : 0 )
	{
	}
	task_proxy(const task_proxy& r)
	{
	    link(r);
	}
	~task_proxy()
	{
	    unlink();
	}
	task_proxy& operator=(const task_proxy& r)
	{
	    if ( this != &r )
	    {
		unlink();
		link(r);
	    }
	    return *this;
	}
	task* operator->() const
	{
	    return *m_t;
	}
	task* operator*() const
	{
	    return *m_t;
	}
	task* get() const
	{
	    return *m_t;
	}
	void set_finished()
	{
	    *m_finished = true;
	}
	bool is_finished() const
	{
	    return *m_finished;
	}
    private:
	void unlink()
	{
	    if ( m_cnt )
	    {
		if ( --*m_cnt == 0 )
		{
		    delete *m_t;
		    delete m_t;
		    delete m_cnt;
		    delete m_finished;
		}
		m_t = 0;
		m_cnt = 0;
		m_finished = 0;
	    }
	}
	void link(const task_proxy& r)
	{
	    m_t = r.m_t;
	    m_cnt = r.m_cnt;
	    m_finished = r.m_finished;
	    if ( m_cnt ) ++*m_cnt;
	}
	task** m_t;
	bool* m_finished;
	unsigned int* m_cnt;
    };
    typedef unsigned int sequence_number;
    explicit task(task_list& task_list_);
    virtual ~task();
    virtual bool ready();
    virtual bool startup();
    bool is_started() const;
    virtual bool depends_on_q(const task* t1) const;
    void depend_on(const task_proxy& t1);
    bool depend_ready();
    virtual void hint() const;
    void print_dependencies(ostream& os) const;
    sequence_number get_sequence_number() const;
protected:
    void _do_depend();
    void _hint() const;
    const sequence_number m_sno;
    list< task_proxy > dependencies;
    bool m_started;
    task_list& m_task_list;
private:
    task(const task&);
    void operator==(const task&);
};

// The list...
class task_list
{
public:
    explicit task_list( );
    ~task_list();
    task::task_proxy add_task(task* t);
    void execute();
	// executes once through the task list, return true if any elements left

    list<task::task_proxy>::const_iterator begin() const;
    list<task::task_proxy>::iterator begin();
    list<task::task_proxy>::const_iterator end() const;
    list<task::task_proxy>::iterator end();
    bool empty() const;
    int size() const;
    task::sequence_number get_then_advance();
    void print_dependencies(ostream& os) const;
private:
    list< task::task_proxy > tasks;
    task::sequence_number seq_no;
    int seq_delta;
    bool verbose;
    static bool def_verbose;
};

class task_copy : public task
{
public:
    task_copy(task_list& tl_, MultiFab& mf, int dgrid,                const MultiFab& smf, int sgrid, const Box& bx);
    task_copy(task_list& tl_, MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb);
    virtual ~task_copy();
    virtual bool ready();
    virtual bool depends_on_q(const task* t) const;
    virtual void hint() const;
    virtual bool startup();
protected:
#ifdef BL_USE_MPI
    MPI_Request m_request;
#endif
    FArrayBox* tmp;
    MultiFab& m_mf;
    const MultiFab& m_smf;
    int m_dgrid;
    const int m_sgrid;
    const Box m_bx;
    const Box m_sbx;
    bool m_local;
};

// TASK_COPY_LOCAL

class task_copy_local : public task
{
public:
    task_copy_local(task_list& tl_, FArrayBox* fab_, int target_proc_id, const Box& bx, const MultiFab& mf_, int grid_);
    virtual ~task_copy_local();
    virtual bool ready();
    virtual void hint() const;
    virtual bool startup();
    virtual bool depends_on_q(const task* t1) const;
private:
    bool m_local;
#ifdef BL_USE_MPI
    MPI_Request m_request;
#endif
    FArrayBox* m_fab;
    const Box m_bx;
    FArrayBox* tmp;
    const MultiFab& m_smf;
    const int m_sgrid;
    const int m_target_proc_id;
};

// TASK_FAB

class task_fab : public task
{
public:
    task_fab(task_list& tl_, const MultiFab&t_, int tt_, const Box& region_, int ncomp_);
    virtual ~task_fab();
    virtual const FArrayBox& fab();
protected:
    int target_proc_id() const;
    const Box region;
    const int ncomp;
    FArrayBox* target;
    int m_target_proc_id;
};

class level_interface;
class amr_boundary_class;

class task_fill_patch : public task_fab
{
public:
    task_fill_patch(task_list& tl_, const MultiFab& t_, int tt_, const Box& region_, const MultiFab& r_, const level_interface& lev_interface_, const amr_boundary_class* bdy_, int idim_ /* = -1*/, int index_ /*= -1*/);
    virtual ~task_fill_patch();
private:
    bool fill_patch_blindly();
    bool fill_exterior_patch_blindly();
    void fill_patch();
private:
    const MultiFab& r;
    const level_interface& lev_interface;
    const amr_boundary_class* bdy;
    const int idim;
    const int index;
};

class task_fec_base : public task
{
public:
    task_fec_base(task_list& tl_, MultiFab& s_, int igrid_);
    virtual ~task_fec_base();
protected:
    void push_back(task_fab* tf);
    bool is_local_target() const;
    FArrayBox& target_fab();
    int grid_number() const;
    const FArrayBox& task_fab_result(int n);
private:
    MultiFab& s;
    const int igrid;
    vector< task::task_proxy > tfvect;
};

#endif

