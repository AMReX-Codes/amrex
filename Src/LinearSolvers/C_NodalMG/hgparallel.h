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
#error Must define BL_USE_MPI in this file
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

#ifdef BL_USE_NEW_HFILES
#include <list>
using namespace std;
#else
#include <list.h>
#endif

class task
{
public:
    typedef unsigned int sequence_number;
    task() : m_sno(~0U), m_started(false) {}
    virtual ~task() {}
    virtual bool ready() = 0;
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	m_sno  = sno;
	m_comm = comm;
	return false;
    }
    virtual void hint() const { };
    void set_sequence_number(sequence_number sno)
    {
	m_sno = sno;
    }
    sequence_number get_sequence_number() const
    {
	return m_sno;
    }
    virtual bool depends_on_q(const task* t1) const
    {
	return false;
    }
    void depend_on(task*& t1)
    {
	dependencies.push_back( &t1 );
    }
    bool depend_ready();
protected:
    sequence_number m_sno;
    list< task** > dependencies;
    bool m_started;
    MPI_Comm m_comm;
};

class task_copy : public task
{
public:
    task_copy(MultiFab& mf, int dgrid,                const MultiFab& smf, int sgrid, const Box& bx);
    task_copy(MultiFab& mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb);
    virtual ~task_copy();
    virtual bool ready();
    virtual bool init(sequence_number sno, MPI_Comm comm);
    virtual void hint() const;
    virtual bool depends_on_q(const task* t) const;
protected:
    void startup();
    MPI_Request m_request;
    FArrayBox* tmp;
    MultiFab& m_mf;
    const MultiFab& m_smf;
    int m_dgrid;
    const int m_sgrid;
    const Box m_bx;
    const Box m_sbx;
    bool m_local;
};

class task_copy_local : public task
{
public:
    task_copy_local(FArrayBox& fab_, const MultiFab& mf_, int grid_, const Box& bx);
    virtual ~task_copy_local();
    virtual bool ready();
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	abort(); return false;
    }
private:
    FArrayBox& m_fab;
    const MultiFab& m_smf;
    const int m_sgrid;
    const Box m_bx;
    const Box m_sbx;
    bool m_local;
};

class task_fab : public task
{
public:
    virtual const FArrayBox& fab() = 0;
    virtual bool init(sequence_number sno, MPI_Comm comm) = 0;
};

class task_fab_get : public task_fab
{
public:
    task_fab_get(const MultiFab& d_, int dgrid_, const MultiFab& s_, int sgrid_, const Box& bx);
    virtual ~task_fab_get();
    virtual const FArrayBox& fab();
    virtual bool ready();
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	abort(); return false;
    }
private:
    const MultiFab& s;
    const MultiFab& d;
    const int dgrid;
    const int sgrid;
    const Box bx;
};

class task_list
{
public:
    task_list(MPI_Comm comm = MPI_COMM_WORLD);
    ~task_list();
    void add_task(task* t);
    void execute();
private:
    list< task** > tasks;
    MPI_Comm comm;
    task::sequence_number seq_no;
    bool verbose;
    static bool def_verbose;
};

class level_interface;
class amr_boundary_class;

class task_fill_patch : public task_fab
{
public:
    task_fill_patch(const Box& region_, const MultiFab& r_, const level_interface& lev_interface_, const amr_boundary_class* bdy_, int idim_ = 0, int index_ = 0);
    virtual ~task_fill_patch();
    virtual const FArrayBox& fab();
    virtual bool ready();
    virtual bool init(sequence_number sno, MPI_Comm comm);
private:
    bool fill_patch_blindly();
    bool fill_exterior_patch_blindly();
    void fill_patch();
    bool newed;
    FArrayBox* target;
    const Box region;
    const MultiFab& r;
    const level_interface& lev_interface;
    const amr_boundary_class* bdy;
    const int idim;
    const int index;
    task_list tl;
};


class task_linked_task : public task
{
public:
    task_linked_task(task* t_) : lcpt(t_) {}
    virtual bool ready() { return lcpt->ready(); }
    virtual bool init(sequence_number no, MPI_Comm comm) { return lcpt->init(no, comm); }
private:
    LnClassPtr<task> lcpt;
};

class task_copy_link : public task
{
public:
    task_copy_link(MultiFab& m_, int jgrid_, int igrid_, const Box& freg_, const task_linked_task& t_)
	: m(m_), jgrid(jgrid_), igrid(igrid_), freg(freg_), t(t_) {}
    virtual bool ready()
    {
	abort();
	if ( t.ready() )
	{
	    m[jgrid].copy(m[igrid], freg);
	    return true;
	}
	return false;
    }
    virtual bool init(sequence_number sno, MPI_Comm comm)
    {
	abort(); return false;
    }
private:
    MultiFab& m;
    const int igrid;
    const int jgrid;
    const Box freg;
    task_linked_task t;
};



#endif

