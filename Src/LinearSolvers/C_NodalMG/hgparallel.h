#ifndef _HGPARALLEL_H_
#define _HGPARALLEL_H_

#include <MultiFab.H>

inline bool is_remote(const MultiFab& r, int igrid)
{
    if (ParallelDescriptor::MyProc() == r.DistributionMap()[igrid]) return false;
    return true;
}

inline bool is_local(const MultiFab& r, int igrid)
{
    return ! is_remote(r, igrid);
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
    virtual ~task() {}
    virtual bool ready() = 0;
};

class task_copy : public task
{
public:
    task_copy(FArrayBox* fab_, const MultiFab& mf_, int grid_, const Box& bx);
    task_copy(MultiFab* mf, int dgrid, const MultiFab& smf_, int sgrid, const Box& bx);
    task_copy(MultiFab* mf, int dgrid, const Box& db, const MultiFab& smf, int sgrid, const Box& sb);
    virtual ~task_copy();
    virtual bool ready();
private:
    FArrayBox* m_fab;
    MultiFab* m_mf;
    const MultiFab& m_smf;
    int m_dgrid;
    const int m_sgrid;
    const Box m_bx;
    const Box s_bx;
    bool m_ready;
};

class task_list
{
public:
    void add_task(task* t);
    void execute();
private:
    list<task*> tasks;
};


#endif

