
#include <winstd.H>

#include <BaseFab.H>
#include <BArena.H>
#include <CArena.H>
#include <Thread.H>

int BoxLib::BF_init::m_cnt = 0;

static ThreadSpecificData<Arena>* arena = 0;

BoxLib::BF_init::BF_init ()
{
    if (m_cnt++ == 0)
        arena = new ThreadSpecificData<Arena>;
}

BoxLib::BF_init::~BF_init ()
{
    if (--m_cnt == 0)
        delete arena;
}

Arena*
BoxLib::ResetArena (Arena* newarena)
{
    BL_ASSERT(newarena != 0);

    Arena* oldarena = arena->get();

    BL_ASSERT(oldarena != 0);

    arena->set(newarena);

    return oldarena;
}

Arena*
BoxLib::The_Arena ()
{
    BL_ASSERT(arena != 0);

    Arena* a = arena->get();

    if (a == 0)
    {
#if defined(BL_THREADS) || defined(BL_COALESCE_FABS)
        arena->set(a = new CArena);
#else
        arena->set(a = new BArena);
#endif
    }

    return a;
}

