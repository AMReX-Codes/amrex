//BL_COPYRIGHT_NOTICE

//
// $Id: tCArena.cpp,v 1.2 1998-02-08 21:13:31 lijewski Exp $
//

#ifndef WIN32
#include <unistd.h>
#endif

#include <REAL.H>
#include <CArena.H>
#include <Utility.H>

#ifdef BL_USE_NEW_HFILES
#include <list>
#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <list.h>
#include <new.h>
#endif

//
// A simple class emulating how we use FABs.
//
class FB
{
public:
    FB ();
    ~FB ();

    bool ok () const;

    enum { CHUNKSIZE = 1024 };

private:
    //
    // Disallowed
    //
    FB (const FB& rhs);
    FB& operator= (const FB&);

    static CArena m_CArena;

    size_t  m_size;
    double* m_data;
};

CArena FB::m_CArena(100*CHUNKSIZE);

FB::FB ()
{
    m_size = CHUNKSIZE*Utility::Random();
    m_CArena.alloc(m_size*sizeof(double), (void**) &m_data);
    //
    // Set specific values in the data.
    //
    for (int i = 0; i < m_size; i++)
        m_data[i] = m_size;
}

FB::~FB ()
{
    ok();
    m_CArena.free(m_data);
}

bool
FB::ok () const
{
    for (int i = 0; i < m_size; i++)
        assert(m_data[i] == m_size);
    return true;
}

int
main ()
{
    set_new_handler(Utility::OutOfMemory);

    list<FB*> fbl;

    for (int j = 0; j < 10; j++)
    {
        cout << "Loop == " << j << endl;

        for (int i = 0; i < 1000; i++)
        {
            fbl.push_back(new FB);
        }

        while (!fbl.empty())
        {
            delete fbl.back();
            fbl.pop_back();
        }
    }

    return 0;
}
