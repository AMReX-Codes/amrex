//BL_COPYRIGHT_NOTICE

#ifndef BL_CARENA
#define BL_CARENA

#ifdef BL_USE_NEW_HFILES
#include <cstddef>
#include <iostream>
using std::ostream;
#else
#include <stddef.h>
#include <iostream.h>
#endif

#include <Boolean.H>
#include <Arena.H>
#include <Array.H>

class CArena
    :
    public Arena
{
public:

    CArena (size_t size = 0);

    ~CArena ();

    void stat (ostream& os, int level = 0x01);

    virtual void* alloc (size_t nbytes, void** ind);
    virtual void free (void* ap);

    void compact ();

    static void verbose (bool _verbose) 
    {
	m_verbose = _verbose;
    }
    static bool verbose () 
    {
	return m_verbose;
    }

protected:

    union Header
    {
	struct
	{
	    Header* frnext;
	    Header* frprev;
	    Header* blnext;
	    Header* blprev;
	    void**  callback;
	    size_t  size;
	} s;
	Word x;
    };

    bool ok (int level);
    Header* freeFree (Header* lower, Header* upper);
    Header* freeUsed (Header* freeblock, Header* busyblock);
    Header* morecore (size_t);
    void dump (const char*, void*);
    //
    // Static data.
    //
    static bool m_verbose;
    //
    // Local data.
    //
    Array<void*> psaved;
    Header       base;
    Header*      highest_mem;
    Header*      lowest_mem;
    Header*      freep;
    size_t       c_hunk_size;

private:
    //
    // Disallowed.
    //
    CArena (const CArena&);
    CArena& operator= (const CArena&);
};

//
// The CArena used by BaseFab code.
//
extern CArena* The_FAB_CArena;

#endif /*BL_CARENA*/
