//BL_COPYRIGHT_NOTICE

#ifndef _SIMPLEDPTR_H_
#define _SIMPLEDPTR_H_

//
// $Id: SimpleDPtr.H,v 1.1 1997-09-12 18:00:16 lijewski Exp $
//

#include <Assert.H>
#include <stddef.h>
#include <DPtr.H>

//
//@Man:
//@Memo: A Simple Smart Data Pointer
/*@Doc:

  This is the simplest class derived from DPtrRep, supplying DPtr function
  using a fixed memory arena provided by BArena.

  This class does NOT provide a copy constructor or assignment operator.
*/

template <class T>
class SimpleDPtr
    : public DPtrRep<T>
{
public:
 
    /*@ManDoc: Construct a new object using the given arena, if supplied,
               or a BArena otherwise.  size is the number of objects of
               type T that the arena will hold.  Calls define(size).
    */
    SimpleDPtr (size_t size,
                Arena* arena = 0);
    //
    //@ManDoc: The virtual destructor -- calls clear().
    //
    virtual ~SimpleDPtr ();
    //
    //@ManDoc: Free the arena space.
    //
    virtual void clear ();
    //
    //@ManDoc: Change the size of the arena to newsize.
    //
    virtual void resize (size_t newsize);

    /*@ManDoc: Set the initial arena size.  The size is specified as the
               number of objects of type T that the arena will hold.
               It is illegal to call this once the size of the arena
               has already been set; i.e. don't call this more than once.
    */
    virtual void define (size_t size);
 
    /*@ManDoc: Smart array access.  This will abort if the pointer is null
               or if the argument is out of range.
    */
    virtual T& operator[] (long) const;
    //
    //@ManDoc: Return the arena size in units the size of type T.
    //
    virtual size_t length ();
    //
    //@ManDoc: Return the arena size in units the size of type T.
    //
    virtual size_t size () const;

protected:
    //
    // The default constructor.  Here so we can be a parent class.
    //
    inline SimpleDPtr ();
    //
    // The arena pointer.
    //
    T* dp;
    //
    // The size of the arena.
    //
    size_t currentsize;

private:
    //
    // Disallow these.
    //
    SimpleDPtr (const SimpleDPtr<T> &);
    SimpleDPtr<T>& operator= (const SimpleDPtr<T>&);
};

//
// Inlines.
//

template <class T>
inline
SimpleDPtr<T>::SimpleDPtr ()
        : dp(0),
          currentsize(0)
{}


#include <BArena.H>

template <class T>
SimpleDPtr<T>::SimpleDPtr (size_t _size,
                           Arena* _arena)
    : DPtrRep<T>(_arena),
      dp(0)
{
    if (arena == 0)
    {
        static BArena _builtin_BArena;
        arena = &_builtin_BArena;
    }
    define(_size);
}

template <class T>
SimpleDPtr<T>::~SimpleDPtr ()
{
    clear();
}

template <class T>
void
SimpleDPtr<T>::clear ()
{
    if (arena)
        arena->free(dp);
    dp = 0;
    currentsize = 0;
}

template <class T>
void
SimpleDPtr<T>::resize (size_t _size)
{
    if (_size > currentsize)
    {
        T *g = (T*) arena->alloc(_size * sizeof(T));
        for (size_t i = 0; i < currentsize; ++i)
            g[i] = dp[i];
        delete [] dp;
        dp = g;
        currentsize = _size;
    }
}

template <class T>
void
SimpleDPtr<T>::define (size_t _size)
{
    assert(dp == 0);
    dp = (T*) arena->alloc(_size*sizeof(T));
    currentsize = _size;
    assert(dp != 0);
}

template <class T>
T&
SimpleDPtr<T>::operator[] (long n) const
{
    assert(n >= 0 && n < long(currentsize));
    return dp[n];
}

template <class T>
size_t
SimpleDPtr<T>::length ()
{
    return currentsize;
}

template <class T>
size_t
SimpleDPtr<T>::size () const
{
    return currentsize;
}

#endif /*_SIMPLEDPTR_H_*/
