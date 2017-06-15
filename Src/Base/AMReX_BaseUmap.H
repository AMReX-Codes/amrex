
#ifndef BL_BASEUMAP_H
#define BL_BASEUMAP_H

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxList.H>
#include <AMReX_CArena.H>
#include <AMReX_Looping.H>
#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_MakeType.H>
#include <unordered_map>

namespace amrex
{

/**
 *
 * A sparse storage container for data in an i,j,k,n,m index
 * space. Interface is intended to make something that behaves
 * like a multi-dimensional array
 */

template <class T>
class BaseUmap
{
public:

    typedef T value_type;
    typedef int key_table_type;

    // Constructor
    // bx is index space covered by this Umap, which
    //      ends up being the size of the array to lookup offsets
    //      into the data array
    // a_nvar is number of components
    // lmax is maximum number of values in a multi valued cell
    BaseUmap (const Box& bx, int a_nvar, int lmax = 5);
    ~BaseUmap ();
    
  



    // Build the key from the IntVect, component, and multi-valued index
    // For array based implementation this is offset into array
    int key (const IntVect& pos, const int n, const int l);

    void set_key (const IntVect& pos, 
              const int n,
              const int l,
              const int k );

    // Store value at (IntVect, ncomp, l)
    void setVal (T val,
                 const IntVect& pos,
                 const int ncomp,
                 const int l);

    T getVal (const IntVect& pos,
              const int ncomp,
              const int l);

protected:
    
    //! Allocates memory for key_table and vector of data
    void define();

    //! Free memory
    void clear ();

    Box domain;
    int nvar;
    std::vector<T> data;
    key_table_type* key_table;
    long kt_size;
    long numpts;
    int max_mv;

    // std::unordered_map<long, T> data;

private:
    bool ptr_owner;

};

template <class T>
BaseUmap<T>::BaseUmap (const Box& bx, int a_nvar, int lmax)
:
domain(bx),
numpts(bx.numPts()*lmax*nvar),
nvar(a_nvar),
key_table(0),
max_mv(lmax),
ptr_owner(false)
{
    define();
}

template <class T>
void
BaseUmap<T>::clear ()
{
    if (key_table)
    {
        //
        // Call T::~T() on the to-be-destroyed memory.
        //
    if (ptr_owner)
    {
        amrex::The_Arena()->free(key_table);
    }

    key_table = 0;
    kt_size = 0;
    }
}


template <class T>
void
BaseUmap<T>::define ()
{
    BL_ASSERT(nvar > 0);
    BL_ASSERT(dptr == 0);
    BL_ASSERT(numpts > 0);

    kt_size  = nvar*numpts*max_mv;
    key_table = static_cast<key_table_type*>(amrex::The_Arena()->alloc(kt_size*sizeof(key_table_type)));
    ptr_owner = true;

    for(int i=0; i<kt_size; ++i)
    {
        key_table[i] = -1;
    }

}

template <class T>
BaseUmap<T>::~BaseUmap ()
{
    clear();
}   

template <class T>
int
BaseUmap<T>::key (const IntVect& pos, 
              const int n,
              const int l )
{
    long loc = domain.index(pos);
    long size = domain.numPts();
    // entry stored in position l*NumComponents + ncomp
    return key_table[loc+(l*nvar + n)*size];

}

template <class T>
void
BaseUmap<T>::set_key (const IntVect& pos, 
              const int n,
              const int l,
              const int k )
{
    long loc = domain.index(pos);
    long size = domain.numPts();
    // entry stored in position l*NumComponents + ncomp
    key_table[loc+(l*nvar + n)*size] = k;

}

template <class T>
void
BaseUmap<T>::setVal (T val,
                 const IntVect& pos,
                 const int ncomp,
                 const int l)
{

    long k = key(pos, ncomp, l);

    if ( k >= 0)
    {
        data[k] = val;
    }
    else
    {
        data.push_back(val);
        set_key(pos, ncomp, l, data.size()-1);
        long k = key(pos, ncomp, l);
        T valout = getVal(pos, ncomp, l);
        std::cout << "Storing value " << val << " in offset " << k << " got back value " << valout << std::endl;
    }
}

template <class T>
T
BaseUmap<T>::getVal (const IntVect& pos,
                     const int ncomp,
                     const int l)
{

    key_table_type k = key(pos, ncomp, l);
    if (k >= 0)
    {
        return data[k];
    }
    else
    {
        //amrex::Abort();
        return -1.0;

    }

}

} // namespace amrex

#endif