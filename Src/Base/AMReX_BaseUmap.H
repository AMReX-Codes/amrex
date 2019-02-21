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
#include <AMReX_REAL.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_MakeType.H>
#include <unordered_map>


// iv - IntVect in absolute coordiantes
// m - variable / component
// l - multivalued index
#define AMREX_BASEUMAP_KIDX(iv,m,l) (l*domain.numPts()*nvar + m*domain.numPts() + domain.index(iv))

namespace amrex
{

    typedef int key_table_type;
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

        struct Tuple
        {
            Tuple () {}
            Tuple (const IntVect& a_pos, int a_ncomp, int a_l)
                : pos(a_pos), ncomp(a_ncomp), l(a_l) {}
            IntVect pos;
            int ncomp;
            int l;
        };

        struct const_iterator
        {
            const_iterator (const BaseUmap& bmap, int idx) : m_bmap(bmap), m_idx(idx) {}
            void operator ++ () {++m_idx;}
            bool operator< (const const_iterator& rhs)
                {
                    BL_ASSERT(&m_bmap == &(rhs.m_bmap));
                    return m_idx<rhs.m_idx;
                }
            const T& operator* () {
                const Tuple& t = tuple();
                return m_bmap(t.pos, t.ncomp, t.l);
            }

            const key_table_type& index () const {return m_idx;}

            const Tuple& tuple () const {
                typename std::map<key_table_type,Tuple>::const_iterator it=m_bmap.data_key_tuple_map.find(m_idx);
                BL_ASSERT(it != m_bmap.data_key_tuple_map.end());
                return it->second;
            }

        protected:
            const BaseUmap<T>& m_bmap;
            int m_idx;
        };

        const const_iterator begin () const {return const_iterator(*this, data_key_tuple_map.begin()->first);}
        const const_iterator end () const {return const_iterator(*this, data_key_tuple_map.end()->first);}
        size_t size() const {return data_key_tuple_map.size();}

        /**
        * \brief Constructor
        *
        * \param bx         index space covered by this Umap, which
        *      ends up being the size of the array to lookup offsets
        *      into the data array
        * \param a_nvar     number of components
        * \param lmax       maximum number of values in a multi valued cell
        */
        BaseUmap (const Box& bx, int a_nvar, int lmax = 5);
        ~BaseUmap ();

        //! Returns the domain (box) where the array is defined
        const Box& box () const { return domain; }

        //! Returns max_mv
        int MaxMV () const { return max_mv; }

        int numPts () const { return data.size();}

        int nComp () const {return nvar;}

        /**
        * \brief Build the key from the IntVect, component, and multi-valued index
        * For array based implementation this is offset into array
        */
        int key (const IntVect& pos, const int n, const int l) const;

        void set_key (const IntVect& pos,
                      const int n,
                      const int l,
                      const int k );

        //! Store value at (IntVect, ncomp, l)
        void setVal (T val,
                     const IntVect& pos,
                     const int ncomp,
                     const int l);

        T& getVal (const IntVect& pos,
                   const int ncomp,
                   const int l);

        const T& getVal (const IntVect& pos,
                         const int ncomp,
                         const int l) const;

        T& operator() (const IntVect& pos,
                       int            ncomp,
                       int            l);

        const T& operator() (const IntVect& pos,
                             int            ncomp,
                             int            l) const;


        T* dataPtr () { return &data[0]; }

        //! Same as above except works on const UMaps.
        const T* dataPtr (int n = 0) const {return &data[0]; }


        key_table_type* keyTablePtr () { BL_ASSERT(!(key_table == 0)); return key_table; }

        //! Same as above except works on const UMaps.
        const key_table_type* keyTablePtr (int n = 0) const { BL_ASSERT(!(key_table == 0)); return key_table; }


    protected:

        //! Allocates memory for key_table and vector of data
        void define();

        //! Free memory
        void clear ();

        Box domain;
        int nvar;
        int numpts;
        std::vector<T> data;
        key_table_type* key_table;
        std::map<key_table_type,Tuple> data_key_tuple_map;
        int kt_size;
        int max_mv;

    private:
        bool ptr_owner;

    };

    template <class T>
    inline
    T&
    BaseUmap<T>::operator() (const IntVect& pos,
                             int            ncomp,
                             int            l)
    {
        return getVal(pos,ncomp,l);
    }

    template <class T>
    inline
    const T&
    BaseUmap<T>::operator() (const IntVect& pos,
                             int            ncomp,
                             int            l) const
    {
        return getVal(pos,ncomp,l);
    }

    template <class T>
    BaseUmap<T>::BaseUmap (const Box& bx, int a_nvar, int lmax)
        :
        domain(bx),
        nvar(a_nvar),
        numpts(bx.numPts()*lmax*nvar),
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
                data_key_tuple_map.clear();
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
                      const int l ) const
    {
//        int loc = domain.index(pos);
  //      int size = domain.numPts();
        // entry stored in position l*NumComponents + ncomp
//        std::cout << "Looking up key for pos:" << pos << " variable: " << n << " multi-value: " << l << std::endl;
//      int idx = l*
        return key_table[AMREX_BASEUMAP_KIDX(pos,n,l)];

    }

    template <class T>
    void
    BaseUmap<T>::set_key (const IntVect& pos,
                          const int n,
                          const int l,
                          const int k )
    {
        BL_ASSERT(l>=0 && l < max_mv);
        BL_ASSERT(domain.contains(pos));
        BL_ASSERT(n>=0 && n<nvar);
        //int loc = domain.index(pos);
        //int size = domain.numPts();
        // entry stored in position l*NumComponents + ncomp
//        key_table[loc+(l*nvar + n)*size] = k;
        key_table[AMREX_BASEUMAP_KIDX(pos,n,l)] = k;
        data_key_tuple_map[k] = Tuple(pos,n,l);

    }

    template <class T>
    void
    BaseUmap<T>::setVal (T val,
                         const IntVect& pos,
                         const int ncomp,
                         const int l)
    {
        BL_ASSERT(l>=0 && l < max_mv);
        BL_ASSERT(domain.contains(pos));
        BL_ASSERT(ncomp>=0 && ncomp<nvar);
        key_table_type k = key(pos, ncomp, l);

        if ( k >= 0)
        {
            data[k] = val;
        }
        else
        {
            data.push_back(val);
            set_key(pos, ncomp, l, data.size()-1);
        }
        k = key(pos, ncomp, l);

    }

    template <class T>
    const T&
    BaseUmap<T>::getVal (const IntVect& pos,
                         const int ncomp,
                         const int l) const
    {


        BL_ASSERT(l >= 0 && l < max_mv);
        BL_ASSERT(domain.contains(pos));

        key_table_type k = key(pos, ncomp, l);
//        if (k<0) amrex::Abort();
        return data[k];
    }

    template <class T>
    T&
    BaseUmap<T>::getVal (const IntVect& pos,
                         const int ncomp,
                         const int l)
    {

        BL_ASSERT(ncomp >= 0 && ncomp < nvar);
        BL_ASSERT(l >= 0 && l < max_mv);
        BL_ASSERT(domain.contains(pos));

        key_table_type k = key(pos, ncomp, l);
        if (k<0 || k>data.size()){
            //        if (k<0) amrex::Abort();
        }

        return data[k];
    }

} // namespace amrex

#endif
