//BL_COPYRIGHT_NOTICE

//
// $Id: Derive.cpp,v 1.12 1999-05-10 18:54:08 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstring>
#else
#include <string.h>
#endif

#include <Derive.H>

DeriveRec::DeriveRec (const aString& name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map,
                      Interpolater*  interp)
    :
    derive_name(name),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(interp),
    bx_map(box_map),
    n_state(0),
    nsr(0),
    rng(0),
    bcr(0)
{}

DeriveRec::~DeriveRec() 
{
   delete [] bcr;
   func     = 0;
   mapper   = 0;
   bx_map   = 0;
   while (rng != 0)
   {
       StateRange* r = rng;
       rng = rng->next;
       delete r;
   }
}

void
DeriveRec::addRange (const DescriptorList& d_list,
                     int                   state_indx,
                     int                   src_comp,
                     int                   num_comp) 
{
    const StateDescriptor& d = d_list[state_indx];

    StateRange* r = new StateRange;

    r->typ  = state_indx;
    r->sc   = src_comp;
    r->nc   = num_comp;
    r->next = 0;
    //
    // Add to end of list.
    //
    if (rng == 0)
    {
        rng = r;
    }
    else
    {
        StateRange* prev = rng;
        while (prev->next != 0)
            prev = prev->next;
        prev->next = r;
    }
    nsr++;
    n_state += num_comp;

    buildBC(d_list);
}

void
DeriveRec::getRange (int  k,
                     int& state_indx,
                     int& src_comp,
                     int& num_comp) const
{
    StateRange* r;

    for (r = rng; r != 0 && k > 0; k--, r = r->next)
        ;
    BL_ASSERT(r != 0);
    state_indx = r->typ;
    src_comp   = r->sc;
    num_comp   = r->nc;
}

void
DeriveRec::buildBC (const DescriptorList& d_list)
{
    BL_ASSERT(nsr > 0);
    delete [] bcr;
    bcr = new int[2*BL_SPACEDIM*n_state];
    int* bci = bcr;
    for (DeriveRec::StateRange* r = rng; r != 0; r = r->next)
    {
        const StateDescriptor& d = d_list[r->typ];

        for (int k = 0; k < r->nc; k++)
        {
            const int* bc = d.getBC(r->sc + k).vect();

            for (int j = 0; j < 2*BL_SPACEDIM; j++)
            {
                bci[j] = bc[j];
            }
            bci += 2*BL_SPACEDIM;
        }
    }
}

bool
DeriveList::canDerive (const aString& name) const 
{
    for (ListIterator<DeriveRec> li(lst); li; ++li)
    {
        if (li().derive_name == name)
            return true;
    }
    return false;
}

const DeriveRec*
DeriveList::get (const aString& name) const
{
    for (ListIterator<DeriveRec> li(lst); li; ++li)
    {
        if (li().derive_name == name)
            return &li();
    }
    return 0;
}

void
DeriveList::addComponent (const aString&        name,
                          const DescriptorList& d_list, 
                          int                   state_indx,
                          int                   s_comp,
                          int                   n_comp)
{
    ListIterator<DeriveRec> li(lst);

    for ( ; li; ++li)
    {
        if (li().derive_name == name)
            break;
    }
    BL_ASSERT (li != 0);
    lst[li].addRange(d_list, state_indx, s_comp, n_comp);
}


