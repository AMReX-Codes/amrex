//BL_COPYRIGHT_NOTICE

//
// $Id: Derive.cpp,v 1.7 1998-03-23 17:23:59 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstring>
#else
#include <string.h>
#endif

#include <Derive.H>

DeriveRec::DeriveRec ()
    :
    bcr(0),
    rng(0)
{}

DeriveRec::DeriveRec (const aString& name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map,
                      IndexType      component_type,
                      Interpolater*  interp)
    :
    derive_name(name),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    func_scr(0),
    mapper(interp),
    rng_type(component_type),
    bx_map(box_map),
    n_state(0),
    nsr(0),
    rng(0),
    bcr(0)
{}

DeriveRec::DeriveRec (const aString& name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFuncScr  der_func_scr,
                      DeriveBoxMap   box_map,
                      IndexType      component_type,
                      Interpolater*  interp)
    :
    derive_name(name),
    der_type(result_type),
    n_derive(nvar_derive),
    func(0),
    func_scr(der_func_scr),
    mapper(interp),
    rng_type(component_type),
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
   func_scr = 0;
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

    assert (d.getType() == rng_type);

    StateRange* r = new StateRange;

    r->typ = state_indx;
    r->sc = src_comp;
    r->nc = num_comp;
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
        StateRange *prev = rng;
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
    assert(r != 0);
    state_indx = r->typ;
    src_comp   = r->sc;
    num_comp   = r->nc;
}

void
DeriveRec::buildBC (const DescriptorList& d_list)
{
    assert(nsr > 0);
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

DeriveList::DeriveList()
    :
    lst()
{}

DeriveList::~DeriveList ()
{
    clear();
    lst.clear();
}

void
DeriveList::clear ()
{
    for (ListIterator<DeriveRec*> li(lst); li; ++li)
    {
        delete lst[li];
        lst[li] = 0;
    }
}

bool
DeriveList::canDerive (const aString& name) const 
{
    for (ListIterator<DeriveRec*> li(lst); li; ++li)
    {
        if (li()->derive_name == name)
            return true;
    }
    return false;
}

const DeriveRec*
DeriveList::get (const aString& name) const
{
    for (ListIterator<DeriveRec*> li(lst); li; ++li)
    {
        if (li()->derive_name == name)
            return li();
    }
    return 0;
}

void
DeriveList::add (const aString& name,
                 IndexType      result_type,
                 int            nvar_derive,
                 DeriveFunc     der_func,
                 DeriveBoxMap   box_map,
                 IndexType      component_type,
                 Interpolater*  interp)
{
    lst.add(new DeriveRec(name,
                          result_type,
                          nvar_derive,
                          der_func,
                          box_map,
                          component_type,
                          interp));
}

void
DeriveList::add (const aString& name,
                 IndexType      result_type,
                 int            nvar_derive,
                 DeriveFuncScr  der_func_scr,
                 DeriveBoxMap   box_map,
                 IndexType      component_type,
                 Interpolater*  interp)
{
    lst.add(new DeriveRec(name,
                          result_type,
                          nvar_derive,
                          der_func_scr,
                          box_map,
                          component_type,
                          interp));
}

void
DeriveList::addComponent (const aString&        name,
                          const DescriptorList& d_list, 
                          int                   state_indx,
                          int                   s_comp,
                          int                   n_comp)
{
    ListIterator<DeriveRec*> li(lst);

    for ( ; li; ++li)
    {
        if (li()->derive_name == name)
            break;
    }
    assert (li != 0);
    lst[li]->addRange(d_list, state_indx, s_comp, n_comp);
}


