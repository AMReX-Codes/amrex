//BL_COPYRIGHT_NOTICE

//
// $Id: Derive.cpp,v 1.3 1997-11-26 20:41:43 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cstring>
#else
#include <string.h>
#endif

#include <Derive.H>

// ------------------------------------------------------------------
DeriveRec::DeriveRec()
    : derive_name(), bcr(0), rng(0)
{}

// ------------------------------------------------------------------
DeriveRec::DeriveRec(const aString &name, IndexType result_type,
                     int nvar_derive, DeriveFunc der_func,
		     DeriveBoxMap box_map, IndexType component_type,
		     Interpolater *interp) :
		     der_type(result_type), n_derive(nvar_derive),
		     func(der_func), bx_map(box_map),
		     rng_type(component_type), mapper(interp)
{
    derive_name = name;
    n_state = 0;
    nsr = 0;
    bcr = 0;
    rng = 0;
}

// ------------------------------------------------------------------
DeriveRec::~DeriveRec() 
{
   delete bcr;
   func = 0;
   mapper = 0;
   bx_map = 0;
   while (rng != 0) {
      StateRange *r = rng;
      rng = rng->next;
      delete r;
   }
}

// ------------------------------------------------------------------
void
DeriveRec::addRange(const DescriptorList &d_list, int state_indx,
		    int src_comp, int num_comp) 
{
    const StateDescriptor &d = d_list[state_indx];
    assert (d.getType() == rng_type);
    StateRange *r = new StateRange;
    if (r == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    r->typ = state_indx;
    r->sc = src_comp;
    r->nc = num_comp;
    r->next = 0;
      // add to end of list
    if (rng == 0) {
	rng = r;
    } else {
	StateRange *prev = rng;
	while (prev->next != 0) prev = prev->next;
	prev->next = r;
    }
    nsr++;
    n_state += num_comp;

    buildBC(d_list);
}

// ------------------------------------------------------------------
void
DeriveRec::getRange(int k, int& state_indx, int& src_comp,
		    int& num_comp) const
{
    StateRange *r;
    for (r = rng; r!=0 && k > 0; k--, r=r->next);
    assert(r != 0);
    state_indx = r->typ;
    src_comp = r->sc;
    num_comp = r->nc;
}

// ------------------------------------------------------------------
void
DeriveRec::buildBC(const DescriptorList& d_list)
{
    assert(nsr > 0);
    delete bcr;
    if ((bcr = new int[2*BL_SPACEDIM*n_state]) == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    int *bci = bcr;
    DeriveRec::StateRange *r;
    for (r = rng; r != 0; r = r->next) {
	const StateDescriptor &d = d_list[r->typ];
        int k;
	for (k = 0; k < r->nc; k++) {
	    int n = r->sc + k;
	    const int* bc = d.getBC(n).vect();
            int j;
	    for (j = 0; j < 2*BL_SPACEDIM; j++) {
		bci[j] = bc[j];
	    }
	    bci += 2*BL_SPACEDIM;
	}
    }
}

// ------------------------------------------------------------------
DeriveList::DeriveList()
    : lst()
{
}

// ------------------------------------------------------------------
DeriveList::~DeriveList()
{
    clear();
    lst.clear();
}

// ------------------------------------------------------------------
void
DeriveList::clear()
{
    ListIterator<DeriveRec*> li(lst);
    while (li) {
	delete lst[li];
	lst[li] = 0;
	++li;
    }
}

// ------------------------------------------------------------------
int
DeriveList::canDerive(const aString &name) const 
{
    ListIterator<DeriveRec*> li(lst);
    while(li) {
	if(li()->derive_name == name)
	    return true;
	++li;
    }
    return false;
}

// ------------------------------------------------------------------
const DeriveRec*
DeriveList::get(const aString &name) const
{
    ListIterator<DeriveRec*> li(lst);
    while(li) {
	if(li()->derive_name == name)
	    return li();
	++li;
    }
    return 0;
}

// ------------------------------------------------------------------
void
DeriveList::add(const aString &name, IndexType result_type,
		int nvar_derive, DeriveFunc der_func,
		DeriveBoxMap box_map, IndexType component_type,
		Interpolater *interp)
{
    DeriveRec* dr = new DeriveRec(name,result_type,nvar_derive,
                                  der_func,box_map,component_type,interp);
    if (dr == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    lst.add(dr);
}

// ------------------------------------------------------------------
void
DeriveList::addComponent(const aString &name,
			 const DescriptorList& d_list, 
			 int state_indx, int s_comp, int n_comp)
{
    ListIterator<DeriveRec*> li(lst);
    while(li) {
	if(li()->derive_name == name) {
	    break;
	}
	++li;
    }
    assert (li != 0);
    lst[li]->addRange(d_list,state_indx,s_comp,n_comp);
}


