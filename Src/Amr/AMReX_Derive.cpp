
#include <cstring>

#include <AMReX_Derive.H>
#include <AMReX_StateDescriptor.H>

namespace amrex {

Box
DeriveRec::TheSameBox (const Box& box) noexcept
{
    return box;
}

Box
DeriveRec::GrowBoxByOne (const Box& box) noexcept
{
    return amrex::grow(box,1);
}


DeriveRec::DeriveRec (const std::string& a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(a_name),
    variable_names(),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(a_interp),
    bx_map(box_map)
{}

DeriveRec::DeriveRec (const std::string& a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc3D   der_func_3d,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(a_name),
    variable_names(),
    der_type(result_type),
    n_derive(nvar_derive),
    func_3d(der_func_3d),
    mapper(a_interp),
    bx_map(box_map)
{}

DeriveRec::DeriveRec (const std::string& a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFuncFab  der_func_fab,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(a_name),
    variable_names(),
    der_type(result_type),
    n_derive(nvar_derive),
    func_fab(der_func_fab),
    mapper(a_interp),
    bx_map(box_map)
{}

// This version doesn't take a Fortran function name, it is entirely defined by the C++
DeriveRec::DeriveRec (const std::string&      a_name,
                      IndexType               result_type,
                      int                     nvar_derive,
                      DeriveRec::DeriveBoxMap box_map)
    :
    derive_name(a_name),
    variable_names(),
    der_type(result_type),
    n_derive(nvar_derive),
    bx_map(box_map)
{}

DeriveRec::DeriveRec (const std::string& a_name,
                      IndexType      result_type,
                      int            nvar_derive,
		      Vector<std::string>& var_names,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(a_name),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(a_interp),
    bx_map(box_map)
{}

DeriveRec::DeriveRec (const std::string& a_name,
                      IndexType      result_type,
                      int            nvar_derive,
		      Vector<std::string>& var_names,
                      DeriveFunc3D   der_func_3d,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(a_name),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func_3d(der_func_3d),
    mapper(a_interp),
    bx_map(box_map)
{}

DeriveRec::DeriveRec (const std::string& a_name,
                      IndexType      result_type,
                      int            nvar_derive,
		      Vector<std::string>& var_names,
                      DeriveFuncFab  der_func_fab,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(a_name),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func_fab(der_func_fab),
    mapper(a_interp),
    bx_map(box_map)
{}

DeriveRec::~DeriveRec () 
{
   delete [] bcr;
   func     = nullptr;
   func_3d  = nullptr;
   func_fab = nullptr;
   mapper   = 0;
   bx_map   = 0;
   while (rng != 0)
   {
       StateRange* r = rng;
       rng = rng->next;
       delete r;
   }
}

const std::string&
DeriveRec::name () const noexcept
{
    return derive_name;
}

IndexType
DeriveRec::deriveType () const noexcept
{
    return der_type;
}

DeriveFunc
DeriveRec::derFunc () const noexcept
{
    return func;
}

DeriveFunc3D
DeriveRec::derFunc3D () const noexcept
{
    return func_3d;
}

DeriveFuncFab
DeriveRec::derFuncFab () const noexcept
{
    return func_fab;
}

DeriveRec::DeriveBoxMap
DeriveRec::boxMap () const noexcept
{
    return bx_map;
}

Interpolater*
DeriveRec::interp () const noexcept
{
    return mapper;
}

int
DeriveRec::numDerive () const noexcept
{
    return n_derive;
}

int
DeriveRec::numRange () const noexcept
{
    return nsr;
}

int
DeriveRec::numState () const noexcept
{
    return n_state;
}

const int*
DeriveRec::getBC () const noexcept
{
    return bcr;
}

void
DeriveRec::addRange (const DescriptorList& d_list,
                     int                   state_indx,
                     int                   src_comp,
                     int                   num_comp) 
{
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
    bcr = new int[2*AMREX_SPACEDIM*n_state];
    int* bci = bcr;
    for (DeriveRec::StateRange* r = rng; r != 0; r = r->next)
    {
        const StateDescriptor& d = d_list[r->typ];

        for (int k = 0; k < r->nc; k++)
        {
            const int* bc = d.getBC(r->sc + k).vect();

            for (int j = 0; j < 2*AMREX_SPACEDIM; j++)
            {
                bci[j] = bc[j];
            }
            bci += 2*AMREX_SPACEDIM;
        }
    }
}

const
std::string&
DeriveRec::variableName(int comp) const noexcept
{
  if (comp < variable_names.size()) 
     return variable_names[comp];

  return derive_name;
}

DeriveList::DeriveList () {}

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveFunc              der_func,
                 DeriveRec::DeriveBoxMap bx_map,
                 Interpolater*           interp)
{
    lst.push_back(DeriveRec(name,result_type,nvar_der,der_func,bx_map,interp));
}

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveFunc3D            der_func_3d,
                 DeriveRec::DeriveBoxMap bx_map,
                 Interpolater*           interp)
{
    lst.push_back(DeriveRec(name,result_type,nvar_der,der_func_3d,bx_map,interp));
}

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveFuncFab           der_func_fab,
                 DeriveRec::DeriveBoxMap bx_map,
                 Interpolater*           interp)
{
    lst.push_back(DeriveRec(name,result_type,nvar_der,der_func_fab,bx_map,interp));
}

// This version doesn't take a Fortran function name, it is entirely defined by the C++
void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveRec::DeriveBoxMap box_map)
{
    lst.push_back(DeriveRec(name,result_type,nvar_der,box_map));
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Vector<std::string>&    vars,
                 DeriveFunc              der_func,
                 DeriveRec::DeriveBoxMap bx_map,
                 Interpolater*           interp)
{
    lst.push_back(DeriveRec(name,res_typ,nvar_der,vars,der_func,bx_map,interp));
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Vector<std::string>&    vars,
                 DeriveFunc3D            der_func_3d,
                 DeriveRec::DeriveBoxMap bx_map,
                 Interpolater*           interp)
{
    lst.push_back(DeriveRec(name,res_typ,nvar_der,vars,der_func_3d,bx_map,interp));
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Vector<std::string>&    vars,
                 DeriveFuncFab           der_func_fab,
                 DeriveRec::DeriveBoxMap bx_map,
                 Interpolater*           interp)
{
    lst.push_back(DeriveRec(name,res_typ,nvar_der,vars,der_func_fab,bx_map,interp));
}

std::list<DeriveRec>&
DeriveList::dlist ()
{
    return lst;
}

bool
DeriveList::canDerive (const std::string& name) const 
{
    for (std::list<DeriveRec>::const_iterator li = lst.begin(), End = lst.end();
         li != End;
         ++li)
    {
        if (li->derive_name == name)
            return true;
    }
    return false;
}

const DeriveRec*
DeriveList::get (const std::string& name) const
{
    for (std::list<DeriveRec>::const_iterator li = lst.begin(), End = lst.end();
         li != End;
         ++li)
    {
        if (li->derive_name == name)
            return &(*li);
    }
    return 0;
}

void
DeriveList::addComponent (const std::string&    name,
                          const DescriptorList& d_list, 
                          int                   state_indx,
                          int                   s_comp,
                          int                   n_comp)
{
    std::list<DeriveRec>::iterator li = lst.begin(), End = lst.end();

    for ( ; li != End; ++li)
    {
        if (li->derive_name == name)
            break;
    }

    BL_ASSERT (li != End);

    li->addRange(d_list, state_indx, s_comp, n_comp);
}

}
