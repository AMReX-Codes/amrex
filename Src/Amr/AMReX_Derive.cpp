
#include <AMReX_Derive.H>
#include <AMReX_StateDescriptor.H>

#include <cstring>

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


DeriveRec::DeriveRec (std::string    a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(std::move(a_name)),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(a_interp),
    bx_map(std::move(box_map))
{}

DeriveRec::DeriveRec (std::string    a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFunc3D   der_func_3d,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(std::move(a_name)),
    der_type(result_type),
    n_derive(nvar_derive),
    func_3d(der_func_3d),
    mapper(a_interp),
    bx_map(std::move(box_map))
{}

DeriveRec::DeriveRec (std::string    a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      DeriveFuncFab  der_func_fab,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(std::move(a_name)),
    der_type(result_type),
    n_derive(nvar_derive),
    func_fab(std::move(der_func_fab)),
    mapper(a_interp),
    bx_map(std::move(box_map))
{}

// This version doesn't take a Fortran function name, it is entirely defined by the C++
DeriveRec::DeriveRec (std::string             a_name,
                      IndexType               result_type,
                      int                     nvar_derive,
                      DeriveRec::DeriveBoxMap box_map)
    :
    derive_name(std::move(a_name)),
    der_type(result_type),
    n_derive(nvar_derive),
    bx_map(std::move(box_map))
{}

DeriveRec::DeriveRec (std::string    a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      Vector<std::string> const& var_names,
                      DeriveFunc     der_func,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(std::move(a_name)),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func(der_func),
    mapper(a_interp),
    bx_map(std::move(box_map))
{}

DeriveRec::DeriveRec (std::string    a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      Vector<std::string> const& var_names,
                      DeriveFunc3D   der_func_3d,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(std::move(a_name)),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func_3d(der_func_3d),
    mapper(a_interp),
    bx_map(std::move(box_map))
{}

DeriveRec::DeriveRec (std::string    a_name,
                      IndexType      result_type,
                      int            nvar_derive,
                      Vector<std::string> const& var_names,
                      DeriveFuncFab  der_func_fab,
                      DeriveBoxMap   box_map,
                      Interpolater*  a_interp)
    :
    derive_name(std::move(a_name)),
    variable_names(var_names),
    der_type(result_type),
    n_derive(nvar_derive),
    func_fab(std::move(der_func_fab)),
    mapper(a_interp),
    bx_map(std::move(box_map))
{}

DeriveRec::~DeriveRec ()
{
   delete [] bcr;
   delete [] bcr3D;
   func     = nullptr;
   func_3d  = nullptr;
   func_fab = nullptr;
   mapper   = nullptr;
   bx_map   = nullptr;
   while (rng != nullptr)
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
      BL_ASSERT( bcr != nullptr );
      return bcr;
}

const int*
DeriveRec::getBC3D () const noexcept
{
      BL_ASSERT( bcr3D != nullptr );
      return bcr3D;
}

void
DeriveRec::addRange (const DescriptorList& d_list,
                     int                   state_indx,
                     int                   src_comp,
                     int                   num_comp)
{
    auto* r = new StateRange;

    r->typ  = state_indx;
    r->sc   = src_comp;
    r->nc   = num_comp;
    r->next = nullptr;
    //
    // Add to end of list.
    //
    if (rng == nullptr)
    {
        rng = r;
    }
    else
    {
        StateRange* prev = rng;
        while (prev->next != nullptr) {
            prev = prev->next;
        }
        prev->next = r;
    }
    nsr++;
    n_state += num_comp;

    buildBC(d_list);
    buildBC3D(d_list);
}

void
DeriveRec::getRange (int  k,
                     int& state_indx,
                     int& src_comp,
                     int& num_comp) const
{
    StateRange* r;

    for (r = rng; r != nullptr && k > 0; k--, r = r->next) {
        ;
    }
    BL_ASSERT(r != nullptr);
    if (r != nullptr) {
        state_indx = r->typ;
        src_comp   = r->sc;
        num_comp   = r->nc;
    }
}

void
DeriveRec::buildBC (const DescriptorList& d_list)
{
    BL_ASSERT(nsr > 0);
    delete [] bcr;
    bcr = new int[2*AMREX_SPACEDIM*n_state];
    int* bci = bcr;
    for (DeriveRec::StateRange* r = rng; r != nullptr; r = r->next)
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

void
DeriveRec::buildBC3D (const DescriptorList& d_list)
{
    BL_ASSERT(nsr > 0);
    delete [] bcr3D;
    bcr3D = new int[2*3*n_state]();
    int* bci = bcr3D;
    for (DeriveRec::StateRange* r = rng; r != nullptr; r = r->next)
    {
        const StateDescriptor& d = d_list[r->typ];

        for (int k = 0; k < r->nc; k++)
        {
            const int* bc = d.getBC(r->sc + k).vect();

            for (int j = 0; j < AMREX_SPACEDIM; j++)
            {
                bci[j] = bc[j];
            }
            bci += 3;
            for (int j = 0; j < AMREX_SPACEDIM; j++)
            {
                bci[j] = bc[AMREX_SPACEDIM+j];
            }
            bci += 3;
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

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveFunc              der_func,
                 const DeriveRec::DeriveBoxMap& bx_map,
                 Interpolater*           interp)
{
    lst.emplace_back(name,result_type,nvar_der,der_func,bx_map,interp);
}

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 DeriveFunc3D            der_func_3d,
                 const DeriveRec::DeriveBoxMap& bx_map,
                 Interpolater*           interp)
{
    lst.emplace_back(name,result_type,nvar_der,der_func_3d,bx_map,interp);
}

void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 const DeriveFuncFab&           der_func_fab,
                 const DeriveRec::DeriveBoxMap& bx_map,
                 Interpolater*           interp)
{
    lst.emplace_back(name,result_type,nvar_der,der_func_fab,bx_map,interp);
}

// This version doesn't take a Fortran function name, it is entirely defined by the C++
void
DeriveList::add (const std::string&      name,
                 IndexType               result_type,
                 int                     nvar_der,
                 const DeriveRec::DeriveBoxMap& box_map)
{
    lst.emplace_back(name,result_type,nvar_der,box_map);
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Vector<std::string> const&    vars,
                 DeriveFunc              der_func,
                 const DeriveRec::DeriveBoxMap& bx_map,
                 Interpolater*           interp)
{
    lst.emplace_back(name,res_typ,nvar_der,vars,der_func,bx_map,interp);
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Vector<std::string> const&    vars,
                 DeriveFunc3D            der_func_3d,
                 const DeriveRec::DeriveBoxMap& bx_map,
                 Interpolater*           interp)
{
    lst.emplace_back(name,res_typ,nvar_der,vars,der_func_3d,bx_map,interp);
}

void
DeriveList::add (const std::string&      name,
                 IndexType               res_typ,
                 int                     nvar_der,
                 Vector<std::string> const&    vars,
                 const DeriveFuncFab&           der_func_fab,
                 const DeriveRec::DeriveBoxMap& bx_map,
                 Interpolater*           interp)
{
    lst.emplace_back(name,res_typ,nvar_der,vars,der_func_fab,bx_map,interp);
}

std::list<DeriveRec>&
DeriveList::dlist ()
{
    return lst;
}

bool
DeriveList::canDerive (const std::string& name) const
{
    for (auto const& li : lst)
    {
        // Can be either a component name ...
        for (int i = 0; i < li.numDerive(); i++) {
           if (li.variableName(i) == name)
               return true;
        }
        // ... or a derive name
        if (li.derive_name == name)
            return true;
    }
    return false;
}

const DeriveRec*
DeriveList::get (const std::string& name) const
{
    for (auto const& li : lst)
    {
        // Can be either a component name ...
        for (int i = 0; i < li.numDerive(); i++) {
           if (li.variableName(i) == name)
               return &(li);
        }
        // ... or a derive name
        if (li.derive_name == name)
            return &(li);
    }
    return nullptr;
}

void
DeriveList::addComponent (const std::string&    name,
                          const DescriptorList& d_list,
                          int                   state_indx,
                          int                   s_comp,
                          int                   n_comp)
{
    auto li = lst.begin(), End = lst.end();

    for ( ; li != End; ++li)
    {
        if (li->derive_name == name) {
            break;
        }
    }

    BL_ASSERT (li != End);

    li->addRange(d_list, state_indx, s_comp, n_comp);
}

}
