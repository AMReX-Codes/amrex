//BL_COPYRIGHT_NOTICE

//
// $Id: StateDescriptor.cpp,v 1.2 1997-11-24 18:52:33 lijewski Exp $
//

#include <StateDescriptor.H>
#include <Interpolater.H>
#include <BCRec.H>

// -------------------------------------------------------------
DescriptorList::DescriptorList()
    : desc(PArrayManage)
{
}

// -------------------------------------------------------------
DescriptorList::~DescriptorList()
{
}

// -------------------------------------------------------------
void
DescriptorList::clear()
{
    desc.clear();
}

// -------------------------------------------------------------
void
DescriptorList::addDescriptor(int indx, IndexType typ, TimeCenter ttyp,
			      int nextra, int num_comp, 
			      Interpolater *interp)
{
    if (indx >= desc.length())
	desc.resize(indx+1);
    desc.set(indx,new StateDescriptor(typ,ttyp,indx,nextra,
				     num_comp,interp));
}  

// -------------------------------------------------------------
void
DescriptorList::setComponent(int indx, int comp, const aString &nm,
                                  const BCRec& bc, BndryFunc func)
{
    desc[indx].setComponent(comp,nm,bc,func);
}  

// -------------------------------------------------------------
void
DescriptorList::resetComponentBCs(int indx, int comp,
                                  const BCRec& bc, BndryFunc func)
{
    desc[indx].resetComponentBCs(comp,bc,func);
}

// -------------------------------------------------------------
// -------------------------------------------------------------
StateDescriptor::StateDescriptor()
    : id(-1), ncomp(0), ngrow(0),
      mapper(0), names(), bc(), bc_func(), t_type(Point)
{
}

// -------------------------------------------------------------
StateDescriptor::StateDescriptor(IndexType btyp, TimeCenter ttyp,
                                 int ident, int nextra, 
                                 int num_comp, Interpolater *interp)
    : type(btyp), t_type(ttyp),
      id(ident), ngrow(nextra),
      ncomp(num_comp), mapper(interp),
      names(), bc(), bc_func()
{
    assert (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);
}

// -------------------------------------------------------------
StateDescriptor::~StateDescriptor()
{
    mapper = 0;
}

// -------------------------------------------------------------
void
StateDescriptor::define(IndexType btyp,  TimeCenter ttyp,
                             int ident, int nextra,
                             int num_comp, Interpolater *interp) 
{
    type   = btyp;
    t_type = ttyp;
    id     = ident;
    ngrow  = nextra;
    ncomp  = num_comp;
    mapper = interp;

    assert (num_comp > 0);
   
    names.resize(num_comp);
    bc.resize(num_comp);
    bc_func.resize(num_comp);

}
		       
// -------------------------------------------------------------
void
StateDescriptor::setComponent(int comp, const aString &nm, const BCRec &bcr,
                                   BndryFunc func)
{
    assert(comp >= 0 && comp < ncomp && names[comp] == "");
    names[comp] = nm;
    bc_func[comp] = func;
    bc[comp] = bcr;
}

// -------------------------------------------------------------
void
StateDescriptor::resetComponentBCs(int comp, const BCRec &bcr,
                                   BndryFunc func)
{
    assert(comp >= 0 && comp < ncomp);
    bc_func[comp] = func;
    bc[comp] = bcr;
}

// -------------------------------------------------------------
void
StateDescriptor::dumpNames(ostream &os, int start_comp,
                                int num_comp) const
{
    assert(start_comp >= 0 && start_comp+num_comp <= ncomp);
    int k;
    for (k = 0; k < num_comp; k++) {
	os << names[start_comp+k] << ' ';
    }
}
