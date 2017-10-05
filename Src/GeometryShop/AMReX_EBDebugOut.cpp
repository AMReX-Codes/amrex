/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "AMReX_EBDebugOut.H"
#include "AMReX_Print.H"
#include "AMReX_BoxArray.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_parstream.H"
#include <iomanip>


namespace amrex
{

  void dumpVVoFs(const Vector<VolIndex>* a_vofs)
  {
    const Vector<VolIndex>& vofs = *a_vofs;
    pout() << "Vector<volindex> contains:" << endl;
    for(int ivof = 0; ivof < vofs.size(); ivof++)
    {
      pout() << vofs[ivof] << std::endl;
    }
  }
  void dumpVFaces(const Vector<FaceIndex>* a_faces)
  {
    const Vector<FaceIndex>& faces = *a_faces;
    pout() << "Vector<FaceIndex> contains:" << endl;
    for(int iface = 0; iface < faces.size(); iface++)
    {
      pout() << iface << ":" << faces[iface] << endl;
    }
  }
  void dumpDBL(const BoxArray* a_fabPtr)
  {
    pout() << "BoxArray contains:"  << *a_fabPtr << endl;
    
  }
  void dumpBA(const BoxArray* a_fabPtr)
  {
    pout() << "BoxArray contains:"  << *a_fabPtr << endl;
    
  }
  void dumpFAB(const FArrayBox* a_fabPtr)
  {
    const FArrayBox& fab = *a_fabPtr;
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        pout() << "\t" << fab(bit(),ivar);
      }
      pout() << "\n";
    }
  }

  void dumpBFR(const BaseFab<Real>* a_fabPtr)
  {
    const BaseFab<Real>& fab = *a_fabPtr;
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        pout() << "\t" << fab(bit(),ivar);
      }
      pout() << "\n";
    }
  }
  void dumpBFI(const BaseFab<int>* a_fabPtr)
  {
    const BaseFab<int>& fab = *a_fabPtr;
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      pout() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        pout() << "\t" << fab(bit(),ivar);
      }
      pout() << "\n";
    }
  }

  void dumpBL(const BoxArray* a_dblInPtr)
  {
    pout() << "BoxArray: ";
    for(int ibox = 0; ibox < a_dblInPtr->size(); ibox++)
    {
      pout() << (*a_dblInPtr)[ibox] << "   ";
    }
    pout() << "\n";
  }

  void dumpIVS(const IntVectSet* a_ivsPtr)
  {
    const IntVectSet& ivs = *a_ivsPtr;
    IVSIterator it(ivs);

    pout() << ": IntVects in the IVS are:" << "\n";

    for (it.begin(); it.ok(); ++it)
    {
      pout() << it() << "  ";
    }

    pout() << "\n";
  }

  void dumpBox(const Box* a_boxPtr)
  {
    pout() << "Box:" << *a_boxPtr << "\n";
  }


  void dumpEBLevel(const FabArray<EBCellFAB>* a_fab)
  {
    pout() << "data for eb level = " << endl;
    for(MFIter mfi(*a_fab); mfi.isValid(); ++mfi)
    {
      dumpEBFAB(&(*a_fab)[mfi]);
    }
  }
    
  void dumpEBFAB(const EBCellFAB* a_fab)
  {
    const EBCellFAB& fab = *a_fab;
    Box box = fab.getRegion();
    IntVectSet ivs(box);
    pout() << "valid and ghost data in ebcellfab over box " << box << endl;

    const EBGraph& ebgraph = a_fab->getEBISBox().getEBGraph();
    const int ncomp = a_fab->nComp();
    for(BoxIterator bit(box); bit.ok(); ++bit)
    {
      vector<VolIndex> vofs = ebgraph.getVoFs(bit());
      for(int ivof = 0; ivof < vofs.size(); ivof++)
      {
        const VolIndex& vof = vofs[ivof];
        pout() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();

        pout() << ";      data=";
        for (int ivar = 0; ivar < ncomp; ivar++)
        {
          pout() << " "
                 << setprecision(8)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << (*a_fab)(vof, ivar);
        }
        pout() << endl;
      }
    }
  }

  void dumpFaceSten(const FaceStencil* a_vectPtr)
  {
    if (a_vectPtr == NULL) return;
    const FaceStencil& sten = *a_vectPtr;
    pout() << "stencil contains:" << endl;
    for (int isten = 0; isten <  sten.size(); isten++)
    {
      const FaceIndex& face = sten.face(isten);
      const Real& weight = sten.weight(isten);
      pout() << "face: ";
      for (SideIterator sit; sit.ok(); ++sit)
      {
        VolIndex vof = face.getVoF(sit());
        pout() <<  vof.gridIndex() << "   " << vof.cellIndex() << "  ";;
      }
      pout() << "weight: " << weight;
      pout() << endl;
    }
  }

  void dumpVoFSten(const VoFStencil* a_vectPtr)
  {
    if (a_vectPtr == NULL) return;
    const VoFStencil& sten = *a_vectPtr;
    pout() << "stencil contains:" << endl;
    for (int isten = 0; isten <  sten.size(); isten++)
    {
      const VolIndex& vof = sten.vof(isten);
      const Real& weight = sten.weight(isten);
      pout() << "vof: (";
      pout() <<  vof.gridIndex() << ", " << vof.cellIndex() << ") ";;
      pout() << "var: " << sten.variable(isten) << " ";
      pout() << "weight: " << weight;
      pout() << endl;
    }
  }

  void dumpIVFAB(const BaseIVFAB<Real>* a_vectPtr)
  {
 
    if(a_vectPtr->numVoFs()==0) 
    {
      pout()<<"empty ";
      return;
    }
    const Vector<VolIndex>& vofs = a_vectPtr->getVoFs();
    for(int i=0; i<vofs.size(); i++)
    {
      const VolIndex& vof=vofs[i];
      pout() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();
      pout() << ";      data=";
      for (int ivar = 0; ivar < a_vectPtr->nComp(); ivar++)
      {
        pout() << " "
               << setprecision(8)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << a_vectPtr->operator()(vof, ivar);
      }
    }
    pout() << "\n";
  }
  void dumpEBFace(const EBFaceFAB* a_fab)
  {
    int idir = a_fab->direction();
    const EBGraph& ebgraph = a_fab->getEBISBox().getEBGraph();
    const int ncomp = a_fab->nComp();
    const EBFaceFAB& fab = *a_fab;
    Box box = fab.getCellRegion();
    box &=  fab.getEBISBox().getDomain();
    IntVectSet ivs(box);

    for (FaceIterator faceit(ivs, ebgraph, idir, FaceStop::SurroundingWithBoundary); faceit.ok(); ++faceit)
    {
      const FaceIndex& face = faceit();

      const VolIndex& voflo = face.getVoF(Side::Lo);
      const VolIndex& vofhi = face.getVoF(Side::Hi);
      pout() << "face= (("
             << voflo.gridIndex() << ", " << voflo.cellIndex() << "), ("
             << vofhi.gridIndex() << ", " << vofhi.cellIndex() << "))" ;

      pout() << ";      data=";
      for (int ivar = 0; ivar < ncomp; ivar++)
      {
        pout() << " "
               << setprecision(8)
               << setiosflags(ios::showpoint)
               << setiosflags(ios::scientific)
               << (*a_fab)(face, ivar);
      }
      pout() << "\n";
    }
  }

}
