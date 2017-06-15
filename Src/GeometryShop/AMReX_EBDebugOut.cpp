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
  void dumpVVoFs(const vector<VolIndex>* a_vofs)
  {
    const vector<VolIndex>& vofs = *a_vofs;
    pout() << "vector<volindex> contains:" << endl;
    for(int ivof = 0; ivof < vofs.size(); ivof++)
    {
      amrex::Print() << vofs[ivof] << std::endl;
    }
  }
  void dumpVFaces(const vector<FaceIndex>* a_faces)
  {
    const vector<FaceIndex>& faces = *a_faces;
    pout() << "vector<FaceIndex> contains:" << endl;
    for(int iface = 0; iface < faces.size(); iface++)
    {
      amrex::Print() << iface << ":" << faces[iface] << endl;
    }
  }
  void dumpDBL(const BoxArray* a_fabPtr)
  {
    amrex::Print() << "BoxArray contains:"  << *a_fabPtr << endl;
    
  }
  void dumpBA(const BoxArray* a_fabPtr)
  {
    amrex::Print() << "BoxArray contains:"  << *a_fabPtr << endl;
    
  }
  void dumpFAB(const FArrayBox* a_fabPtr)
  {
    const FArrayBox& fab = *a_fabPtr;
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      amrex::Print() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        amrex::Print() << "\t" << fab(bit(),ivar);
      }
      amrex::Print() << "\n";
    }
  }

  void dumpBFR(const BaseFab<Real>* a_fabPtr)
  {
    const BaseFab<Real>& fab = *a_fabPtr;
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      amrex::Print() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        amrex::Print() << "\t" << fab(bit(),ivar);
      }
      amrex::Print() << "\n";
    }
  }
  void dumpBFI(const BaseFab<int>* a_fabPtr)
  {
    const BaseFab<int>& fab = *a_fabPtr;
    BoxIterator bit(fab.box());
    for (bit.reset(); bit.ok(); ++bit)
    {
      amrex::Print() << "\t" << bit() ;
      for (int ivar = 0; ivar < fab.nComp(); ivar++)
      {
        amrex::Print() << "\t" << fab(bit(),ivar);
      }
      amrex::Print() << "\n";
    }
  }

  void dumpBL(const BoxArray* a_dblInPtr)
  {
    amrex::Print() << "BoxArray: ";
    for(int ibox = 0; ibox < a_dblInPtr->size(); ibox++)
    {
      amrex::Print() << (*a_dblInPtr)[ibox] << "   ";
    }
    amrex::Print() << "\n";
  }

  void dumpIVS(const IntVectSet* a_ivsPtr)
  {
    const IntVectSet& ivs = *a_ivsPtr;
    IVSIterator it(ivs);

    amrex::Print() << ": IntVects in the IVS are:" << "\n";

    for (it.begin(); it.ok(); ++it)
    {
      amrex::Print() << it() << "  ";
    }

    amrex::Print() << "\n";
  }

  void dumpBox(const Box* a_boxPtr)
  {
    amrex::Print() << "Box:" << *a_boxPtr << "\n";
  }


  void dumpEBFAB(const EBCellFAB* a_fab)
  {
    const EBCellFAB& fab = *a_fab;
    Box box = fab.getRegion();
    IntVectSet ivs(box);
    amrex::Print() << "valid and ghost data in ebcellfab" << "\n";

    const EBGraph& ebgraph = a_fab->getEBISBox().getEBGraph();
    const int ncomp = a_fab->nComp();
    for (VoFIterator vofit(ivs, ebgraph); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      amrex::Print() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();

      amrex::Print() << ";      data=";
      for (int ivar = 0; ivar < ncomp; ivar++)
      {
        amrex::Print() << " "
                       << setprecision(8)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << (*a_fab)(vof, ivar);
      }
      amrex::Print() << "\n";
    }
  }


  void dumpIVFAB(const BaseIVFAB<Real>* a_vectPtr)
  {
 
    if(a_vectPtr->numVoFs()==0) 
    {
      amrex::Print()<<"empty ";
      return;
    }
    const std::vector<VolIndex>& vofs = a_vectPtr->getVoFs();
    for(int i=0; i<vofs.size(); i++)
    {
      const VolIndex& vof=vofs[i];
      amrex::Print() << "vof= " << vof.gridIndex() << ", " << vof.cellIndex();
      amrex::Print() << ";      data=";
      for (int ivar = 0; ivar < a_vectPtr->nComp(); ivar++)
      {
        amrex::Print() << " "
                       << setprecision(8)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << a_vectPtr->operator()(vof, ivar);
      }
    }
    amrex::Print() << "\n";
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
      amrex::Print() << "face= (("
                     << voflo.gridIndex() << ", " << voflo.cellIndex() << "), ("
                     << vofhi.gridIndex() << ", " << vofhi.cellIndex() << "))" ;

      amrex::Print() << ";      data=";
      for (int ivar = 0; ivar < ncomp; ivar++)
      {
        amrex::Print() << " "
                       << setprecision(8)
                       << setiosflags(ios::showpoint)
                       << setiosflags(ios::scientific)
                       << (*a_fab)(face, ivar);
      }
      amrex::Print() << "\n";
    }
  }

}
