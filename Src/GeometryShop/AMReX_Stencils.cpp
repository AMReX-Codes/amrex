#include "AMReX_Stencils.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_parstream.H"

namespace amrex
{
  ///
  void 
  VoFStencil::
  outputToPout() const
  {
    pout() << "num vof \t  weight \t variable = " << endl;
    for(int ivof = 0; ivof < this->size(); ivof++)
    {
      pout() << "(" << ivof << "):" << vofs[ivof]  << "\t" << weights[ivof]<< "\t" << variables[ivof] << endl;
    }
  }

  ///
  void 
  FaceStencil::
  outputToPout() const
  {
    pout() << " face \t  weight = " << endl;
    for(int ivof = 0; ivof < this->size(); ivof++)
    {
      pout() << faces[ivof]  << "\t" << weights[ivof]<< endl;
    }
  }


  /**************/
  int
  VoFStencil::size() const
  {
    return weights.size();
  }
  /**************/
  /**************/
  const  VolIndex&
  VoFStencil::vof(int i) const
  {
    return vofs[i];
  }
  /**************/
  /**************/
  const Real&
  VoFStencil::weight(int i) const
  {
    return weights[i];
  }

  Real&
  VoFStencil::weight(int i)
  {
    return weights[i];
  }

  Real applyVoFStencil(const VoFStencil& a_sten, const EBCellFAB& a_fab)
  {
    Real retval = 0.;
    for (int isten = 0; isten < a_sten.size(); isten++)
      {
        retval += (a_sten.weight(isten))*(a_fab((a_sten.vof(isten)), a_sten.variable(isten)));
      }

    return retval;
  }
  ///
  /**
     shift all entries by a_shift
  */
  ///
  void 
  VoFStencil:: 
  shift(const IntVect& a_shift)
  {
    for(int ivof = 0; ivof < this->size(); ivof++)
      {
        VolIndex& vf = vofs[ivof];
        vf.shift(a_shift);
      }
  }

  ///
  void 
  FaceStencil:: 
  shift(const IntVect& a_shift)
  {
    for(int iface = 0; iface < this->size(); iface++)
      {
        FaceIndex& f = faces[iface];
        f.shift(a_shift);
      }

  }


  /**************/
  /**************/
  VoFStencil::VoFStencil()
    : vofs(0),
      weights(0)
  {
  }
  /**************/
  /**************/
  void
  VoFStencil::clear()
  {
    vofs.resize(0);
    variables.resize(0);
    weights.resize(0);
  }
  /**************/
  /**************/
  void
  FaceStencil::clear()
  {
    faces.resize(0);
    weights.resize(0);
    variables.resize(0);
  }

  /**************/
  /**************/
  VoFStencil::VoFStencil(const VoFStencil&  stenin)
    : vofs(stenin.vofs),
      weights(stenin.weights),
      variables(stenin.variables)
  {
  }
  /**************/
  /**************/
  VoFStencil::~VoFStencil()
  {
  }
  /**************/
  /**************/
  void
  VoFStencil::add(const VolIndex& a_vof,const Real& wgt, int ivar)
  {
    bool alreadyhere = false;
    for (int ivof = 0; ivof < vofs.size(); ivof++)
      {
        if ((vofs[ivof] == a_vof) && (variables[ivof] == ivar))
          {
            alreadyhere = true;
            weights[ivof] += wgt;
          }
      }
    if (!alreadyhere)
      {
        vofs.push_back(a_vof);
        weights.push_back(wgt);
        variables.push_back(ivar);
      }
  }
  /**************/
  /**************/
  VoFStencil&
  VoFStencil::operator+=(const VoFStencil& vofstenin)
  {
    for (int ivof = 0; ivof < vofstenin.size(); ivof++)
      {
        add(vofstenin.vof(ivof), vofstenin.weight(ivof), vofstenin.variable(ivof));
      }
    return *this;
  }
  /**************/
  /**************/
  void
  VoFStencil::operator*=(const Real& a_scaling)
  {
    for (int ivof = 0; ivof < size(); ivof++)
      {
        weights[ivof] *= a_scaling;
      }
  }
  /**************/
  /**************/
  void
  FaceStencil::operator*=(const Real& a_scaling)
  {
    for (int iface = 0; iface < size(); iface++)
      {
        weights[iface] *= a_scaling;
      }
  }
  /**************/
  /**************/
  VoFStencil&
  VoFStencil::operator=(const VoFStencil& vofstenin)
  {
    clear();
    this->operator+=(vofstenin);
    return *this;
  }
  /**************/

  /**************/
  /**************/
  FaceStencil::FaceStencil()
    : faces(0),
      weights(0)
  {
  }
  /**************/
  /**************/
  FaceStencil::FaceStencil(const FaceStencil& facestenin)
    : faces(facestenin.faces),
      weights(facestenin.weights),
      variables(facestenin.variables)
  {
  }
  /**************/
  /**************/
  // destructor
  FaceStencil::~FaceStencil()
  {
  }
  /**************/
  /**************/
  void
  FaceStencil::add(const FaceIndex& a_face,const Real& a_weight, int ivar)
  {
    bool alreadyhere = false;
    for (int iface = 0; iface < faces.size(); iface++)
      {
        if ( (faces[iface] == a_face) && (variables[iface] == ivar))
          {
            alreadyhere = true;
            weights[iface] += a_weight;
          }
      }
    if (!alreadyhere)
      {
        faces.push_back(a_face);
        weights.push_back(a_weight);
        variables.push_back(ivar);
      }
  }
  /**************/
  /**************/
  int
  FaceStencil::size() const
  {
    return weights.size();
  }
  /**************/
  /**************/
  const FaceIndex&
  FaceStencil::face(int i) const
  {
    return faces[i];
  }
  /**************/
  /**************/
  const Real&
  FaceStencil::weight(int i) const
  {
    return weights[i];
  }
  const int&
  FaceStencil::variable(int i) const
  {
    return variables[i];
  }
  int&
  FaceStencil::variable(int i)
  {
    return variables[i];
  }
  const int&
  VoFStencil::variable(int i) const
  {
    return variables[i];
  }
  int&
  VoFStencil::variable(int i)
  {
    return variables[i];
  }
  /**************/
  /**************/
  FaceStencil&
  FaceStencil::operator+=(const FaceStencil& Facestenin)
  {
    for (int iFace = 0; iFace < Facestenin.size(); iFace++)
      {
        add(Facestenin.face(iFace), Facestenin.weight(iFace));
      }
    return *this;
  }
  /**************/
  /**************/
  FaceStencil&
  FaceStencil::operator=(const FaceStencil& Facestenin)
  {
    clear();
    this->operator+=(Facestenin);
    return *this;
  }
}

