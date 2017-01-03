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

#include <cmath>

#include "AMReX_SphereIF.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"

int checkSphere()
{
  int eekflag =  0;
#if 0
  //check to see that the sum of the volume fractions
  //comes to close to exactly the total volume
  //First: calculate volume of domain

  const IntVect& ivSize = a_domain.size();
  RealVect hiCorn;
  RealVect domLen;
  Real cellVolume = 1.0;
  Real totalVolume = 1.0;

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      hiCorn[idir] = a_origin[idir] + a_dx*Real(ivSize[idir]);
      cellVolume *= a_dx;
      domLen[idir] = hiCorn[idir] - a_origin[idir];
      totalVolume *= domLen[idir];
    }

  // Calculate the exact volume of the regular region
  Real volExact;
  Real bndryExact;

  if (SpaceDim == 2)
  {
    volExact = acos(-1.0) * a_radius*a_radius;
    bndryExact = 2*acos(-1.0) * a_radius;
  }

  if (SpaceDim == 3)
  {
    volExact = 4.0/3.0 * acos(-1.0) * a_radius*a_radius*a_radius;
    bndryExact = 4.0 * acos(-1.0) * a_radius*a_radius;
  }

  if (a_insideRegular == false)
  {
    volExact = totalVolume - volExact;
  }

  // Now calculate the volume of approximate sphere
  Real tolerance = 0.01;
  Real volTotal = 0.0;
  Real bndryTotal = 0.0;

  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisBox = a_ebisl[dit()];
      for (BoxIterator bit(a_grids.get(dit())); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          Vector<VolIndex> vofs = ebisBox.getVoFs(iv);

          if (vofs.size() > 1)
            {
              eekflag = 2;
              return eekflag;
            }

          for (int ivof = 0; ivof < vofs.size(); ivof++)
            {
              Real volFrac = ebisBox.volFrac(vofs[ivof]);
              volTotal += volFrac * cellVolume;

              Real bndryArea = ebisBox.bndryArea(vofs[ivof]);
              bndryTotal += bndryArea * cellVolume / a_dx;

              if (volFrac > tolerance && volFrac < 1.0 - tolerance)
                {
                  if (!ebisBox.isIrregular(iv))
                    {
                      eekflag = 4;
                      return eekflag;
                    }
                }
            }
        }
    }
//                   RealVect normal = ebisBox.normal(vofs[ivof]);
//                   for (int idir = 0; idir < SpaceDim; idir++)
//                     {
//                       if (Abs(normal[idir] - correctNorm[idir]) > tolerance)
//                         {
//                           eekflag = 5;
//                           return eekflag;
//                         }
//                     }
//                   RealVect centroid = ebisBox.centroid(vofs[ivof]);
//                   for (int idir = 0; idir < SpaceDim; idir++)
//                     {
//                       if (Abs(centroid[idir]) > 0.5)
//                         {
//                           eekflag = 6;
//                           return eekflag;
//                         }
//                     }
  //check how close is the answer
  pout() << endl;

#ifdef CH_MPI
  Real vt = volTotal;
  MPI_Allreduce(&vt, &volTotal, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  Real bt = bndryTotal;
  MPI_Allreduce(&bt, &bndryTotal, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
#endif
  pout() << "volTotal = "<< volTotal << endl;
  pout() << "volExact = "<< volExact << endl;
  pout() << "error = "   << volTotal - volExact << endl;
  pout() << endl;

  pout() << "bndryTotal = "<< bndryTotal << endl;
  pout() << "bndryExact = "<< bndryExact << endl;
  pout() << "error = "   << bndryTotal - bndryExact << endl;
  pout() << endl;

  if (Abs(volTotal - volExact) / volExact > tolerance)
  {
    MayDay::Error("Approx. volume not within tolerance");
  }

  if (Abs(bndryTotal - bndryExact) / bndryExact > tolerance)
  {
    MayDay::Error("Approx. boundary area not within tolerance");
  }

#endif
  return eekflag;
}


using std::cout;
using std::endl;
int
main(int argc,char **argv)
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif

  const char* in_file = "sphere.inputs";
  //parse input file
  amrex::ParmParse pp;
  pp.Initialize(argc, argv,in_file);

  // check volume and surface area of approximate sphere
  int eekflag = checkSphere();
  if (eekflag != 0)
    {
      cout << "non zero eek detected = " << eekflag << endl;
      cout << "sphere test failed" << endl;
    }
  else
    {
      cout << "sphere test passed" << endl;
    }

  return 0;
}


