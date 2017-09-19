
#include <AMReX_EBLevel.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MultiFab.H>

#include <type_traits>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex 
{
  static const IntVect eblev_debiv(D_DECL(0,0,0));

  Box
  getLevelDomain (const MultiFab& mf)
  {
    const auto& eblevel = amrex::getEBLevel(mf);
    return eblevel.getDomain();
  }

  const EBLevel&
  getEBLevel (const MultiFab& mf)
  {
    const EBFArrayBoxFactory* factory = dynamic_cast<EBFArrayBoxFactory const*>(&(mf.Factory()));
    BL_ASSERT(factory);
    return factory->getEBLevel();
  }

  const FabArray<EBCellFlagFab>&
  getMultiEBCellFlagFab (const MultiFab& mf)
  {
    const auto& eblevel = amrex::getEBLevel(mf);
    return eblevel.getMultiEBCellFlagFab();    
  }

// const FabArray<EBFaceFlagFab>&
// getMultiEBFaceFlagFab (const MultiFab& mf)
// {
//     const auto& eblevel = amrex::getEBLevel(mf);
//     return eblevel.getMultiEBFaceFlagFab();    
// }

  EBLevel::EBLevel ()
  {
  }

  EBLevel::~EBLevel ()
  {
  }

  EBLevel::EBLevel (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng)
    : EBLevelGrid(ba,dm,domain,ng),
      m_cellflags(std::make_shared<FabArray<EBCellFlagFab> >(ba, dm, 1, ng)),
//      m_faceflags(std::make_shared<FabArray<EBFaceFlagFab> >(ba, dm, AMREX_SPACEDIM, ng)),
      m_ebisl(std::make_shared<EBISLayout>(EBLevelGrid::getEBISL()))
  {
    defineDoit(ba, dm, domain, ng);
  }


  void
  EBLevel::define (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng)
  {
    static_assert(sizeof(EBCellFlag) == 4, "sizeof EBCellFlag != 4");
    static_assert(std::is_standard_layout<EBCellFlag>::value == true, "EBCellFlag is not pod");

    static_assert(sizeof(EBFaceFlag) == 4, "sizeof EBFaceFlag != 4");
    static_assert(std::is_standard_layout<EBFaceFlag>::value == true, "EBFaceFlag is not pod");

    EBLevelGrid::define(ba, dm, domain, ng);

    m_cellflags = std::make_shared<FabArray<EBCellFlagFab> >(ba, dm, 1, ng);
    //      m_faceflags = std::make_shared<FabArray<EBFaceFlagFab> >(ba, dm, AMREX_SPACEDIM, ng);
    m_ebisl = std::make_shared<EBISLayout>(EBLevelGrid::getEBISL());

    defineDoit(ba, dm, domain, ng);
  }

  void 
  EBLevel::setIrregularVolFab (EBCellFlagFab & fab, const EBISBox& ebis, const Box& bx)
  {
    BL_PROFILE("EBLevel::setIrregularVoFab()");
    int nregular=0, nsingle=0, nmulti=0, ncovered=0;
    int ncells = bx.numPts();
    Box domain = ebis.getDomain();
    for (BoxIterator bi(bx); bi.ok(); ++bi)
    {
      const IntVect& iv = bi();
      auto& cellflag = fab(iv);
      if (ebis.isRegular(iv)) 
      {
        cellflag.setRegular();
        ++nregular;
      } 
      else if (ebis.isCovered(iv)) 
      {
        cellflag.setCovered();
        ++ncovered;
      } 
      else if (ebis.isMultiValued(iv)) 
      {
        cellflag.setMultiValued(ebis.numVoFs(iv));
        ++nmulti;
        amrex::Abort("EBLevel: multi-value cell not supported yet");
      } 
      else 
      {
        cellflag.setSingleValued();
        ++nsingle;
      }
    }
        
    if (nregular == ncells) {
      fab.setType(FabType::regular);
    } else if (ncovered == ncells) {
      fab.setType(FabType::covered);
    } else if (nmulti > 0) {
      fab.setType(FabType::multivalued);
    } else {
      fab.setType(FabType::singlevalued);
    }

    const Box& ibx = amrex::grow(fab.box(),-1) & domain;
    for (BoxIterator bi(ibx); bi.ok(); ++bi)
    {
      const IntVect& iv = bi();
      auto& cellflag = fab(iv);
      cellflag.setDisconnected();
      cellflag.setConnected(IntVect::TheZeroVector());

      const auto& vofs = ebis.getVoFs(iv);
      for (const auto& vi : vofs)
      {
        std::array<int,AMREX_SPACEDIM> dirs = {AMREX_D_DECL(0,1,2)};
        do
        {
          IntVect offset_0 = IntVect::TheZeroVector();
          for (SideIterator sit_0; sit_0.ok(); ++sit_0)
          {
            offset_0[dirs[0]] = amrex::sign(sit_0());

            const auto& vofs_0 = ebis.getVoFs(vi, dirs[0], sit_0(), 1);
            for (const auto& vi_0 : vofs_0)
            {
              cellflag.setConnected(offset_0);

#if (AMREX_SPACEDIM >= 2)
              IntVect offset_1 = offset_0;
              for (SideIterator sit_1; sit_1.ok(); ++sit_1)
              {
                offset_1[dirs[1]] = amrex::sign(sit_1());

                const auto& vofs_1 = ebis.getVoFs(vi_0, dirs[1], sit_1(), 1);
                for (const auto& vi_1 : vofs_1)
                {
                  cellflag.setConnected(offset_1);

#if (AMREX_SPACEDIM == 3)
                  IntVect offset_2 = offset_1;
                  for (SideIterator sit_2; sit_2.ok(); ++sit_2)
                  {
                    offset_2[dirs[2]] = amrex::sign(sit_2());

                    const auto& vofs_2 = ebis.getVoFs(vi_1, dirs[2], sit_2(), 1);
                    for (const auto& vi_2 : vofs_2)
                    {
                      cellflag.setConnected(offset_2);
                    }
                  }
#endif
                }
              }
#endif
            }
          }
        } while (std::next_permutation(dirs.begin(), dirs.end()));
      }
    } 
  }
///
  void 
  EBLevel::setIrregularFaceFab (EBFaceFlagFab & fabs, const EBISBox& ebis, const Box& ccbx)
  {
    BL_PROFILE("EBLevel::setIrregularFaceFab()");
    for (BoxIterator bi(ccbx); bi.ok(); ++bi)
    {
      const IntVect& iv = bi();
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
      {
        {
          const auto& lo_faces = ebis.getAllFaces(iv, dir, Side::Lo);
          auto& flag = fabs(dir,iv);
          flag.setDisconnected();
          flag.setConnected(dir, IntVect::TheZeroVector());
          if (lo_faces.size() == 0) {
            flag.setCovered();
          } else if (lo_faces.size() == 1) {
            if (ebis.areaFrac(lo_faces[0]) == 1.0) {
              flag.setRegular();
            } else {
              flag.setSingleValued();
            }
            setFaceConnection(flag, lo_faces[0], ebis);
          } else {
            flag.setMultiValued(lo_faces.size());
            amrex::Abort("EBLevel: multi-value face not supported yet");
          }
        }
                
        if (iv[dir] == ccbx.bigEnd(dir))
        {
          const IntVect& ivhi = iv + IntVect::TheDimensionVector(dir);
          const auto& hi_faces = ebis.getAllFaces(iv, dir, Side::Hi);
          auto& flag = fabs(dir,ivhi);
          flag.setDisconnected();
          flag.setConnected(dir, IntVect::TheZeroVector());
          if (hi_faces.size() == 0) {
            flag.setCovered();
          } else if (hi_faces.size() == 1) {
            if (ebis.areaFrac(hi_faces[0]) == 1.0) {
              flag.setRegular();
            } else {           
              flag.setSingleValued();
            }
            setFaceConnection(flag, hi_faces[0], ebis);
          } else {
            flag.setMultiValued(hi_faces.size());
            amrex::Abort("EBLevel: multi-value face not supported yet");
          }
        }
      }
    }
  }
///
  void
  EBLevel::defineDoit (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng)
  {
    BL_PROFILE("EBLevel::defineDoit()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

      for (MFIter mfi(*m_cellflags); mfi.isValid(); ++mfi)
      {
        auto& fab = (*m_cellflags)[mfi];
        const EBISBox& ebis = (*m_ebisl)[mfi];
        const Box& bx = fab.box() & domain;

        //change this to if 0 if you want this to work the old way
#if 1
        fab.copy(ebis.getEBGraph().getEBCellFlagFab());
        fab.setType(ebis.getEBGraph().getEBCellFlagFab().getType());
#else

        if(ebis.isRegular(bx))
        {
          EBCellFlag flag;
          flag.setRegular();
          fab.setVal(flag);
          fab.setType(FabType::regular);
        }
        else if(ebis.isCovered(bx))
        {
          EBCellFlag flag;
          flag.setCovered();
          fab.setVal(flag);
          fab.setType(FabType::covered);
        }
        else
        {
          setIrregularVolFab(fab, ebis, bx);
        }
#endif        
      }

      if (m_faceflags) 
      {
        for (MFIter mfi(*m_faceflags); mfi.isValid(); ++mfi)
        {
          const EBISBox& ebis = (*m_ebisl)[mfi];
          auto& fabs = (*m_faceflags)[mfi];
          const Box& ccbx = amrex::grow(fabs.box(),-1) & domain;
          if(ebis.isRegular(ccbx))
          {
            EBFaceFlag flag;
            flag.setRegular();
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              BaseFab<EBFaceFlag>& fab = fabs.getFaceFlagFab(idir);
              fab.setVal(flag);
            }
          }
          else if(ebis.isCovered(ccbx))
          {
            EBFaceFlag flag;
            flag.setCovered();
            for(int idir = 0; idir < SpaceDim; idir++)
            {
              BaseFab<EBFaceFlag>& fab = fabs.getFaceFlagFab(idir);
              fab.setVal(flag);
            }
          }
          else 
          {
            setIrregularFaceFab(fabs, ebis, ccbx);
          }
        }
      }
    }
  }

  void
  EBLevel::setFaceConnection (EBFaceFlag& flag, const FaceIndex& face, const EBISBox& ebis)
  {
    const auto& vof_lo = face.getVoF(Side::Lo);
    const auto& vof_hi = face.getVoF(Side::Hi);

    // at domain boundary
    if (vof_lo.cellIndex() < 0 || vof_hi.cellIndex() < 0) return;

    const int dir = face.direction();
    Array<int> tdirs;  // transverse directions
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (idim != dir) {
        tdirs.push_back(idim);
      }
    }

    do
    {
      IntVect offset_0 = IntVect::TheZeroVector();
      for (SideIterator sit_0; sit_0.ok(); ++sit_0)
      {
        offset_0[tdirs[0]] = amrex::sign(sit_0());
        const auto& vofs_0_lo = ebis.getVoFs(vof_lo, tdirs[0], sit_0(), 1);
        const auto& vofs_0_hi = ebis.getVoFs(vof_hi, tdirs[0], sit_0(), 1);
        if (vofs_0_lo.size() == vofs_0_hi.size() && vofs_0_lo.size() == 1)
        {
          if (ebis.isConnected(vofs_0_lo[0], vofs_0_hi[0])) {
            // wz. This is known to break in some cases.
            flag.setConnected(dir, offset_0);
#if (AMREX_SPACEDIM == 3)
            IntVect offset_1 = offset_0;
            for (SideIterator sit_1; sit_1.ok(); ++sit_1)
            {
              offset_1[tdirs[1]] = amrex::sign(sit_1());
              const auto& vofs_1_lo = ebis.getVoFs(vofs_0_lo[0], tdirs[1], sit_1(), 1);
              const auto& vofs_1_hi = ebis.getVoFs(vofs_0_hi[0], tdirs[1], sit_1(), 1);
              if (vofs_1_lo.size() == vofs_1_hi.size() && vofs_1_lo.size() == 1)
              {
                if (ebis.isConnected(vofs_1_lo[0], vofs_1_hi[0])) {
                  // wz. This is known to break in some cases.
                  flag.setConnected(dir, offset_1);
                }
              }
            }
#endif
          }
        }
      }

    } while (std::next_permutation(tdirs.begin(), tdirs.end()));
  }

}
