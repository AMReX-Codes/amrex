
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>

namespace amrex {

EBFArrayBox::EBFArrayBox () = default;

EBFArrayBox::EBFArrayBox (Arena* ar)
    : FArrayBox(ar)
{
}

EBFArrayBox::EBFArrayBox (const EBCellFlagFab& ebflag, const Box& bx, int ncomps, Arena* ar,
                          const EBFArrayBoxFactory* factory, int box_index)
    : FArrayBox(bx, ncomps, ar),
      m_ebcellflag(&ebflag),
      m_factory(factory),
      m_box_index(box_index)
{
    BL_ASSERT(ebflag.box().contains(amrex::enclosedCells(bx)));
    const Box& ccbx = amrex::enclosedCells(bx);
    m_type = ebflag.getType(ccbx);
}

EBFArrayBox::EBFArrayBox (EBFArrayBox const& rhs, MakeType make_type, int scomp, int ncomp)
    : FArrayBox(rhs, make_type, scomp, ncomp),
      m_ebcellflag(rhs.m_ebcellflag),
      m_factory(rhs.m_factory),
      m_box_index(rhs.m_box_index)
{
    m_type = rhs.m_type;
}

EBFArrayBox::~EBFArrayBox () = default;

const FArrayBox*
EBFArrayBox::getLevelSetData () const
{
    if (m_factory && m_box_index >= 0) {
        MultiFab const& mf = m_factory->getLevelSet();
        return &(mf[m_box_index]);
    } else {
        return nullptr;
    }
}

const FArrayBox*
EBFArrayBox::getVolFracData () const
{
    if (m_factory && m_box_index >= 0) {
        MultiFab const& mf = m_factory->getVolFrac();
        return &(mf[m_box_index]);
    } else {
        return nullptr;
    }
}

const FArrayBox*
EBFArrayBox::getCentroidData () const
{
    if (m_factory && m_box_index >= 0) {
        MultiCutFab const& mf = m_factory->getCentroid();
        if (mf.ok(m_box_index)) {
            return &(mf[m_box_index]);
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

const FArrayBox*
EBFArrayBox::getBndryCentData () const
{
    if (m_factory && m_box_index >= 0) {
        MultiCutFab const& mf = m_factory->getBndryCent();
        if (mf.ok(m_box_index)) {
            return &(mf[m_box_index]);
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

const FArrayBox*
EBFArrayBox::getBndryNormalData () const
{
    if (m_factory && m_box_index >= 0) {
        MultiCutFab const& mf = m_factory->getBndryNormal();
        if (mf.ok(m_box_index)) {
            return &(mf[m_box_index]);
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

const FArrayBox*
EBFArrayBox::getBndryAreaData () const
{
    if (m_factory && m_box_index >= 0) {
        MultiCutFab const& mf = m_factory->getBndryArea();
        if (mf.ok(m_box_index)) {
            return &(mf[m_box_index]);
        } else {
            return nullptr;
        }
    } else {
        return nullptr;
    }
}

Array<const FArrayBox*, AMREX_SPACEDIM>
EBFArrayBox::getAreaFracData () const
{
    if (m_factory && m_box_index >= 0) {
        Array<MultiCutFab const*, AMREX_SPACEDIM> const& mfs = m_factory->getAreaFrac();
        if (mfs[0]->ok(m_box_index)) {
            return {AMREX_D_DECL(&((*mfs[0])[m_box_index]),
                                 &((*mfs[1])[m_box_index]),
                                 &((*mfs[2])[m_box_index]))};
        } else {
            return {AMREX_D_DECL(nullptr,nullptr,nullptr)};
        }
    } else {
        return {AMREX_D_DECL(nullptr,nullptr,nullptr)};
    }
}

Array<const FArrayBox*, AMREX_SPACEDIM>
EBFArrayBox::getFaceCentData () const
{
    if (m_factory && m_box_index >= 0) {
        Array<MultiCutFab const*, AMREX_SPACEDIM> const& mfs = m_factory->getFaceCent();
        if (mfs[0]->ok(m_box_index)) {
            return {AMREX_D_DECL(&((*mfs[0])[m_box_index]),
                                 &((*mfs[1])[m_box_index]),
                                 &((*mfs[2])[m_box_index]))};
        } else {
            return {AMREX_D_DECL(nullptr,nullptr,nullptr)};
        }
    } else {
        return {AMREX_D_DECL(nullptr,nullptr,nullptr)};
    }
}

Array<const FArrayBox*, AMREX_SPACEDIM>
EBFArrayBox::getEdgeCentData () const
{
    if (m_factory && m_box_index >= 0) {
        Array<MultiCutFab const*, AMREX_SPACEDIM> const& mfs = m_factory->getEdgeCent();
        if (mfs[0]->ok(m_box_index)) {
            return {AMREX_D_DECL(&((*mfs[0])[m_box_index]),
                                 &((*mfs[1])[m_box_index]),
                                 &((*mfs[2])[m_box_index]))};
        } else {
            return {AMREX_D_DECL(nullptr,nullptr,nullptr)};
        }
    } else {
        return {AMREX_D_DECL(nullptr,nullptr,nullptr)};
    }
}

const EBCellFlagFab&
getEBCellFlagFab (const FArrayBox& fab)
{
    const auto* ebfab = static_cast<EBFArrayBox const*>(&fab);
    BL_ASSERT(ebfab);
    return ebfab->getEBCellFlagFab();
}

}
