
#include <AMReX_EBFArrayBox.H>

namespace amrex {

EBFArrayBox::EBFArrayBox ()
{
}

EBFArrayBox::EBFArrayBox (const EBISBox& ebisBox, const Box& box, int ncomps)
    : FArrayBox(box, ncomps),
      m_ebisbox(ebisBox)
{
    if (ebisBox.isRegular(box))
    {
        m_type = FabType::regular;
        return;  // nothing to do
    }
    else if (ebisBox.isCovered(box))
    {
        m_type = FabType::allcovered;
        return;
    }
    else if (ebisBox.getMultiCells(box).isEmpty())
    {
        m_type = FabType::singlevalue;
        return;
    }
    else
    {
        m_type = FabType::multivalue;
        return;
    }
}

EBFArrayBox::~EBFArrayBox ()
{

}

}
