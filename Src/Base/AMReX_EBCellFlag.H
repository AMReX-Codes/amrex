#ifndef AMREX_EBCELLFLAG_H_
#define AMREX_EBCELLFLAG_H_

namespace amrex {

class EBCellFlag
{
public:
    EBCellFlag () = default;
    EBCellFlag (const EBCellFlag& rhs) noexcept = default;
    EBCellFlag& operator= (const EBCellFlag& rhs) noexcept = default;
    EBCellFlag& operator= (int i) { flag = i; return *this; }
private:
    int flag = 0;
};

}

#endif
