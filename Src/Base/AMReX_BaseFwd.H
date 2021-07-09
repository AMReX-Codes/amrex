#ifndef AMREX_BASE_FWD_H_
#define AMREX_BASE_FWD_H_

namespace amrex {

class MultiFab;
class iMultiFab;
template <class FAB> class FabArray;
template <typename T> class LayoutData;
class FabArrayBase;

class FArrayBox;
class IArrayBox;
template <class T> class BaseFab;
template <typename T> struct Array4;

template <class T, unsigned int N> struct GpuArray;

class BoxArray;
class BoxList;
class DistributionMapping;
class Geometry;

class Box;
class IntVect;
class IndexType;
struct Dim3;
struct XDim3;
class RealBox;

class MFIter;
struct MFItInfo;

template <class FAB> class FabFactory;

class AMReX;
class Arena;
class BCRec;
class ParmParse;
struct RandomEngine;

template <int N> struct ParserExecutor;
class Parser;

namespace Gpu {
    class Elixir;
}

}

#endif
