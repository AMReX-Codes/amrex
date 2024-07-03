#ifndef AMREX_BASE_FWD_H_
#define AMREX_BASE_FWD_H_
#include <AMReX_Config.H>

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

template<int dim>
class BoxND;
using Box = BoxND<AMREX_SPACEDIM>;
template<int dim>
class IntVectND;
using IntVect = IntVectND<AMREX_SPACEDIM>;
template<int dim>
class IndexTypeND;
using IndexType = IndexTypeND<AMREX_SPACEDIM>;
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
