#ifndef _MFTower_
#define _MFTower_


#include <winstd.H>

#include <BCRec.H>
#include <Layout.H>

struct MFTower
{
  MFTower(const Layout&    layout,
          const IndexType& t,
          int              nComp,
          int              nGrow,
          int              numLevs);

  MFTower(Layout&           _layout,
          PArray<MultiFab>& pamf,
          int               _nLevs,
          int               _bComp,
          int               _nComp);

  int NComp() const;
  int NGrow() const;
  int BaseComp() const;
  int NumLevels() const;
  MultiFab& operator[](int i);
  const MultiFab& operator[](int i) const;
  bool IsCompatible(const MFTower& rhs,
                    int            bComp,
                    int            nComp) const;
  Real norm(int numLevs, int comp) const;
  const Layout& GetLayout() const;
  const IndexType& ixType() const;

  void Write(const std::string& fileName) const;
 
  static void CCtoECgrad(PArray<MFTower>& mft_ec,
                         const MFTower&   mft_cc,
                         Real             mult,
                         int              sComp,
                         int              dComp,
                         int              nComp,
                         int              numLevs);
    
  static void ECtoCCdiv(MFTower&               mft_cc,
                        const PArray<MFTower>& mft_ec,
                        Real                   mult,
                        int                    sComp,
                        int                    dComp,
                        int                    nComp,
                        int                    numLevs);
    
  static void AverageDown(MFTower& mft,
                          int      sComp,
                          int      nComp,
                          int      numLevs);
    
  void SetVal(Real     val,
              int      sComp,
              int      nComp,
              int      numLevs);
  
  static IndexType EC[BL_SPACEDIM], CC, NC;

protected:
  void define_alloc();
  void define_noalloc(PArray<MultiFab>& pamf);
  const Layout& layout;
  PArray<MultiFab> mft;
  IndexType iType;
  int nComp, nGrow, nLevs, bComp;
};


struct MFTFillPatch
{
  MFTFillPatch(Layout& layout);
  int NumLevels() const;
  Layout& GetLayout();

  void BuildStencil(const BCRec& bc,
                    int          maxorder);

  void DoCoarseFineParallelInterp(MFTower& mft,
                                  int      sComp,
                                  int      nComp,
                                  int      numLevs) const;
  void FillGrowCells(MFTower& mft,
                     int      sComp,
                     int      nComp,
                     bool     do_piecewise_constant,
                     int      numLevs) const;

  void FillGrowCellsSimple(MFTower& mft,
                           int      sComp,
                           int      nComp,
                           int      numLevs) const;

  const Array<Array<IVSMap> >& GrowCellStencil() const {return perpInterpStencil;} //[lev][dir][iv] = stencil to fill grow cells

  struct MyHash // Unique characteristics
  {
    MyHash()
      : maxorder(-1) {}
    MyHash(const BCRec& _bc, int maxorder);
    BCRec bc;
    int maxorder;
  };

protected:
  void BuildCFParallelInterpStencil();


  MyHash myhash;
  int nLevs;
  Layout& layout;
  Array<Array<IVSMap> > parallelInterpStencil; //[lev][dir][iv] = stencil to perform parallel bc Interp
  Array<Array<IVSMap> > perpInterpStencil; //[lev][dir][iv] = stencil to perform perpendicular bc Interp
};

inline bool operator==(const MFTFillPatch::MyHash& lhs, const MFTFillPatch::MyHash& rhs);
inline bool operator!=(const MFTFillPatch::MyHash& lhs, const MFTFillPatch::MyHash& rhs);

#endif
