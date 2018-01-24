
#ifndef BL_LOOPING_H
#define BL_LOOPING_H

#if (BL_SPACEDIM == 1)


/*
  The macro ForAllThisCBNN(T,b,ns,nc) is intended to facilitate efficient
  looping over the contents of BaseFabs and objects derived from BaseFab.
  Special attention has been paid to make it work efficiently on vector
  supercomputers.

  This is the constant version of ForAllThisBNN; i.e. it works when the
  underlying BaseFab is constant.

  Looping macro mnemonics:

    This stands for the current object
    C for a const
    X stands for a BaseFab
    B for a Box
    N for an int
*/
#define ForAllThisCBNN(T,b,ns,nc)                                       \
{                                                                       \
    BL_ASSERT(contains(b));                                             \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    const int *_th_plo = loVect();                                      \
    const int *_th_plen = length();                                     \
    const int *_b_lo = (b).loVect();                                    \
    IntVect b_length = (b).size();                                      \
    const int *_b_len = b_length.getVect();                             \
    const T* _th_p = dptr;                                              \
    for(int _n = (ns); _n < (ns)+(nc); ++_n) {                          \
        const T *_th_pp = _th_p                                         \
            + ((_b_lo[0] - _th_plo[0])                                  \
               + _n * _th_plen[0]);                                     \
        for(int _i = 0; _i < _b_len[0]; ++_i) {                         \
            const int iR = _i + _b_lo[0];                               \
            const T &thisR = _th_pp[_i];


/*
  The macro EndForTX must be used to end all ForAllThisBNNXC,
  ForAllThisBNNXCBN and ForAllThisBNNXCBNYCBN looping constructs.
*/
#define EndForTX }}}}

/*
  The macro EndFor must be used to end all ForAllXBNN, ForAllXCBNN,
  ForAllThisBNN, and ForAllThisCBNN looping constructs.
*/
#define EndFor }}}

/*
  The macro EndForPencil must be used to end ForAll*Pencil looping constructs.
*/
#define EndForPencil }}

#elif (BL_SPACEDIM == 2)



#define ForAllThisCBNN(T,b,ns,nc)                                       \
{                                                                       \
    BL_ASSERT(contains(b));                                             \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    const int *_th_plo = loVect();                                      \
    const int *_th_plen = length();                                     \
    const int *_b_lo = (b).loVect();                                    \
    IntVect b_length = (b).size();                                      \
    const int *_b_len = b_length.getVect();                             \
    const T* _th_p = dptr;                                              \
    for(int _n = (ns); _n < (ns)+(nc); ++_n) {                          \
        for(int _j = 0; _j < _b_len[1]; ++_j) {                         \
            const int jR = _j + _b_lo[1];                               \
            const T *_th_pp = _th_p                                     \
                + ((_b_lo[0] - _th_plo[0])                              \
                   + _th_plen[0]*(                                      \
                       (_j + _b_lo[1] - _th_plo[1])                     \
                       + _n * _th_plen[1]));                            \
            for(int _i = 0; _i < _b_len[0]; ++_i) {                     \
                const int iR = _i + _b_lo[0];                           \
                const T &thisR = _th_pp[_i];



#define EndFor }}}}
#define EndForTX }}}}}
#define EndForPencil }}}

#elif (BL_SPACEDIM == 3)



#define ForAllThisCBNN(T,b,ns,nc)                                       \
{                                                                       \
    BL_ASSERT(contains(b));                                             \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    const int *_th_plo = loVect();                                      \
    const int *_th_plen = length();                                     \
    const int *_b_lo = (b).loVect();                                    \
    IntVect b_length = (b).size();                                      \
    const int *_b_len = b_length.getVect();                             \
    const T* _th_p = dptr;                                              \
    for(int _n = (ns); _n < (ns)+(nc); ++_n) {                          \
        for(int _k = 0; _k < _b_len[2]; ++_k) {                         \
            const int kR = _k + _b_lo[2];                               \
            for(int _j = 0; _j < _b_len[1]; ++_j) {                     \
                const int jR = _j + _b_lo[1];                           \
                const T *_th_pp = _th_p                                 \
                    + ((_b_lo[0] - _th_plo[0])                          \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _th_plen[1]*(                              \
                               (kR - _th_plo[2])                        \
                               + _n * _th_plen[2])));                   \
                for(int _i = 0; _i < _b_len[0]; ++_i) {                 \
                    const int iR = _i + _b_lo[0];                       \
                    const T &thisR = _th_pp[_i];




#define EndFor }}}}}
#define EndForTX }}}}}}
#define EndForPencil }}}}

#endif


#endif /*BL_LOOPING_H*/
