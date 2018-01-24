
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
  The macro ForAllThisBNNXC(T,b,ns,nc,x,nss) is intended to facilitate
  efficient looping over the contents of BaseFabs and objects derived from
  BaseFab.  Special attention has been paid to make it work efficiently on
  vector supercomputers.

  This macro acts upon the BaseFab *this and in addition is able to utiliize
  values in the const BaseFab x.  The loop runs over the points in the Box b
  and over components starting at ns and ending at ns+nc-1.  The reference
  variables are thisR and xR, respectively.  As usual the x in xR is replaced
  by the macro's fifth argument.  The sixth argument nss is the number of the
  argument in x that corresponds to the ns argument in *this.

  Looping macro mnemonics:

    This stands for the current object
    C for a const
    X stands for a BaseFab
    B for a Box
    N for an int
*/
#define ForAllThisBNNXC(T,b,ns,nc,x,nss)                                \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= box();                                                  \
    _subbox_ &= b;                                                      \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n; nR += 0;                                       \
            T *_th_pp = _th_p                                           \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_subbox_lo[0] - _x_plo[0])                          \
                   + _n * _x_plen[0]);                                  \
            for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {      \
                int iR = _i + _subbox_lo[0]; iR += 0;                   \
                T &thisR = * _th_pp; const T & x##R = _x_pp[_i];


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

#define ForAllThisBNNXC(T,b,ns,nc,x,nss)                                \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= box();                                                  \
    _subbox_ &= b;                                                      \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n; nR += 0;                                       \
            for(int _j = 0; _j < _subbox_len[1]; ++_j) {                \
                const int jR = _j + _subbox_lo[1];                      \
                T *_th_pp = _th_p                                       \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_subbox_lo[0] - _x_plo[0])                      \
                       + _x_plen[0]*(                                   \
                           (jR - _x_plo[1])                             \
                           + _n * _x_plen[1]));                         \
                for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {  \
                    int iR = _i + _subbox_lo[0]; iR += 0;               \
                    T &thisR = * _th_pp; const T & x##R = _x_pp[_i];



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

#define ForAllThisBNNXC(T,b,ns,nc,x,nss)                                \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_((x).box());                                            \
    _subbox_ &= box();                                                  \
    _subbox_ &= b;                                                      \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n; nR += 0;                                       \
            for(int _k = 0; _k < _subbox_len[2]; ++_k) {                \
                const int kR = _k + _subbox_lo[2];                      \
                for(int _j = 0; _j < _subbox_len[1]; ++_j) {            \
                    const int jR = _j + _subbox_lo[1];                  \
                    T *_th_pp = _th_p                                   \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_subbox_lo[0] - _x_plo[0])                  \
                           + _x_plen[0]*(                               \
                               (jR - _x_plo[1])                         \
                               + _x_plen[1]*(                           \
                                   (kR - _x_plo[2])                     \
                                   + _n * _x_plen[2])));                \
                    for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) { \
                        int iR = _i + _subbox_lo[0]; iR += 0;           \
                        T &thisR = * _th_pp; const T & x##R = _x_pp[_i];



#define EndFor }}}}}
#define EndForTX }}}}}}
#define EndForPencil }}}}

#endif

#define ForAllThisXC(T,x)       ForAllThisBNNXC(T,domain,0,nComp(),x,0)

#endif /*BL_LOOPING_H*/
