
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
  The macro ForAllThisBNNXCBN(T,b,ns,nc,x,bx,nss) is intended to facilitate
  efficient looping over the contents of BaseFabs and objects derived from
  BaseFab.  Special attention has been paid to make it work efficiently on
  vector supercomputers.

  This macro acts upon the BaseFab *this and in addition is able to utiliize
  values in the const BaseFab x.  The loop runs over the points in the
  Box b with components starting at ns and ending at ns+nc-1.  The reference
  variables are thisR and xR, respectively.  As usual the x in xR is replaced
  by the macro's fifth argument.  The sixth argument nss is the number of the
  argument in x that corresponds to the ns argument in *this.  Box bx must
  be the same size as this->box() intersected with b.

  Looping macro mnemonics:

    This stands for the current object
    C for a const
    X stands for a BaseFab
    B for a Box
    N for an int
*/
#define ForAllThisBNNXCBN(T,b,ns,nc,x,bx,nss)                           \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_ = box();                                               \
    _subbox_ &= b;                                                      \
    BL_ASSERT(bx.sameSize(_subbox_));                                   \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        const int *_bx_lo = (bx).loVect();                              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nss; n##x##R += 0;                       \
            T *_th_pp = _th_p                                           \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_bx_lo[0] - _x_plo[0])                              \
                   + _n * _x_plen[0]);                                  \
            for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {      \
                int iR = _i + _subbox_lo[0]; iR += 0;                   \
                int i##x##R = _i + _bx_lo[0]; i##x##R += 0;             \
                T &thisR = * _th_pp; const T & x##R = _x_pp[_i];

/*
  The macro ForAllThisBNNXCBNYCBN(T,b,ns,nc,x,bx,nsx,y,by,nsy) is intended to
  facilitate efficient looping over the contents of BaseFabs and objects
  derived from BaseFab.  Special attention has been paid to make it work
  efficiently on vector supercomputers.

  This macro acts upon the BaseFab *this and in addition is able to utiliize
  values in the const BaseFab x and const BaseFab y.  The loop runs over the
  points in the intersection of Box b with components starting at ns and
  ending at ns+nc-1.  The reference variables are thisR, xR, and yR
  respectively. As usual the x in xR is replaced by the macro's fifth argument
  and likewise for the y in yR.  The seventh argument nsx is the number of the
  argument in x that corresponds to the ns argument in *this, and the eighth
  argument nsy is the number of the argument in y that corresponds to the ns
  argument in *this.  Boxes bx and by must be the same size as this->box()
  intersected with b.

  Looping macro mnemonics:

    This stands for the current object
    C for a const
    X stands for a BaseFab
    B for a Box
    N for an int
*/
#define ForAllThisBNNXCBNYCBN(T,b,ns,nc,x,bx,nsx,y,by,nsy)              \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    BL_ASSERT((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_ = box();                                               \
    _subbox_ &= b;                                                      \
    BL_ASSERT((bx).sameSize(_subbox_));                                 \
    BL_ASSERT((by).sameSize(_subbox_));                                 \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_y_plo = (y).loVect();                               \
        const int *_y_plen = (y).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        const int *_bx_lo = (bx).loVect();                              \
        const int *_by_lo = (by).loVect();                              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nsx);                              \
        const T* _y_p  = (y).dataPtr(nsy);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR      = _n + ns;  nR      += 0;                       \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            T *_th_pp = _th_p                                           \
                + ((_subbox_lo[0] - _th_plo[0])                         \
                   + _n * _th_plen[0]);                                 \
            const T *_x_pp = _x_p                                       \
                + ((_bx_lo[0] - _x_plo[0])                              \
                   + _n * _x_plen[0]);                                  \
            const T *_y_pp = _y_p                                       \
                + ((_by_lo[0] - _y_plo[0])                              \
                   + _n * _y_plen[0]);                                  \
            for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {      \
                int iR = _i + _subbox_lo[0];  iR += 0;                  \
                int i##x##R = _i + _bx_lo[0]; i##x##R += 0;             \
                int i##y##R = _i + _by_lo[0]; i##y##R += 0;             \
                T &thisR = * _th_pp;                                    \
                const T & x##R = _x_pp[_i];                             \
                const T & y##R = _y_pp[_i];

#define ForAllRevXBNYCBNNN(T,x,bx,nsx,y,by,nsy,nc,ri)                   \
{                                                                       \
    BL_ASSERT((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    BL_ASSERT((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    BL_ASSERT((x).contains(bx));                                        \
    BL_ASSERT((y).contains(by));                                        \
    BL_ASSERT((bx).sameSize(by));                                       \
    const int *_x_plo = (x).loVect();                                   \
    const int *_x_plen = (x).length();                                  \
    const int *_y_plo = (y).loVect();                                   \
    const int *_y_plen = (y).length();                                  \
    IntVect bx_length = (bx).size();                                    \
    const int *_len = bx_length.getVect();                              \
    const int *_bx_lo = (bx).loVect();                                  \
    const int *_by_lo = (by).loVect();                                  \
    T* _x_p  = (x).dataPtr(nsx);                                        \
    const T* _y_p  = (y).dataPtr(nsy);                                  \
    for(int _n = 0; _n < (nc); ++_n) {                                  \
        int n##x##R = _n + nsx; n##x##R += 0;                           \
        int n##y##R = _n + nsy; n##y##R += 0;                           \
        int _ix = 0;                                                    \
        T *_x_pp = _x_p                                                 \
            + ((_bx_lo[0] - _x_plo[0]) + _len[0] - 1                    \
                + _n * _x_plen[0]);                                     \
        const T *_y_pp = _y_p                                           \
            + ((_by_lo[0] - _y_plo[0])                                  \
                + _n * _y_plen[0]);                                     \
        for(int _i = 0; _i < _len[0]; ++_i, --_ix) {                    \
            T & x##R = _x_pp[_ix];                                      \
            const T & y##R = _y_pp[_i];

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

#define ForAllThisBNNXCBN(T,b,ns,nc,x,bx,nss)                           \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_ = box();                                               \
    _subbox_ &= b;                                                      \
    BL_ASSERT(bx.sameSize(_subbox_));                                   \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        const int *_bx_lo = (bx).loVect();                              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nss; n##x##R += 0;                       \
            for(int _j = 0; _j < _subbox_len[1]; ++_j) {                \
                const int jR = _j + _subbox_lo[1];                      \
                const int j##x##R = _j + _bx_lo[1];                     \
                T *_th_pp = _th_p                                       \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_bx_lo[0] - _x_plo[0])                          \
                       + _x_plen[0]*(                                   \
                           (j##x##R - _x_plo[1])                        \
                           + _n * _x_plen[1]));                         \
                for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {  \
                    int iR = _i + _subbox_lo[0]; iR += 0;               \
                    int i##x##R = _i + _bx_lo[0]; i##x##R += 0; \
                    T &thisR = * _th_pp; const T & x##R = _x_pp[_i];

#define ForAllThisBNNXCBNYCBN(T,b,ns,nc,x,bx,nsx,y,by,nsy)              \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    BL_ASSERT((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_ = box();                                               \
    _subbox_ &= b;                                                      \
    BL_ASSERT((bx).sameSize(_subbox_));                                 \
    BL_ASSERT((by).sameSize(_subbox_));                                 \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_y_plo = (y).loVect();                               \
        const int *_y_plen = (y).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        const int *_bx_lo = (bx).loVect();                              \
        const int *_by_lo = (by).loVect();                              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nsx);                              \
        const T* _y_p  = (y).dataPtr(nsy);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            for(int _j = 0; _j < _subbox_len[1]; ++_j) {                \
                const int jR = _j + _subbox_lo[1];                      \
                const int j##x##R = _j + _bx_lo[1];                     \
                const int j##y##R = _j + _by_lo[1];                     \
                T *_th_pp = _th_p                                       \
                    + ((_subbox_lo[0] - _th_plo[0])                     \
                       + _th_plen[0]*(                                  \
                           (jR - _th_plo[1])                            \
                           + _n * _th_plen[1]));                        \
                const T *_x_pp = _x_p                                   \
                    + ((_bx_lo[0] - _x_plo[0])                          \
                       + _x_plen[0]*(                                   \
                           (j##x##R - _x_plo[1])                        \
                           + _n * _x_plen[1]));                         \
                const T *_y_pp = _y_p                                   \
                    + ((_by_lo[0] - _y_plo[0])                          \
                       + _y_plen[0]*(                                   \
                           (j##y##R - _y_plo[1])                        \
                           + _n * _y_plen[1]));                         \
                for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) {  \
                    int iR = _i + _subbox_lo[0];  iR += 0;              \
                    int i##x##R = _i + _bx_lo[0]; i##x##R += 0;         \
                    int i##y##R = _i + _by_lo[0]; i##y##R += 0;         \
                    T &thisR = * _th_pp;                                \
                    const T & x##R = _x_pp[_i];                         \
                    const T & y##R = _y_pp[_i];

#define ForAllRevXBNYCBNNN(T,x,bx,nsx,y,by,nsy,nc,ir)                   \
{                                                                       \
    BL_ASSERT((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    BL_ASSERT((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    BL_ASSERT((ir) >= 0 && (ir) < BL_SPACEDIM);                         \
    BL_ASSERT((x).contains(bx));                                        \
    BL_ASSERT((y).contains(by));                                        \
    BL_ASSERT((bx).sameSize(by));                                       \
    const int *_x_plo = (x).loVect();                                   \
    const int *_x_plen = (x).length();                                  \
    const int *_y_plo = (y).loVect();                                   \
    const int *_y_plen = (y).length();                                  \
    const int *_bx_lo = (bx).loVect();                                  \
    const int *_by_lo = (by).loVect();                                  \
    IntVect bx_length = (bx).size();                                    \
    const int *_len = bx_length.getVect();                              \
    T* _x_p  = (x).dataPtr(nsx);                                        \
    const T* _y_p  = (y).dataPtr(nsy);                                  \
    for(int _n = 0; _n < (nc); ++_n) {                                  \
        int n##x##R = _n + nsx; n##x##R += 0;                           \
        int n##y##R = _n + nsy; n##y##R += 0;                           \
        for(int _j = 0; _j < _len[1]; ++_j) {                           \
            const int j##x##R = _j + _bx_lo[1];                 \
            const int jrev##x##R = _len[1]-1-_j + _bx_lo[1];            \
            const int j##y##R = _j + _by_lo[1];                 \
            T *_x_pp;                                                   \
            int _ix = 0;                                                \
            int _istrd;                                                 \
            if (ir == 0) {                                              \
                _x_pp = _x_p                                            \
                    + ((_bx_lo[0] - _x_plo[0]) + _len[0] - 1            \
                       + _x_plen[0]*(                                   \
                           (j##x##R - _x_plo[1])                        \
                           + _n * _x_plen[1]));                         \
                _istrd = -1;                                            \
            } else {                                                    \
                _x_pp = _x_p                                            \
                    + ((_bx_lo[0] - _x_plo[0])                          \
                       + _x_plen[0]*(                                   \
                           (jrev##x##R - _x_plo[1])                     \
                           + _n * _x_plen[1]));                         \
                _istrd = 1;                                             \
            }                                                           \
            const T *_y_pp = _y_p                                       \
                    + ((_by_lo[0] - _y_plo[0])                          \
                       + _y_plen[0]*(                                   \
                           (j##y##R - _y_plo[1])                        \
                           + _n * _y_plen[1]));                         \
            int _x_rev = _len[0]-1; _x_rev += 0;                        \
            for(int _i = 0; _i < _len[0]; ++_i, _ix+=_istrd) {          \
                T & x##R = _x_pp[_ix];                                  \
                const T & y##R = _y_pp[_i];


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

#define ForAllThisBNNXCBN(T,b,ns,nc,x,bx,nss)                           \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nss) >= 0 && (nss) + (nc) <= (x).nComp());               \
    Box _subbox_(box());                                                \
    _subbox_ &= b;                                                      \
    BL_ASSERT((bx).sameSize(_subbox_));                                 \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        const int *_bx_lo = (bx).loVect();                              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nss);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nss; n##x##R += 0;                       \
            for(int _k = 0; _k < _subbox_len[2]; ++_k) {                \
                const int kR = _k + _subbox_lo[2];                      \
                const int k##x##R = _k + _bx_lo[2];                     \
                for(int _j = 0; _j < _subbox_len[1]; ++_j) {            \
                    const int jR = _j + _subbox_lo[1];                  \
                    const int j##x##R = _j + _bx_lo[1];         \
                    T *_th_pp = _th_p                                   \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                    for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) { \
                        int iR = _i + _subbox_lo[0]; iR += 0;           \
                        int i##x##R = _i + _bx_lo[0]; i##x##R += 0;     \
                        T &thisR = * _th_pp; const T & x##R = _x_pp[_i];

#define ForAllThisBNNXCBNYCBN(T,b,ns,nc,x,bx,nsx,y,by,nsy)              \
{                                                                       \
    BL_ASSERT((ns) >= 0 && (ns) + (nc) <= nComp());                     \
    BL_ASSERT((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    BL_ASSERT((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    Box _subbox_(box());                                                \
    _subbox_ &= b;                                                      \
    BL_ASSERT((bx).sameSize(_subbox_));                                 \
    BL_ASSERT((by).sameSize(_subbox_));                                 \
    if(_subbox_.ok()) {                                                 \
        const int *_th_plo = loVect();                                  \
        const int *_th_plen = length();                                 \
        const int *_x_plo = (x).loVect();                               \
        const int *_x_plen = (x).length();                              \
        const int *_y_plo = (y).loVect();                               \
        const int *_y_plen = (y).length();                              \
        const int *_subbox_lo = _subbox_.loVect();                      \
        IntVect _subbox_length = _subbox_.size();                       \
        const int *_subbox_len = _subbox_length.getVect();              \
        const int *_bx_lo = (bx).loVect();                              \
        const int *_by_lo = (by).loVect();                              \
        T* _th_p = dataPtr(ns);                                         \
        const T* _x_p  = (x).dataPtr(nsx);                              \
        const T* _y_p  = (y).dataPtr(nsy);                              \
        for(int _n = 0; _n < (nc); ++_n) {                              \
            int nR = _n + ns; nR += 0;                                  \
            int n##x##R = _n + nsx; n##x##R += 0;                       \
            int n##y##R = _n + nsy; n##y##R += 0;                       \
            for(int _k = 0; _k < _subbox_len[2]; ++_k) {                \
                const int kR = _k + _subbox_lo[2];                      \
                const int k##x##R = _k + _bx_lo[2];                     \
                const int k##y##R = _k + _by_lo[2];                     \
                for(int _j = 0; _j < _subbox_len[1]; ++_j) {            \
                    const int jR = _j + _subbox_lo[1];                  \
                    const int j##x##R = _j + _bx_lo[1];                 \
                    const int j##y##R = _j + _by_lo[1];                 \
                    T *_th_pp = _th_p                                   \
                        + ((_subbox_lo[0] - _th_plo[0])                 \
                           + _th_plen[0]*(                              \
                               (jR - _th_plo[1])                        \
                               + _th_plen[1]*(                          \
                                   (kR - _th_plo[2])                    \
                                   + _n * _th_plen[2])));               \
                    const T *_x_pp = _x_p                               \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                    const T *_y_pp = _y_p                               \
                        + ((_by_lo[0] - _y_plo[0])                      \
                           + _y_plen[0]*(                               \
                               (j##y##R - _y_plo[1])                    \
                               + _y_plen[1]*(                           \
                                   (k##y##R - _y_plo[2])                \
                                   + _n * _y_plen[2])));                \
                    for(int _i = 0; _i < _subbox_len[0]; ++_i, ++_th_pp) { \
                        int iR = _i + _subbox_lo[0];  iR += 0;          \
                        int i##x##R = _i + _bx_lo[0]; i##x##R += 0;     \
                        int i##y##R = _i + _by_lo[0]; i##y##R += 0;     \
                        T &thisR = * _th_pp;                            \
                        const T & x##R = _x_pp[_i];                     \
                        const T & y##R = _y_pp[_i];

#define ForAllRevXBNYCBNNN(T,x,bx,nsx,y,by,nsy,nc,ir)                   \
{                                                                       \
    BL_ASSERT((ir) >= 0 && (ir) < BL_SPACEDIM);                         \
    BL_ASSERT((nsx) >= 0 && (nsx) + (nc) <= (x).nComp());               \
    BL_ASSERT((nsy) >= 0 && (nsy) + (nc) <= (y).nComp());               \
    BL_ASSERT((x).contains(bx));                                        \
    BL_ASSERT((y).contains(by));                                        \
    BL_ASSERT((bx).sameSize(by));                                       \
    const int *_x_plo = (x).loVect();                                   \
    const int *_x_plen = (x).length();                                  \
    const int *_y_plo = (y).loVect();                                   \
    const int *_y_plen = (y).length();                                  \
    const int *_bx_lo = (bx).loVect();                                  \
    const int *_by_lo = (by).loVect();                                  \
    IntVect bx_length = (bx).size();                                    \
    const int *_len = bx_length.getVect();                              \
    T* _x_p  = (x).dataPtr(nsx);                                        \
    const T* _y_p  = (y).dataPtr(nsy);                                  \
    for(int _n = 0; _n < (nc); ++_n) {                                  \
        int n##x##R = _n + nsx; n##x##R += 0;                           \
        int n##y##R = _n + nsy; n##y##R += 0;                           \
        for(int _k = 0; _k < _len[2]; ++_k) {                           \
            const int k##x##R = _k + _bx_lo[2];                         \
            const int krev##x##R = _len[2]-1-_k + _bx_lo[2];            \
            const int k##y##R = _k + _by_lo[2];                         \
            for(int _j = 0; _j < _len[1]; ++_j) {                       \
                const int j##x##R = _j + _bx_lo[1];                     \
                const int jrev##x##R = _len[1]-1-_j + _bx_lo[1];        \
                const int j##y##R = _j + _by_lo[1];                     \
                T *_x_pp;                                               \
                int _ix = 0;                                            \
                int _istrd = 1;                                         \
                if (ir == 0) {                                          \
                    _x_pp = _x_p                                        \
                        + ((_bx_lo[0] - _x_plo[0]) + _len[0]-1          \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                    _istrd = -1;                                        \
                } else if (ir == 1) {                                   \
                    _x_pp = _x_p                                        \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (jrev##x##R - _x_plo[1])                 \
                               + _x_plen[1]*(                           \
                                   (k##x##R - _x_plo[2])                \
                                   + _n * _x_plen[2])));                \
                } else {                                                \
                    _x_pp = _x_p                                        \
                        + ((_bx_lo[0] - _x_plo[0])                      \
                           + _x_plen[0]*(                               \
                               (j##x##R - _x_plo[1])                    \
                               + _x_plen[1]*(                           \
                                   (krev##x##R - _x_plo[2])             \
                                   + _n * _x_plen[2])));                \
                }                                                       \
                const T *_y_pp = _y_p                                   \
                    + ((_by_lo[0] - _y_plo[0])                          \
                       + _y_plen[0]*(                                   \
                           (j##y##R - _y_plo[1])                        \
                           + _y_plen[1]*(                               \
                               (k##y##R - _y_plo[2])                    \
                               + _n * _y_plen[2])));                    \
                for(int _i = 0; _i < _len[0]; ++_i, _ix += _istrd) {    \
                    T & x##R = _x_pp[_ix];                              \
                    const T & y##R = _y_pp[_i];

#define EndFor }}}}}
#define EndForTX }}}}}}
#define EndForPencil }}}}

#endif

/*
  The macro ForAllThisXC(T,x) is a shortened form of ForAllThisBNNXC
  where the Box defaults to the domain of *this and the components run over
  all the components of *this.
*/
#define ForAllThisXC(T,x)       ForAllThisBNNXC(T,domain,0,nComp(),x,0)

#endif /*BL_LOOPING_H*/
