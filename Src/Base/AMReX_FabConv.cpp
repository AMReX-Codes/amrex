
#include <iostream>
#include <cstdlib>
#include <limits>
#include <cstring>

#include <AMReX.H>
#include <AMReX_FabConv.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FPC.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

namespace amrex {

bool RealDescriptor::bAlwaysFixDenormals (false);
int  RealDescriptor::writeBufferSize(262144);  // ---- these are number of reals,
int  RealDescriptor::readBufferSize(262144);   // ---- not bytes

IntDescriptor::IntDescriptor () {}

IntDescriptor::IntDescriptor (long     nb,
                              Ordering o)
    : numbytes(nb),
      ord(o)
{}

IntDescriptor::Ordering
IntDescriptor::order () const
{
    return ord;
}

int
IntDescriptor::numBytes () const
{
    return numbytes;
}

bool
IntDescriptor::operator== (const IntDescriptor& id) const
{
    return ord == id.ord && numbytes == id.numbytes;
}

bool
IntDescriptor::operator!= (const IntDescriptor& id) const
{
    return !operator==(id);
}

std::ostream&
operator<< (std::ostream& os,
            const IntDescriptor& id)
{
    amrex::StreamRetry sr(os, "opRD", 4);

    while(sr.TryOutput()) {
        os << "(";
        os << id.numBytes();
        os << ',';
        os << id.order();
        os << ")";
    }
    return os;
}

std::istream&
operator>> (std::istream& is,
            IntDescriptor& id)
{
    char c;
    is >> c;
    if (c != '(')
        amrex::Error("operator>>(istream&,RealDescriptor&): expected a \'(\'");
    int numbytes;
    is >> numbytes;
    id.numbytes = numbytes;
    is >> c;
    if (c != ',')
        amrex::Error("operator>>(istream&,RealDescriptor&): expected a \',\'");
    int ord;
    is >> ord;
    id.ord = (IntDescriptor::Ordering) ord;
    is >> c;
    if (c != ')')
        amrex::Error("operator>>(istream&,RealDescriptor&): expected a \')\'");
    return is;
}

RealDescriptor::RealDescriptor ()
{}

RealDescriptor::RealDescriptor (const long* fr_,
                                const int*  ord_,
                                int         ordl_)
    : fr(fr_, fr_+8),
      ord(ord_, ord_+ordl_)
{}

const long*
RealDescriptor::format () const &
{
    BL_ASSERT(fr.size() != 0);
    return fr.dataPtr();
}

const Vector<long>&
RealDescriptor::formatarray () const &
{
    BL_ASSERT(fr.size() != 0);
    return fr;
}

const int*
RealDescriptor::order () const &
{
    BL_ASSERT(ord.size() != 0);
    return ord.dataPtr();
}

const Vector<int>&
RealDescriptor::orderarray () const &
{
    BL_ASSERT(ord.size() != 0);
    return ord;
}

int
RealDescriptor::numBytes () const
{
    BL_ASSERT(fr.size() != 0);
    return (fr[0] + 7 ) >> 3;
}

bool
RealDescriptor::operator== (const RealDescriptor& rd) const
{
    return fr == rd.fr && ord == rd.ord;
}

bool
RealDescriptor::operator != (const RealDescriptor& rd) const
{
    return !operator==(rd);
}

void
RealDescriptor::SetFixDenormals()
{
    bAlwaysFixDenormals = true;
}

void
RealDescriptor::SetReadBufferSize(int rbs)
{
    BL_ASSERT(rbs > 0);
    readBufferSize = rbs;
}

void
RealDescriptor::SetWriteBufferSize(int wbs)
{
    BL_ASSERT(wbs > 0);
    writeBufferSize = wbs;
}

RealDescriptor*
RealDescriptor::clone () const
{
    return new RealDescriptor(*this);
}

//
// This exists solely to support reading "old" FABs.
//

static
const int*
selectOrdering (int prec,
                int ordering)
{
    switch (prec)
    {
    case FABio::FAB_FLOAT:
        switch (ordering)
        {
        case FABio::FAB_NORMAL_ORDER:
            return FPC::normal_float_order;
        case FABio::FAB_REVERSE_ORDER:
            return FPC::reverse_float_order;
        case FABio::FAB_REVERSE_ORDER_2:
            return FPC::reverse_float_order_2;
        default:
            amrex::Error("selectOrdering(): Crazy ordering");
        }
        break;
    case FABio::FAB_DOUBLE:
        switch (ordering)
        {
        case FABio::FAB_NORMAL_ORDER:
            return FPC::normal_double_order;
        case FABio::FAB_REVERSE_ORDER:
            return FPC::reverse_double_order;
        case FABio::FAB_REVERSE_ORDER_2:
            return FPC::reverse_double_order_2;
        default:
            amrex::Error("selectOrdering(): Crazy ordering");
        }
        break;
    default:
        amrex::Error("selectOrdering(): Crazy precision");
    }
    return 0;
}

//
// This is here solely to support reading "old" FABs.
//

RealDescriptor*
RealDescriptor::newRealDescriptor (int         iot,
                                   int         prec,
                                   const char* sys,
                                   int         ordering)
{
    RealDescriptor* rd = 0;

    switch (iot)
    {
    case FABio::FAB_IEEE:
    {
        const int* ord = selectOrdering(prec, ordering);
        switch (prec)
        {
        case FABio::FAB_FLOAT:
            rd = new RealDescriptor(FPC::ieee_float, ord, 4);
            return rd;
        case FABio::FAB_DOUBLE:
            rd = new RealDescriptor(FPC::ieee_double, ord, 8);
            return rd;
        }
    }
    case FABio::FAB_NATIVE:
    default:
        amrex::Error("RealDescriptor::newRealDescriptor(): Crazy precision");
    }
    rd = new RealDescriptor;
    return rd;
}

inline
void
ONES_COMP_NEG (long& n,
               int   nb,
               long  incr)
{
    if (nb == 8*sizeof(long))
        n = ~n + incr;
    else
    {
        const long MSK = (1L << nb) - 1L;
        n = (~n + incr) & MSK;
    }
}

//
// Return bit specified as on offset from the given pointer.
//

inline
int
_pd_get_bit (char*      base,
             int        offs,
             int        nby,
             const int* ord)
{
    int n      = offs >> 3;
    int nbytes = n % nby;
    n     -= nbytes;
    offs   = offs % 8;

    if (ord == NULL)
        base += (n + nbytes);
    else
        base += (n + (ord[nbytes] - 1));

    int mask = (1 << (7 - offs));

    return (*base & mask) != 0;
}

//
// Make a copy of the bit field specified by the starting bit, OFFS and
// the number of bits, NBI, from the byte array pointed to by IN.
// All indexing is 0 based.  The copy is to be put in a long and returned.
// This imposes a 32 bit limit (minimum) so repeated calls - must be made
// for longer fields
//

static
long
_pd_extract_field (char*      in,
                   int        offs,
                   int        nbi,
                   int        nby,
                   const int* ord)
{
    int ind;
    long bit_field = 0L;
    //
    // Move past the apropriate number of bytes so that the start bit is
    // in the first byte.  OFFY is the offset of the byte containing the
    // bit OFFS
    //
    long n   = offs >> 3;
    int offy = int(n % nby);
    n   -= offy;
    offs = offs % 8;
    //
    // Advance the pointer past the unneeded items.
    //
    in += n;
    unsigned char bpb = 8 - offs;

    if (ord == NULL)
        ind = offy++;
    else
    {
        if (offy >= nby)
        {
            offy -= nby;
            in   += nby;
        }
        ind = (ord[offy++] - 1);
    }

    int tgt  = in[ind];
    unsigned char mask = (1 << bpb) - 1;
    bit_field = ((bit_field << bpb) | (tgt & mask));
    nbi -= bpb;
    if (nbi < 0)
        bit_field = bit_field >> (-nbi);
    else
    {
        for (; nbi > 0; nbi -= bpb)
        {
            //
            // ind  = (ord == NULL) ? offy++ : (ord[offy++] - 1);
            //
            if (ord == NULL)
                ind = offy++;
            else
            {
                if (offy >= nby)
                {
                    offy -= nby;
                    in += nby;
                }
                ind = (ord[offy++] - 1);
            }

            tgt  = in[ind];
            bpb  = nbi > 8 ? 8 : nbi;
            mask = (1 << bpb) - 1;
            bit_field = ((bit_field << bpb) | ((tgt >> (8 - bpb)) & mask));
        }
    }

    return bit_field;
}

//
// Byte reverse nitems words.  Each word is nb bytes long where nb is even.
//

static
void
_pd_btrvout (char* out,
             long  nb,
             long  nitems)
{
    for (long jl = 0, nbo2 = nb >> 1; jl < nbo2; jl++)
    {
        long jh  = nb - jl - 1;
        char* p1 = out + jh;
        char* p2 = out + jl;
        for (long i = 0L; i < nitems; i++)
        {
            char tmp = *p1;
            *p1 = *p2;
            *p2 = tmp;
            p1 += nb;
            p2 += nb;
        }
    }
}

const int BitsMax       = 8*sizeof(long);
const int REVERSE_ORDER = 2;

//
// Copy the least significant NB bits from the given long into the byte array
// pointed to by OUT.  All indexing is 0 based.  OFFS is the offset from the
// beginning of OUT in bits.  This assumes that the output bit array is
// initialized to all zeros after offs.
//

inline
void
_pd_insert_field (long  in_long,
                  int   nb,
                  char* out,
                  int   offs,
                  int   l_order,
                  int   l_bytes)
{
    int  dm;
    long longmask;

    char* in = (char *) &in_long;
    //
    // If the output start bit is not in the first byte move past the
    // apropriate number of bytes so that the start bit is in the first byte.
    //
    if (offs > 7)
    {
        out  += (offs >> 3);
        offs %= 8;
    }
    //
    // If mi is less than offs, copy the first dm bits over, reset offs to 0,
    // Advance mi by dm, and handle the rest as if mi >= offs.
    //
    int mi = BitsMax - nb;
    if (mi < offs)
    {
        dm = BitsMax - (8 - offs);
        if (nb == BitsMax)
            longmask = ~((1L << dm) - 1L);
        else
            longmask = ((1L << nb) - 1L) ^ ((1L << dm) - 1L);

        unsigned char fb = ((in_long&longmask)>>dm)&((1L<<(nb-dm))-1L);
        *(out++) |= fb;

        mi  += 8 - offs;
        offs = 0;
    }
    //
    // Assuming mi >= offs, left shift the input so that it is bit aligned
    // with the output.
    //
    dm       = mi - offs;
    longmask = ~((1L << dm) - 1L);
    in_long  = (in_long << dm) & longmask;
    //
    // Reorder the bytes apropriately.
    //
    if (l_order == REVERSE_ORDER)
        _pd_btrvout(in, l_bytes, 1L);
    //
    // Copy the remaining aligned bytes over.
    //
    for (int n = (offs+nb+7)/8; n > 0; n--, *(out++) |= *(in++))
        ;
}

//
// Set the bit specified as on offset from the given pointer.
//

inline
void
_pd_set_bit (char* base, int offs)
{
    int nbytes = offs >> 3;

    base += nbytes;
    offs -= 8*nbytes;

    int mask = (1 << (7 - offs));

    *base  |= mask;
}

//
// Given a pointer to an array ARR with NITEMS of NBYTES each put them
// in the order defined by ORD.  This assumes they're in the order 1 .. n
// on input.
//

static
void
_pd_reorder (char*      arr,
             long       nitems,
             int        nbytes,
             const int* ord)
{
    const int MAXLINE = 16;
    char local[MAXLINE];

    for (int j; nitems > 0; nitems--)
    {
        arr--;
        for (j = 0; j < nbytes; local[j] = arr[ord[j]], j++);
        arr++;
        for (j = 0; j < nbytes; *(arr++) = local[j++]);
    }
}

//
// This should only be called with two arrays of Reals.
// It maps the in array into the out array, changing the ordering
// from inord to outord.
//

static
void
permute_real_word_order (void*       out,
                         const void* in,
                         long        nitems,
                         const int*  outord,
                         const int*  inord, 
                         int         REALSIZE)
{
    BL_PROFILE("permute_real_word_order");
    
    char* pin  = (char*) in;
    char* pout = (char*) out;

    pin--; pout--;

    for (; nitems > 0; nitems--, pin += REALSIZE, pout += REALSIZE)
    {
        for (int i = 0; i < REALSIZE; i++)
            pout[outord[i]] = pin[inord[i]];
    }
}

//
// Parametrized Data Conversion Method
//
// Floating point formats are characterized by a set of parameters which
// describe the fundamental elements of a floating point number. These are
//
//  Sign     - always assumed to be a single bit
//           - requires bit offset
//  Exponent - assumed to be a biased integer smaller than 32 bits
//           - (this allows the conversion to use a long on all known
//           - platforms - an exponent greater than 32 bits long would
//           - allow much larger numbers than should be needed for
//           - scientific computations)
//           - requires a bit offset, a bit length, and a bias
// Mantissa  - assumed to be a bitstream of arbitrary length
//           - requires a bit offset and a bit length
// HMB       - in all floating point representations the mantissa is
//           - normalized so that the most significant bit is one.
//           - in some formats the one is explicitly included in the
//           - representation and in others it is only implicit
//           - this gives some formats an extra bit of precision.
//           - requires a flag which is TRUE if the HMB is explicit
//
// Two other factors involved are: the byte order which could be
// mixed with the bit layout of the numbers but isn't in actual practice
// on current machines; and whether one's complement or two's complement
// arithmetic is used. Modern machines all use two's complement arithmetic
// and the model used here and now is that data from one's complement
// machines is to be read only.  This restriction is relatively easy
// to relax, but there is no evidence that it should be.
//
// An issue which is not a problem in the current implementation is that
// old machines with byte sizes other than 8 bits can be accomodated
// because the conversions treat the input and output as bitstreams
// instead of bytestreams.
//
// The conversion process is summarized as follows:
//   1) Extract the sign bit and exponent field from the input number
//   2) Subtract the bias of the source format and add the bias
//      of the target format
//   3) Check for overflow in the exponent
//   4) Insert the new exponent and the sign bit in the target stream
//   5) Copy the mantissa bits from the source to the target
//      compensating for differences in the HMB between the two
//      formats
//   6) Take care of any known anomalies - e.g. CRAY format is
//      inconsistent in that the HMB is explicitly on for all numbers
//      with the exception of 0.0
//   7) Reorder the bytes of the target stream appropriately
//
// The floating point formats for a variety of platforms are supplied by
// PDBLib and are defined at the top of this file

// _pd_FCONVERT - general floating point conversion routine
//              - convert from floating point format specified by infor
//              - to format specified by outfor
//              -
//              - floating point format specification:
//              -
//              -   format[0] = # of bits per number
//              -   format[1] = # of bits in exponent
//              -   format[2] = # of bits in mantissa
//              -   format[3] = start bit of sign
//              -   format[4] = start bit of exponent
//              -   format[5] = start bit of mantissa
//              -   format[6] = high order mantissa bit (CRAY needs this)
//              -   format[7] = bias of exponent
//

void
PD_fconvert (void*       out,
             const void* in,
             long        nitems,
             int         boffs,
             const long* outfor,
             const int*  outord,
             const long* infor,
             const int*  inord,
             int         l_order,
             int         l_bytes,
             int         onescmp)
{
    BL_PROFILE("PD_fconvert");
    long i, expn, expn_max, hexpn, mant, DeltaBias, hmbo, hmbi;
    int nbits, inbytes, outbytes, sign;
    int indxin, indxout, inrem, outrem, dindx;
    int bi_sign, bo_sign, bi_exp, bo_exp, bi_mant, bo_mant;
    int nbi_exp, nbo_exp, nbi, nbo;
    char *lout, *lin;
    unsigned char *rout;

    nbi     = int(infor[0]);
    nbo     = int(outfor[0]);
    nbi_exp = int(infor[1]);
    nbo_exp = int(outfor[1]);
    bi_sign = int(infor[3] + boffs);
    bo_sign = int(outfor[3]);
    bi_exp  = int(infor[4] + boffs);
    bo_exp  = int(outfor[4]);
    bi_mant = int(infor[5] + boffs);
    bo_mant = int(outfor[5]);

    hmbo    = (outfor[6] & 1L);
    hmbi    = (infor[6] & 1L);

    inbytes   = (nbi + 7) >> 3;
    outbytes  = (nbo + 7) >> 3;
    DeltaBias = outfor[7] + hmbo - infor[7] - hmbi;
    hexpn     = 1L << (outfor[1] - 1L);
    expn_max  = (1L << outfor[1]) - 1L;

    size_t number = size_t(nitems);
    BL_ASSERT(int(number) == nitems);
    memset(out, 0, number*outbytes);

    lout = (char*)out;
    lin  = (char*)in;

    for (i = 0L; i < nitems; i++)
    {
        //
        // Move the exponent over.
        //
        expn = _pd_extract_field(lin, bi_exp, nbi_exp, inbytes, inord);
        sign = _pd_get_bit(lin, bi_sign, inbytes, inord);
        //
        // If we have a negative number and ones complement arithmetic on the
        // input side (won't have it on the output side with modern data).
        // Take the complement of the exponent and mantissa.
        //
        if (onescmp)
        {
            if (sign)
            {
                ONES_COMP_NEG(expn, nbi_exp, 1L);
            }
            else
                expn += (expn < hexpn);
        }
        if (expn != 0)
            expn += DeltaBias;
        if ((0 <= expn) && (expn < expn_max))
        {
            _pd_insert_field(expn, nbo_exp, lout, bo_exp, l_order, l_bytes);

            if (sign)
                _pd_set_bit(lout, bo_sign);

            indxin  = bi_mant;
            inrem   = int(infor[2]);
            indxout = bo_mant;
            outrem  = int(outfor[2]);
            //
            // If input high mantissa bit (HMB) is assumed 1 and not written
            // (e.g. IEEE) but output HMB is assumed 0 (e.g. CRAY) write the
            // input starting at the output HMB+1 and set the HMB.
            //
            dindx = int(hmbo - hmbi);
            if (dindx > 0)
            {
                _pd_set_bit(lout, indxout);
                indxout += dindx;
                outrem  -= dindx;
            }
            //
            // If input HMB is assumed 0 (e.g. CRAY) but output HMB is
            // assumed 1 and not written (e.g. IEEE) take the input from
            // HMB+1 and write it to output HMB.
            //
            else if (dindx < 0)
            {
                indxin -= dindx;
                inrem  += dindx;
            }
            //
            // Move the mantissa over in sizeof(long) packets.
            //
            while ((inrem > 0) && (outrem > 0))
            {
                nbits = BitsMax > inrem ? inrem : BitsMax;
                nbits = nbits > outrem ? outrem : nbits;
                mant  = _pd_extract_field(lin, indxin, nbits, inbytes, inord);
                //
                // Do complement for negative ones complement data.
                //
                if (onescmp && sign)
                    ONES_COMP_NEG(mant, nbits, 0L);

                _pd_insert_field(mant, nbits, lout, indxout, l_order, l_bytes);

                indxin  += nbits;
                indxout += nbits;
                inrem   -= nbits;
                outrem  -= nbits;
            }
        }
        //
        // In case of overflow use 1.0e+(expn_max).
        //
        else if (expn_max <= expn)
        {
            _pd_insert_field(expn_max, nbo_exp, lout, bo_exp, l_order, l_bytes);

            if (_pd_get_bit(lin, bi_sign, inbytes, inord))
                _pd_set_bit(lout, bo_sign);
        }
        bi_sign += nbi;
        bi_exp  += nbi;
        bi_mant += nbi;
        bo_sign += nbo;
        bo_exp  += nbo;
        bo_mant += nbo;
    }
    //
    // Handle CRAY inconsistency which has zero as the only floating point
    // number with a 0 in the HMB.  Also problem for IEEE 96 bit float - fixed
    // by Dave Munro.
    //
    if (hmbo)
    {
        int j, mask = (1 << (7 - bo_mant % 8));

        indxout = int(outfor[5]/8);
        rout    = (unsigned char *) out;
        for (i = 0L; i < nitems; i++, rout += outbytes)
        {
            for (j = 0; j < outbytes; j++)
                if ((j == indxout) ? (rout[j] != mask) : rout[j])
                    break;
            if (j == outbytes)
                rout[indxout] = 0;
        }
    }
    //
    // Put the output bytes into the specified order.
    //
    _pd_reorder((char*)out, nitems, outbytes, outord);
}

static
void
PD_fixdenormals (void*       out,
                 long        nitems,
                 const long* outfor,
                 const int*  outord)
{
    BL_PROFILE("PD_fixdenormals");
    const int nbo = int(outfor[0]);

    int nbo_exp  = int(outfor[1]);
    int bo_exp   = int(outfor[4]);
    int outbytes = (nbo + 7) >> 3;

    char* lout = (char*) out;

    for (long i = 0L; i < nitems; i++)
    {
        if (_pd_extract_field(lout, bo_exp, nbo_exp, outbytes, outord) == 0)
        {
            //
            // Set the word to zero.
            //
            char* loutoffset = lout+(i*outbytes);
            memset(loutoffset, '\0', outbytes);
        }
        bo_exp  += nbo;
    }
}

//
// It's really sad that I need to do this ...
//

#undef  GETARRAY
#define GETARRAY(TYPE)                                             \
static                                                             \
void                                                               \
getarray (std::istream&  is,                                       \
          Vector< TYPE >& ar)                                       \
{                                                                  \
    char c;                                                        \
    is >> c;                                                       \
    if (c != '(')                                                  \
        amrex::Error("getarray(istream&): expected a \'(\'");     \
    int size;                                                      \
    is >> size;                                                    \
    is >> c;                                                       \
    if ( c != ',')                                                 \
        amrex::Error("getarray(istream&): expected a \',\'");     \
    is >> c;                                                       \
    if (c != '(')                                                  \
        amrex::Error("getarray(istream&): expected a \'(\'");     \
    ar.resize(size);                                               \
    for(int i = 0; i < size; ++i)                                  \
        is >> ar[i];                                               \
    is >> c;                                                       \
    if (c != ')')                                                  \
        amrex::Error("getarray(istream&): expected a \')\'");     \
    is >> c;                                                       \
    if (c != ')')                                                  \
        amrex::Error("getarray(istream&): expected a \')\'");     \
}
GETARRAY(int)
GETARRAY(long)
#undef GETARRAY

#undef  PUTARRAY
#define PUTARRAY(TYPE)                 \
static                                 \
void                                   \
putarray (std::ostream&        os,     \
          const Vector< TYPE >& ar)     \
{                                      \
    int i;                             \
    os << '(';                         \
    os << ar.size() << ", (";          \
    for (i = 0; i < ar.size(); ++i)    \
    {                                  \
        os << ar[i];                   \
        if (i != ar.size() - 1)        \
            os << ' ';                 \
    }                                  \
    os << "))";                        \
}
PUTARRAY(int)
PUTARRAY(long)
#undef PUTARRAY

std::ostream&
operator<< (std::ostream&         os,
            const RealDescriptor& rd)
{
  amrex::StreamRetry sr(os, "opRD", 4);

  while(sr.TryOutput()) {
    os << "(";
    putarray(os, rd.formatarray());
    os << ',';
    putarray(os, rd.orderarray());
    os << ")";
  }
  return os;
}

std::istream&
operator>> (std::istream&   is,
            RealDescriptor& rd)
{
    char c;
    is >> c;
    if (c != '(')
        amrex::Error("operator>>(istream&,RealDescriptor&): expected a \'(\'");
    Vector<long> fmt;
    getarray(is, fmt);
    is >> c;
    if (c != ',')
        amrex::Error("operator>>(istream&,RealDescriptor&): expected a \',\'");
    Vector<int> ord;
    getarray(is, ord);
    is >> c;
    if (c != ')')
        amrex::Error("operator>>(istream&,RealDescriptor&): expected a \')\'");
    rd = RealDescriptor(fmt.dataPtr(),ord.dataPtr(),ord.size());
    return is;
}

static
void
PD_convert (void*                 out,
            const void*           in,
            long                  nitems,
            int                   boffs,
            const RealDescriptor& ord,
            const RealDescriptor& ird,
            const IntDescriptor&  iid,
            int                   onescmp = 0)
{
    BL_PROFILE("PD_convert");
    if (ord == ird && boffs == 0)
    {
        size_t n = size_t(nitems);
        BL_ASSERT(int(n) == nitems);
        memcpy(out, in, n*ord.numBytes());
    }
    else if (ord.formatarray() == ird.formatarray() && boffs == 0 && ! onescmp) {
        permute_real_word_order(out, in, nitems,
                                ord.order(), ird.order(), ord.numBytes());
    }
    else if (ird == FPC::NativeRealDescriptor() && ord == FPC::Native32RealDescriptor()) {
      const Real *rIn = static_cast<const Real *>(in);
      float *rOut= static_cast<float *>(out);
      for(long i(0); i < nitems; ++i) {
        rOut[i] = rIn[i];
      }
    }
    else
    {
        PD_fconvert(out, in, nitems, boffs, ord.format(), ord.order(),
                    ird.format(), ird.order(), iid.order(), iid.numBytes(),
                    onescmp);
        PD_fixdenormals(out, nitems, ord.format(), ord.order());
    }
}

//
// Convert nitems in RealDescriptor format to native Real format.
//

void
RealDescriptor::convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       void*                 in,
                                       const RealDescriptor& id)
{
    BL_PROFILE("RD:convertToNativeFormat_vp");

    PD_convert(out,
               in,
               nitems,
               0,
               FPC::NativeRealDescriptor(),
               id,
               FPC::NativeLongDescriptor());

    if(bAlwaysFixDenormals) {
      PD_fixdenormals(out, nitems, FPC::NativeRealDescriptor().format(),
		      FPC::NativeRealDescriptor().order());
    }
}

//
// Read nitems from istream in RealDescriptor format to native Real format.
//

void
RealDescriptor::convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       std::istream&         is,
                                       const RealDescriptor& id)
{
    BL_PROFILE("RD:convertToNativeFormat_is");

    long buffSize(std::min(long(readBufferSize), nitems));
    char *bufr = new char[buffSize * id.numBytes()];

    while (nitems > 0)
    {
        int get = int(nitems) > readBufferSize ? readBufferSize : int(nitems);
        is.read(bufr, id.numBytes()*get);
        PD_convert(out,
                   bufr,
                   get,
                   0,
                   FPC::NativeRealDescriptor(),
                   id,
                   FPC::NativeLongDescriptor());

        if(bAlwaysFixDenormals) {
          PD_fixdenormals(out, get, FPC::NativeRealDescriptor().format(),
			  FPC::NativeRealDescriptor().order());
        }
        nitems -= get;
        out    += get;
    }

    if(is.fail()) {
      amrex::Error("convert(Real*,long,istream&,RealDescriptor&) failed");
    }

    delete [] bufr;
}

//
// Convert nitems Reals in native format to RealDescriptor format.
//

void
RealDescriptor::convertFromNativeFormat (void*                 out,
                                         long                  nitems,
                                         const Real*           in,
                                         const RealDescriptor& od)
{
    BL_PROFILE("RD:convertFromNativeFormat_vp");
    PD_convert(out,
               in,
               nitems,
               0,
               od,
               FPC::NativeRealDescriptor(),
               FPC::NativeLongDescriptor());
}

//
// Convert nitems Reals in native format to RealDescriptor format
// and write them to the ostream.
//

void
RealDescriptor::convertFromNativeFormat (std::ostream&         os,
                                         long                  nitems,
                                         const Real*           in,
                                         const RealDescriptor& od)
{
  BL_PROFILE("RD:convertFromNativeFormat_os");
  long nitemsSave(nitems);
  long buffSize(std::min(long(writeBufferSize), nitems));
  const Real *inSave(in);
  amrex::StreamRetry sr(os, "RD_cFNF", 4);

  while(sr.TryOutput()) {
    nitems = nitemsSave;
    in = inSave;

    char *bufr = new char[buffSize * od.numBytes()];

    while (nitems > 0)
    {
        int put = int(nitems) > writeBufferSize ? writeBufferSize : int(nitems);
        PD_convert(bufr,
                   in,
                   put,
                   0,
                   od,
                   FPC::NativeRealDescriptor(),
                   FPC::NativeLongDescriptor());
        os.write(bufr, od.numBytes()*put);
        nitems -= put;
        in     += put;
    }

    delete [] bufr;
  }
}

//
// Convert nitems floats in native format to RealDescriptor format
// and write them to the ostream.
//

void
RealDescriptor::convertFromNativeFloatFormat (std::ostream&         os,
                                              long                  nitems,
                                              const float*          in,
                                              const RealDescriptor& od)
{
  BL_PROFILE("RD:convertFromNativeFloatFormat");
  long nitemsSave(nitems);
  long buffSize(std::min(long(writeBufferSize), nitems));
  const float *inSave(in);
  amrex::StreamRetry sr(os, "RD_cFNF", 4);

  while(sr.TryOutput()) {
    nitems = nitemsSave;
    in = inSave;

    char *bufr = new char[buffSize * od.numBytes()];

    while (nitems > 0)
    {
        int put = int(nitems) > writeBufferSize ? writeBufferSize : int(nitems);
        PD_convert(bufr,
                   in,
                   put,
                   0,
                   od,
                   FPC::Native32RealDescriptor(),
                   FPC::NativeLongDescriptor());
        os.write(bufr, od.numBytes()*put);
        nitems -= put;
        in     += put;
    }

    delete [] bufr;
  }
}

//
// Convert nitems doubles in native format to RealDescriptor format
// and write them to the ostream.
//

void
RealDescriptor::convertFromNativeDoubleFormat (std::ostream&         os,
                                               long                  nitems,
                                               const double*         in,
                                               const RealDescriptor& od)
{
  BL_PROFILE("RD:convertFromNativeDoubleFormat");
  long nitemsSave(nitems);
  long buffSize(std::min(long(writeBufferSize), nitems));
  const double *inSave(in);
  amrex::StreamRetry sr(os, "RD_cFNF", 4);

  while(sr.TryOutput()) {
    nitems = nitemsSave;
    in = inSave;

    char *bufr = new char[buffSize * od.numBytes()];

    while (nitems > 0)
    {
        int put = int(nitems) > writeBufferSize ? writeBufferSize : int(nitems);
        PD_convert(bufr,
                   in,
                   put,
                   0,
                   od,
                   FPC::Native64RealDescriptor(),
                   FPC::NativeLongDescriptor());
        os.write(bufr, od.numBytes()*put);
        nitems -= put;
        in     += put;
    }

    delete [] bufr;
  }
}

//
// Read nitems from istream in RealDescriptor format to native float format.
//

void
RealDescriptor::convertToNativeFloatFormat (float*                out,
                                            long                  nitems,
                                            std::istream&         is,
                                            const RealDescriptor& id)
{
    BL_PROFILE("RD:convertToNativeFloatFormat");

    long buffSize(std::min(long(readBufferSize), nitems));
    char *bufr = new char[buffSize * id.numBytes()];

    while (nitems > 0)
    {
        int get = int(nitems) > readBufferSize ? readBufferSize : int(nitems);
        is.read(bufr, id.numBytes()*get);
        PD_convert(out,
                   bufr,
                   get,
                   0,
                   FPC::Native32RealDescriptor(),
                   id,
                   FPC::NativeLongDescriptor());

        if(bAlwaysFixDenormals) {
          PD_fixdenormals(out, get, FPC::Native32RealDescriptor().format(),
			  FPC::Native32RealDescriptor().order());
        }
        nitems -= get;
        out    += get;
    }

    if(is.fail()) {
      amrex::Error("convert(Real*,long,istream&,RealDescriptor&) failed");
    }

    delete [] bufr;
}

//
// Read nitems from istream in RealDescriptor format to native double format.
//

void
RealDescriptor::convertToNativeDoubleFormat (double*               out,
                                             long                  nitems,
                                             std::istream&         is,
                                             const RealDescriptor& id)
{
    BL_PROFILE("RD:convertToNativeDoubleFormat");

    long buffSize(std::min(long(readBufferSize), nitems));
    char *bufr = new char[buffSize * id.numBytes()];

    while (nitems > 0)
    {
        int get = int(nitems) > readBufferSize ? readBufferSize : int(nitems);
        is.read(bufr, id.numBytes()*get);
        PD_convert(out,
                   bufr,
                   get,
                   0,
                   FPC::Native64RealDescriptor(),
                   id,
                   FPC::NativeLongDescriptor());

        if(bAlwaysFixDenormals) {
          PD_fixdenormals(out, get, FPC::Native64RealDescriptor().format(),
			  FPC::Native64RealDescriptor().order());
        }
        nitems -= get;
        out    += get;
    }

    if(is.fail()) {
      amrex::Error("convert(Real*,long,istream&,RealDescriptor&) failed");
    }

    delete [] bufr;
}

}
