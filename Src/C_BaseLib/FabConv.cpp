//BL_COPYRIGHT_NOTICE

//
// $Id: FabConv.cpp,v 1.2 1997-09-18 20:12:48 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <cstdlib>
#include <climits>
#include <cstring>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#endif

#include <BoxLib.H>
#include <FabConv.H>
#include <FArrayBox.H>
#include <FPC.H>
#include <Misc.H>
#include <REAL.H>

//
// Declarations of Cray-specific FP format to IEEE routines.
//

#if defined(BL_ARCH_CRAY)
#define FORT_IEG2CRAY  IEG2CRAY
#define FORT_CRAY2IEG  CRAY2IEG
extern "C"
{
    void FORT_IEG2CRAY(int& type, int& num, char* forn, int& bitoff,
                       char* cry, int& stride, char& craych);
    void FORT_CRAY2IEG(int& type, int& num, char* forn, int& bitoff,
                       char* cry, int& stride, char& craych);
}
#endif /*defined(BL_ARCH_CRAY)*/

//
// This is not inlined as it's an inherited virtual.
//

RealDescriptor*
RealDescriptor::clone () const
{
    RealDescriptor* rd = new RealDescriptor(*this);
    if (rd == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
    return rd;
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
            BoxLib::Error("selectOrdering(): Crazy ordering");
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
            BoxLib::Error("selectOrdering(): Crazy ordering");
        }
        break;
    default:
        BoxLib::Error("selectOrdering(): Crazy precision");
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
            if (rd == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            return rd;
        case FABio::FAB_DOUBLE:
            rd = new RealDescriptor(FPC::ieee_double, ord, 8);
            if (rd == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            return rd;
        }
    }
    case FABio::FAB_NATIVE:
        if (sys != 0 && strncmp(sys, "CRAY", 4) == 0)
        {
            rd = new RealDescriptor(FPC::cray_float, FPC::cray_float_order, 8);
            if (rd == 0)
                BoxLib::OutOfMemory(__FILE__, __LINE__);
            return rd;
        }
    default:
        BoxLib::Error("RealDescriptor::newRealDescriptor(): Crazy precision");
    }
    rd = new RealDescriptor;
    if (rd == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);
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
_PD_get_bit (char*      base,
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
_PD_extract_field (char*      in,
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
_PD_btrvout (char* out,
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
_PD_insert_field (long  in_long,
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
	_PD_btrvout(in, l_bytes, 1L);
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
_PD_set_bit (char* base, int offs)
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
_PD_reorder (char*      arr,
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
// It maps the `in' array into the `out' array, changing the ordering
// from inord to outord.
//

static
void
permute_real_word_order (void*       out,
                         const void* in,
                         long        nitems,
                         const int*  outord,
                         const int*  inord)
{
    const int REALSIZE = sizeof(Real);
    
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

// _PD_FCONVERT - general floating point conversion routine
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
    assert(number == nitems);
    memset(out, 0, number*outbytes);

    lout = (char*)out;
    lin  = (char*)in;

    for (i = 0L; i < nitems; i++)
    {
        //
        // Move the exponent over.
        //
	expn = _PD_extract_field(lin, bi_exp, nbi_exp, inbytes, inord);
	sign = _PD_get_bit(lin, bi_sign, inbytes, inord);
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
	    _PD_insert_field(expn, nbo_exp, lout, bo_exp, l_order, l_bytes);

	    if (sign)
		_PD_set_bit(lout, bo_sign);

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
		_PD_set_bit(lout, indxout);
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
		mant  = _PD_extract_field(lin, indxin, nbits, inbytes, inord);
                //
                // Do complement for negative ones complement data.
                //
		if (onescmp && sign)
		    ONES_COMP_NEG(mant, nbits, 0L);

		_PD_insert_field(mant, nbits, lout, indxout, l_order, l_bytes);

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
	    _PD_insert_field(expn_max, nbo_exp, lout, bo_exp, l_order, l_bytes);

	    if (_PD_get_bit(lin, bi_sign, inbytes, inord))
		_PD_set_bit(lout, bo_sign);
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
    _PD_reorder((char*)out, nitems, outbytes, outord);
}

//
// These routines are specialized to do important floating point
// conversions on the alpha.  They make use of architectural knowledege to
// run faster.  Please add more of these as needed.
//

#if defined(__alpha) && !defined(BL_USE_FLOAT)
static
void
cray64toalpha64_fconvert (void*       out,
                          const void* in,
                          long        nitems)
{
    const long DeltaBias = 0x3FFL - 0x4000L - 1L;
    const long expn_max  = (1L << 11L) - 1L;

    memset(out, 0, nitems*8);

    long *lin  = (long *)in;
    long *lout = (long *)out;

    for (long i = 0; i < nitems; i++)
    {
        //
        // Step 1: change ordering to match alpha.
        //
      long ordinp = 0, input = *lin;
      for (size_t j = 0; j < sizeof(long); j++)
      {
          ordinp <<= 8;
          ordinp |= (input & 0xff);
          input >>= 8;
      }
      //
      // Step 2: extract sign, exponent as longs.
      //
      long sign   = (ordinp>>63) & 1;
      long expn   = (ordinp>>48) & 0x7FFF;
      long ordout = 0;
      //
      // Step 3: add biases.
      //
      if (expn != 0)
          expn += DeltaBias;
      if (0 <= expn && expn < expn_max)
      {
          ordout |= (sign<<63);
          ordout |= (expn<<52);
          //
          // Step 4 get the mantissa, keeping in mind cray HSB convention.
          //
          long mant = (ordinp) & 0x7FFFFFFFFFFF;
          mant <<= 5;
          ordout |= mant;
      }
      else if (expn_max <= expn)
      {
          //
          // Overflow.  Make something big.
          //
          ordout = 0x7ff0000000000000L;
          ordout |= (sign<<63);
      }
      else
          //
          // Denorm?
          //
          ordout = 0;
      //
      // Step last: store results and update pointers.
      //
      lin++;
      *lout++ = ordout;
    }
}
#endif /*defined(__alpha) && !defined(BL_USE_FLOAT)*/

#if defined(__alpha) && defined(BL_USE_FLOAT)
static
void
ieee32toalpha32_fconvert (void*       out,
                          const void* in,
                          long        nitems)
{
    const long expn_max = (1L << 8L) - 1L;

    memset(out, 0, nitems*4);

    int* iin  = (int*)in;
    int* iout = (int*)out;

    for (long i = 0; i < nitems; i++)
    {
        //
        // Step 1: change ordering to match alpha.
        //
        int ordinp = 0;
        for (int j = 0, input = *iin; j < sizeof(int); j++)
        {
            ordinp <<= 8;
            ordinp |= (input & 0xff);
            input >>= 8;
        }
        //
        // Step 2: extract exponent.
        //
        long expn = (ordinp>>23) & 0xFF;
        if (expn_max <= expn)
        {
            //
            // Overflow.  Make something big.
            //
            int sign = ordinp & 0x80000000;
            ordinp = 0x7f800000;
            ordinp |= sign;
        }
        else if (expn <= 0)
            //
            // Denorm?
            //
            ordinp = 0;
        //
        // Step last: store results and update pointers.
        //
        iin++;
        *iout++ = ordinp;
    }
}
#endif /*defined(__alpha) && defined(BL_USE_FLOAT)*/

static
void
PD_fixdenormals (void*       out,
                 long        nitems,
                 const long* outfor,
                 const int*  outord)
{
    const int nbo = int(outfor[0]);

    int nbo_exp  = int(outfor[1]);
    int bo_exp   = int(outfor[4]);
    int outbytes = (nbo + 7) >> 3;

    char* lout = (char*) out;

    for (long i = 0L; i < nitems; i++)
    {
        if (_PD_extract_field(lout, bo_exp, nbo_exp, outbytes, outord) == 0)
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
getarray (istream&       is,                                       \
          Array< TYPE >& ar)                                       \
{                                                                  \
    char c;                                                        \
    is >> c;                                                       \
    if (c != '(')                                                  \
        BoxLib::Error("getarray(istream&): expected a \'(\'");     \
    int size;                                                      \
    is >> size;                                                    \
    is >> c;                                                       \
    if ( c != ',')                                                 \
        BoxLib::Error("getarray(istream&): expected a \',\'");     \
    is >> c;                                                       \
    if (c != '(')                                                  \
        BoxLib::Error("getarray(istream&): expected a \'(\'");     \
    ar.resize(size);                                               \
    for(int i = 0; i < size; ++i)                                  \
        is >> ar[i];                                               \
    is >> c;                                                       \
    if (c != ')')                                                  \
        BoxLib::Error("getarray(istream&): expected a \')\'");     \
    is >> c;                                                       \
    if (c != ')')                                                  \
        BoxLib::Error("getarray(istream&): expected a \')\'");     \
}
GETARRAY(int)
GETARRAY(long)
#undef GETARRAY

#undef  PUTARRAY
#define PUTARRAY(TYPE)                 \
static                                 \
void                                   \
putarray (ostream&             os,     \
          const Array< TYPE >& ar)     \
{                                      \
    int i;                             \
    os << '(';                         \
    os << ar.length() << ", (";        \
    for (i = 0; i < ar.length(); ++i)  \
    {                                  \
        os << ar[i];                   \
        if (i != ar.length() - 1)      \
            os << " ";                 \
    }                                  \
    os << "))";                        \
}
PUTARRAY(int)
PUTARRAY(long)
#undef PUTARRAY

ostream&
operator<< (ostream&              os,
            const RealDescriptor& id)
{
    os << "(";
    putarray(os, id.fr);
    os << ',';
    putarray(os, id.ord);
    os << ")";
    if (os.fail())
        BoxLib::Error("operator<<(ostream&,RealDescriptor&) failed");
    return os;
}

istream&
operator>> (istream&        is,
            RealDescriptor& rd)
{
    char c;
    is >> c;
    if (c != '(')
        BoxLib::Error("operator>>(istream&,RealDescriptor&): expected a \'(\'");
    getarray(is, rd.fr);
    is >> c;
    if (c != ',')
        BoxLib::Error("operator>>(istream&,RealDescriptor&): expected a \',\'");
    getarray(is, rd.ord);
    is >> c;
    if (c != ')')
        BoxLib::Error("operator>>(istream&,RealDescriptor&): expected a \')\'");
    return is;
}

static
void
PD_convert (void*                 out,
            const void*           in,
            long                  nitems,
            int                   boffs,
            const RealDescriptor& od,
            const RealDescriptor& id,
            const IntDescriptor&  ld,
            int                   onescmp = 0)
{
#if defined(__alpha) && !defined(BL_USE_FLOAT)
    if (id == FPC::CrayRealDescriptor() && od == FPC::NativeRealDescriptor())
        cray64toalpha64_fconvert(out, in, nitems);
    else
#endif /*defined(__alpha) && !defined(BL_USE_FLOAT)*/

#if defined(__alpha) & defined( BL_USE_FLOAT)
  if (id==FPC::Ieee32NormalRealDescriptor() && od==FPC::NativeRealDescriptor())
      ieee32toalpha32_fconvert(out, in, nitems);
  else
#endif

#if defined(BL_ARCH_CRAY)
    if (id==FPC::NativeRealDescriptor() && od==FPC::Ieee32NormalRealDescriptor())
    {
        char craych;
        int wdsize = 4, conv_typ = 2, stride = 1, offset = 0, len = nitems;
        assert(len == nitems);
        FORT_CRAY2IEG(conv_typ,len,(char*)out,offset,(char*)in,stride,craych);
    }
    else
    if (id==FPC::NativeRealDescriptor() && od==FPC::Ieee64NormalRealDescriptor())
    {
        //
        // Note that there is currently no way to specify that we want to
        // write out IEEE64 on a Cray.  In fact, I don't believe there is
        // any need for this.  Hence this block can never be reached.  I'm
        // leaving it in solely for symmetry.  Of course, there may be some
        // use for it in the future.
        //
        char craych;
        int wdsize = 8, conv_typ = 8, stride = 1, offset = 0, len = nitems;
        assert(len == nitems);
        FORT_CRAY2IEG(conv_typ,len,(char*)out,offset,(char*)in,stride,craych);
    }
    else
    if (id==FPC::Ieee32NormalRealDescriptor() && od==FPC::NativeRealDescriptor())
    {
        char craych;
        int wdsize = 4, conv_typ = 2, stride = 1, offset = 0, len = nitems;
        assert(len == nitems);
        FORT_IEG2CRAY(conv_typ,len,(char*)in,offset,(char*)out,stride,craych);
    }
    else
    if (id==FPC::Ieee64NormalRealDescriptor() && od==FPC::NativeRealDescriptor())
    {
        char craych;
        int wdsize = 8, conv_typ = 8, stride = 1, offset = 0, len = nitems;
        assert(len == nitems);
        FORT_IEG2CRAY(conv_typ,len,(char*)in,offset,(char*)out,stride,craych);
    }
    else
#endif /*defined(BL_ARCH_CRAY)*/

    if (od == id && boffs == 0)
    {
        size_t n = size_t(nitems);
        assert(n == nitems);
        memcpy(out, in, n*od.numBytes());
    }
    else if (od.formatarray() == id.formatarray() && boffs == 0 && !onescmp)
        permute_real_word_order(out, in, nitems, od.order(), id.order());
    else
    {
        PD_fconvert(out, in, nitems, boffs, od.format(), od.order(),
                    id.format(), id.order(), ld.order(), ld.numBytes(),
                    onescmp);
        PD_fixdenormals(out, nitems, od.format(), od.order());
    }
}

//
// These routines aren't currently used.  Eventually, we may define
// convert() functions for integral types at which time these'll be used.
//

#if 0
//
// Convert ones complement integers to twos complement.
// Modern machines use twos complement arithmetic and therefore
// this is a one way conversion.
//

static
void
_PD_ones_complement (char* out,
                     long  nitems,
                     int   nbo)
{
    //
    // The next two lines are to get around a KCC warning message.
    // We've got to use signed chars here, but KCC won't let us
    // simply cast out to lout directly.
    //
    void*        vout = out;
    signed char* lout = (signed char *) vout;

    for (int i = 0L; i < nitems; i++)
    {
        if (*lout < 0)
        {
            unsigned int carry = 1;
            for (int j = nbo-1; (j >= 0) && (carry > 0); j--)
            {
                carry  += lout[j];
                lout[j] = carry & 0xFF;
                carry   = (carry > 0xFF);
            }
        }
        lout += nbo;
    }
}

//
// Convert integers of nbi bytes to integers of nbo bytes.
// The number of bytes for each integer are given.
//

static
void
PD_iconvert (void* out,
             void* in,
             long  nitems,
             long  nbo,
             int   ordo,
             long  nbi,
             int   ordi,
             int   onescmp = 0)
{
    long i;
    int j;
    char *lout, *lin, *po, *pi;

    lin  = (char*) in;
    lout = (char*) out;
    //
    // Convert nitems integers test sign bit to properly convert
    // negative integers.
    //
    if (nbi < nbo)
    {
        if (ordi == REVERSE_ORDER)
        {
            for (j = nbi; j < nbo; j++)
            {
                po = lout + j - nbi;
                pi = lin + nbi - 1;
                for (i = 0L; i < nitems; i++)
                {
                    *po = (*pi & 0x80) ? 0xff : 0;
                    po += nbo;
                    pi += nbi;
                }
            }
            for (j = nbi; j > 0; j--)
            {
                po = lout + nbo - j;
                pi = lin + j - 1;
                for (i = 0L; i < nitems; i++)
                {
                    *po = *pi;
                    po += nbo;
                    pi += nbi;
                }
            }
        }
        else
        {
            for (j = nbi; j < nbo; j++)
            {
                po = lout + j - nbi;
                pi = lin;
                for (i = 0L; i < nitems; i++)
                {
                    *po = (*pi & 0x80) ? 0xff : 0;
                    po += nbo;
                    pi += nbi;
                }
            }
            for (j = 0; j < nbi; j++)
            {
                po = lout + j + nbo - nbi;
                pi = lin + j;
                for (i = 0L; i < nitems; i++)
                {
                    *po = *pi;
                    po += nbo;
                    pi += nbi;
                }
            }
        }
    }
    else if (nbi >= nbo)
    {
        if (ordi == REVERSE_ORDER)
        {
            for (j = nbo; j > 0; j--)
            {
                po = lout + nbo - j;
                pi = lin + j - 1;
                for (i = 0L; i < nitems; i++)
                {
                    *po = *pi;
                    po += nbo;
                    pi += nbi;
                }
            }
        }
        else
        {
            for (j = nbi - nbo; j < nbi; j++)
            {
                po = lout + j - nbi + nbo;
                pi = lin + j;
                for (i = 0L; i < nitems; i++)
                {
                    *po = *pi;
                    po += nbo;
                    pi += nbi;
                }
            }
        }
    }
    //
    // If input used ones complement arithmetic convert to twos complement.
    //
    if (onescmp)
        _PD_ones_complement((char*)out, nitems, nbo);

    if (ordo == REVERSE_ORDER)
        _PD_btrvout((char*)out, nbo, nitems);
}

static
void
PD_convert (void*                out,
            void*                in,
            long                 nitems,
            const IntDescriptor& od,
            const IntDescriptor& id,
            int                  onescmp = 0)
{
    if (od == id)
        memcpy(out, in, nitems*od.numBytes());
    else
    {
        PD_iconvert(out,
                    in,
                    nitems,
                    od.numBytes(),
                    od.order(),
                    id.numBytes(),
                    id.order(),
                    onescmp);
    }
}
#endif /*#if 0*/

//
// Convert nitems in RealDescriptor format to native Real format.
//

void
RealDescriptor::convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       void*                 in,
                                       const RealDescriptor& id)
{
    PD_convert(out,
               in,
               nitems,
               0,
               FPC::NativeRealDescriptor(),
               id,
               FPC::NativeLongDescriptor());
}

//
// Read nitems from istream in ReadDescriptor format to native Real format.
//

void
RealDescriptor::convertToNativeFormat (Real*                 out,
                                       long                  nitems,
                                       istream&              is,
                                       const RealDescriptor& id)
{
    const int SHOULDREAD = 8192;

    char* bufr = new char[SHOULDREAD * id.numBytes()];

    if (bufr == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    while (nitems > 0)
    {
        int get = int(nitems) > SHOULDREAD ? SHOULDREAD : int(nitems);
        is.read(bufr, id.numBytes()*get);
        PD_convert(out,
                   bufr,
                   get,
                   0,
                   FPC::NativeRealDescriptor(),
                   id,
                   FPC::NativeLongDescriptor());
        nitems -= get;
        out    += get;
    }

    if (is.fail())
        BoxLib::Error("convert(Real*,long,istream&,RealDescriptor&) failed");

    delete [] bufr;
}

//
// Convert nitems Reals in native format to RealDescriptor format.
//

void
RealDescriptor::convertFromNativeFormat (void*                 out,
                                         long                  nitems,
                                         Real*                 in,
                                         const RealDescriptor& od)
{
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
RealDescriptor::convertFromNativeFormat (ostream&              os,
                                         long                  nitems,
                                         const Real*           in,
                                         const RealDescriptor& od)
{
    const int SHOULDWRITE = 8192;

    char* bufr = new char[SHOULDWRITE * od.numBytes()];

    if (bufr == 0)
        BoxLib::OutOfMemory(__FILE__, __LINE__);

    while (nitems > 0)
    {
        int put = int(nitems) > SHOULDWRITE ? SHOULDWRITE : int(nitems);
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

    if (os.fail())
        BoxLib::Error("convert(ostream&,long,Real*,RealDescriptor&): failed");

    delete [] bufr;
}
