//BL_COPYRIGHT_NOTICE

#ifndef BL_FARRAYBOX_H
#define BL_FARRAYBOX_H

//
// $Id: FArrayBox.H,v 1.8 1997-12-14 01:05:52 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#else
#include <iostream.h>
#endif

#include <Pointers.H>
#include <Box.H>
#include <BaseFab.H>
#include <ArithFab.H>
#include <NormedFab.H>
#include <OrderedFab.H>
#include <REAL.H>
#include <SPACE.H>

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
#include <PointDomain.H>
#include <PointFab.H>
#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

class FArrayBox;

//
//@Man:
//@Memo: A Class Facilitating I/O for Fabs
/*@Doc:
       This data-less class aids I/O for FABs and encapsulates information
       about the floating point format being used in output.
       Note that the "new" format for writing out FABs is self-describing;
       i.e. we can always read in a FAB written in the "new" format.  For this
       reason, it is usually preferable to write FABs out in the native
       format on the machine, unless you're doing computations in 64 bit and
       only want to write out 32 bit FABs.

       With the exception of the enumeration constants, this class is
       primarily for FArrayBox implementors; i.e. user's shouldn't
       call any of the member functions in this class directly.
*/

class FABio
{
public:
    /*@ManDoc: An enum which controls precision of FAB output.
               Valid values are FAB\_FLOAT and FAB\_DOUBLE.  This
               is deprecated; i.e. please don't use it except
               for reading old FABs as it will probably be going
               away in a later release.
    */
    enum Precision
    {
        FAB_FLOAT = 0,
        FAB_DOUBLE
    };
    /*@ManDoc: An enum which controls format of FAB output.
               Valid values are FAB\_ASCII, FAB\_IEEE, FAB\_NATIVE,
               FAB\_8BIT and FAB\_IEEE\_32;

               FAB\_ASCII: write the FAB out in ASCII format.

               FAB\_8BIT: write the FAB out with all floating-point
               values scaled to range 0 - 255.

               FAB\_NATIVE: write out floating-point values in the native
               format.  This is usually the "best" choice of formats.

               FAB\_IEEE\_32: write out floating-point values in IEEE 32
               bit normal format.  This is recommended for use when your
               internal computations are done in 64 bits and you want to save
               space when writing out the FABs.

               FAB\_IEEE: this is deprecated.  It is identical to
               FAB\_IEEE\_32.
    */
    enum Format
    {
        FAB_ASCII = 0,
        FAB_IEEE,
        FAB_NATIVE,
        //
        // This is set to four so that when reading in an old FAB,
        // we don't get confused when we see an old FAB\_8BITRLE file.
        //
        FAB_8BIT = 4,
        FAB_IEEE_32
    };
    /*@ManDoc: An enum which controls byte ordering of FAB output.
               Valid values are FAB\_NORMAL\_ORDER, FAB\_REVERSE\_ORDER,
               and FAB\_REVERSE\_ORDER\_2.  This is deprecated; i.e. please
               don't use it except for reading old FABs as it will probably
               be going away in a later release.  These exist solely to
               describe the ordering of "old" FABs that you want to read.
    */
    enum Ordering
    {
        FAB_NORMAL_ORDER,
        FAB_REVERSE_ORDER,
        FAB_REVERSE_ORDER_2
    };

    //
    // None of the member functions here are documented with doc++
    // as this class really isn't meant for users.

    //
    //The virtual destructor.
    //
    virtual ~FABio ();
    //
    // Pure virtual function.  Derived classes MUST override this
    // function to read an FArrayBox from the istream, under the
    // assumption that the header has alread been read.
    //
    virtual void read (istream&   is,
                       FArrayBox& fb) const = 0;
    //
    // Pure virtual function.  Derived classes MUST override this
    // function to write the FArrayBox to the ostream, under the
    // assumption that the header for the FAB has already been
    // written.  Write it out as if it only had num_comp components
    // with component comp being the first one.
    //
    virtual void write (ostream&         os,
                        const FArrayBox& fb,
                        int              comp,
                        int              num_comp) const = 0;
    //
    // Pure virtual function.  Derived classes MUST override this
    // function to skip over the next FAB f in the istream, under the
    // assumption that the header for the FAB f has already been
    // skpped over.
    //
    virtual void skip (istream&   is,
                       FArrayBox& f) const = 0;
    //
    // Write out a header describing FArrayBox f that contains
    // nvar components.  It must be the case that nvar <= f.nComp().
    //
    virtual void write_header (ostream&         os,
                               const FArrayBox& f,
                               int              nvar) const;
    //
    // Read in the header from the istream.
    // Returns a new'd FABio of the written-out type.
    // Complements write_header.  The user is responsible
    // for delete'ing the returned FABio*.  The FArrayBox f is
    // resized to be on the Box and number of components read
    // in from the header file.  This is in preparation for
    // next doing a read.  This is split up so that we can make
    // the read functions virtual, while having all the code for
    // detailing the type of FArrayBox that was written out in one place.
    //
    static FABio* read_header (istream&   is,
                               FArrayBox& f);
};

//@Man:
//@Memo: A Fortran Array of REALs
/*@Doc: 

  Fortran Array Box's (generally called FAB's) are objects constructed
  to emulate the `FORTRAN' array.  Useful operations can be performed
  upon FAB's in C++, and they provide a convenient interface to
  `FORTRAN' when it is necessary to retreat into that language.

  FArrayBox is derived from ArithFab<Real>, NormedFab<Real> and
  OrderedFab<Real> and inherits many member functions from them.
  FArrayBox adds additional useful capabilities which make sense
  for Real types, such as I/O and L**p norms.

  FArrayBox's may be output in various format's with various
  precisions. Available formats are ascii, `IEEE 32', native, and
  8bit.  Ascii is self explanatory.  IEEE 32 is the standard on
  most microprocessor based computers; it specifies that output should be in
  IEEE 32 bit floating point format in normal order.  Native is the most
  efficient, but not necessarily most portable, I/O format for a
  given computer.  8bit is a 8 bit format, suitable for graphics output.
  In general, you probably want to use NATIVE, unless you need to write
  your data out with less precision than your computation in order to
  save space.

  The format and precision may be set in a file read by the ParmParse
  class by the "fab.format" variable.  Allowed values are NATIVE, ASCII,
  8BIT and IEEE32.

  FABs written using `operator<<' are always written in ASCII.
  FABS written using `writOn' use the FABio::Format specified with
  `setFormat' or the FABio::Format specified in the ParmParse file
  read by `init'. If the FABio::Format is not set explicitly by either
  of these two methods, then it defaults to NATIVE.

  The C pre-processor macro `BL\_SPACEDIM' must be defined to use
  this class.  The internal precision of FARRAYBOX objects is
  set by defining either `BL\_USE\_FLOAT' or `BL\_USE\_DOUBLE' 

  This class does NOT provide a copy constructor or assignment operator.  */

class FArrayBox
    :
    public ArithFab<Real>, 
    public OrderedFab<Real>,
    public NormedFab<Real>
{
    //
    // FABio is a friend of ours.
    //
    friend class FABio;
public:
    //
    //@ManDoc: Construct an invalid FAB with no memory.
    //
    FArrayBox ();
    //
    /*@ManDoc: Construct an initial FAB with the data space allocated but
               not inititialized. ncomp is the number of components
               (variables) at each data point in the Box.
    */
    explicit FArrayBox (const Box& b,
                        int        ncomp=1);

#ifdef BL_USE_POINTLIB
#ifndef BL_CRAY_BUG_DEFARG
 
    /*@ManDoc: Build a Fab From a PointFab.  This is defined only if 
               BL\_CRAY\_BUG\_DEFARG is not defined and if BL\_USE\_POINTLIB
               is defined.
    */
    FArrayBox (const PointFab<PointDomain>& pf,
               Real                         val);

#endif /*BL_CRAY_BUG_DEFARG*/
#endif /*BL_USE_POINTLIB*/

    //
    //@ManDoc: Constructs FAB by reading it from istream.
    //
    explicit FArrayBox (istream& is);
    //
    //@ManDoc:  The (virtual) destructor.
    //
    virtual ~FArrayBox ();

    /*@ManDoc: Compute the Lp-norm of this FAB using components
               (comp : comp+numcomp-1).  p < 0  -> ERROR.
               p = 0  -> infinity norm (max norm).
               p = 1  -> sum of ABS(FAB)
               p > 1  -> Lp-norm
    */
    virtual Real norm (int p=2,
                       int comp=0,
                       int numcomp=1) const;
    //
    //@ManDoc: Same as above except only on given subbox.
    //
    virtual Real norm (const Box& subbox,
                       int        p=2,
                       int        comp=0,
                       int        numcomp=1) const;
    //
    //@ManDoc: Print subrange of FAB in formatted form to output stream.
    //
    friend void printRange (ostream&,
                            const FArrayBox&,
                            const Box&,
                            int              comp,
                            int              numcomp = 1);
#if   BL_SPACEDIM==1

    /*@ManDoc: Print subrange defined by (ilo:ihi)
               using components (comp : comp+numcomp-1).
    */
    friend void printRange (ostream&,
                            const FArrayBox&,
                            int               ilo,
                            int               ihi,
                            int               comp,
                            int               numcomp=1);

    /*@ManDoc: Same as above but print to standard output.
               This is useful for debugging.
    */
    friend void printRange (const FArrayBox&,
                            int               ilo,
                            int               ihi,
                            int               comp,
                            int               numcomp=1);
#elif BL_SPACEDIM==2

    /*@ManDoc: Print subrange defined by (ilo:ihi) X (jlo:jhi)
               using components (comp : comp+numcomp-1).
    */
    friend void printRange (ostream&,
                            const FArrayBox&,
                            int               ilo,
                            int               ihi,
                            int               jlo,
                            int               jhi,
                            int               comp,
                            int               numcomp=1);
 
    /*@ManDoc: Same as above but print to standard output.
               This is useful for debugging.
    */
    friend void printRange (const FArrayBox&,
                            int               ilo,
                            int               ihi,
                            int               jlo,
                            int               jhi,
                            int               comp,
                            int               numcomp=1);
#elif BL_SPACEDIM==3

    /*@ManDoc: Print subrange defined by (ilo:ihi) X (jlo:jhi) X (klo:khi)
               using components (comp : comp+numcomp-1).
    */
    friend void printRange (ostream&,
                            const FArrayBox&,
                            int               ilo,
                            int               ihi,
                            int               jlo,
                            int               jhi,
                            int               klo,
                            int               khi,
                            int               comp,
                            int               numcomp=1);
 
    /*@ManDoc: Same as above but print to standard output.
               This is useful for debugging.
    */
    friend void printRange (const FArrayBox&,
                            int               ilo,
                            int               ihi,
                            int               jlo,
                            int               jhi,
                            int               klo,
                            int               khi,
                            int               comp,
                            int               numcomp=1);
#endif
    //
    //@ManDoc: Write FABs in ASCII form.
    //
    friend ostream& operator<< (ostream&         os,
                                const FArrayBox& fb);
    //
    //@ManDoc: Read FABs in ASCII form.
    //
    friend istream& operator>> (istream&   is,
                                FArrayBox& fb);

    /*@ManDoc: Writes out the FAB in whatever format you've set.
               The default format is ASCII.
    */
    void writeOn (ostream& of) const;
    //
    /*@ManDoc: Write only selected range of components.  comp specifies
               from which component (starting at 0) to write at each
               point in space.  num\_comp specifies how many data points
               to write out at each point is space -- it defaults to 1.
               It must be the case the comp >= 0 \&\& num\_comp >= 1 \&\&
               (comp+num\_comp) <= nComp().  The FAB is written out in
               whatever format you've set, with the default format being
               ASCII.  The FAB that is written to disk will be an
               num\_comp component FAB.
    */
    void writeOn (ostream& of,
                  int      comp,
                  int      num_comp=1) const;
    //
    //@ManDoc: Read FAB from istream.  Format is as it was written out.
    //
    void readFrom (istream& is);

    /*@ManDoc: Skip over the next FAB from the input stream.
               Return the Box defining the domain of the FAB and the
               number of components.
    */
    static Box skipFAB (istream& is,
                        int&     num_comp);
    //
    //@ManDoc: Skip over the next FAB from the input stream.
    //
    static void skipFAB (istream& is);

    /*@ManDoc: Set the FABio::Format in the program.
               This is the preferred way to set the output format
               in "new" FABs.  When desiging new programs, this should
               be the only function that needs to be called in order
               to set the format.
    */
    static void setFormat (FABio::Format fmt);
    //
    //@ManDoc: Gets the FABio::Format set in the program.
    //
    static FABio::Format getFormat ();

    /*@ManDoc: Set the FABio::Ordering for reading old FABs.  It does
               NOT set the ordering for output.
               This is deprecated.  It exists only to facilitate
               reading old FABs.  When you're reading in an "old" FAB,
               you must set the Ordering, before attempting
               to read it in.  This is because FABs written out in the
               "old" format weren't self-describing; i.e. information
               such as the Ordering was lost when the "old" FAB was
               written out.
    */
    static void setOrdering (FABio::Ordering ordering);

    /*@ManDoc: Gets the FABio::Ordering set in the program.  This is
               deprecated.  It does NOT do the right thing with the
               new FAB I/O format.
     */
    static FABio::Ordering getOrdering ();

    /*@ManDoc: Set the FABio::Precision.  This is deprecated.  It
               is not useful with the "new" FAB I/O format.
    */
    static void setPrecision (FABio::Precision precision);

    /*@ManDoc: Returns the FABio::Precision.  This is deprecated.  It
               is not useful with the "new" FAB I/O format.  Always
               returns FABio::Float.
    */
    static FABio::Precision getPrecision ();
    //
    //@ManDoc: Returns reference to the FABio object used by the program.
    //
    static const FABio& getFABio ();

    /*@ManDoc: Sets the FABio object used by the program.  It is an error
               if the passed pointer `rd' is the null pointer.
    */
    static void setFABio (FABio* rd);

protected:
    //
    // Initialize from ParmParse with "fab" prefix.
    //
    static void init ();
    //
    // Format and ordering for all FAB output.
    // This stuff exists solely to support reading old FABs.
    //
    static FABio::Format   format;
    static FABio::Ordering ordering;
    //
    // Has init() been called?
    //
    static bool Initialized;
    //
    // The FABio pointer describing our output format.
    //
    static FABio* fabio;

private:
    //
    // These are disallowed.
    //
    FArrayBox (const FArrayBox&);
    FArrayBox& operator= (const FArrayBox&);
};

//
//@ManDoc: A typedef for backward compatibility with old code.
//
typedef FArrayBox FARRAYBOX;

//
// Inlines.
//

inline
FABio::Format
FArrayBox::getFormat ()
{
    return format;
}

inline
const FABio&
FArrayBox::getFABio ()
{
    return *fabio;
}

inline
void
FArrayBox::setFABio (FABio* rd)
{
    assert(rd != 0);
    delete fabio;
    fabio = rd;
}

inline
Real
FArrayBox::norm (int p,
                 int comp,
                 int numcomp) const
{
    return norm(domain,p,comp,numcomp);
}

inline
void
FArrayBox::writeOn (ostream& os) const
{
    writeOn(os,0,nComp());
}

inline
void
FArrayBox::skipFAB (istream& is)
{
    int ignore = 0;
    skipFAB(is, ignore);
}

#endif /*BL_FARRAYBOX_H*/
