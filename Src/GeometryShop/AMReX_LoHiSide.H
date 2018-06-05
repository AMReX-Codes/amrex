#ifndef AMREX_LOHISIDE_H_
#define AMREX_LOHISIDE_H_

namespace amrex
{
  /// Encapsulation of high or low side.
  /**
     Encapsulation of high or low side.
  */
  class Side
  {
  public:
    ///
    /**
       Lo is the low side.<br>
       Hi is the high side.
    */
    enum LoHiSide
      {
        Invalid=-1,
        Lo=0,
        Hi,
        NUMSIDES
      };

    /**
       Constructor.
    */
    Side()
    {
    }

    /**
       Destructor.
    */
    ~Side()
    {
    }

    ///
    /**
       Returns the other side.  That is, swap Hi for Lo and vice versa.
    */
    friend Side::LoHiSide flip(const Side::LoHiSide& a_side);

    ///
    static Side::LoHiSide flip(const Side::LoHiSide& a_side);

    ///
    /**
       Returns -1 for Lo Side, +1 for Hi Side.
    */
    friend
    int sign(const Side::LoHiSide& a_side);
  };

  // needs a declaration outside the class (Stroustrup 11.5.1)
  int sign(const Side::LoHiSide& a_side);

  /// Iterator for low and high side.
  /**
     SideIterator iterates over the values Side::Lo and Side::Hi.
     There are no adjustable parameters whatsoever.
  */
  class SideIterator
  {
  public:
    ///
    /**
       Default, and only, constructor.
    */
    SideIterator();

    ///
    /**
       Default destructor.
    */
    ~SideIterator()
    {}

    /**
       \name Modification functions
    */
    /*@{*/

    ///
    /**
       Set the iterator to the first side.
    */
    void begin();

    ///
    /**
       Set the iterator to the first side.
    */
    void reset()
    {
      begin();
    }

    ///
    /**
       Advance the iterator to the next face.
    */
    void next();

    ///
    /**
       Advance the iterator to the next face.
    */
    void operator ++ ();

    /*@}*/

    /**
       \name Access functions
    */
    /*@{*/

    ///
    /**
       Returns the value of the side for the current iteration.
    */
    Side::LoHiSide operator () () const;

    ///
    /**
       Returns <tt>true</tt> if there is a next side, <tt>false</tt> if not.
    */
    bool ok() const;


  private:
    int m_current;
  };
}
#endif
