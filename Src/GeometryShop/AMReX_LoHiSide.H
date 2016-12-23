

/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */


#ifndef _LOHISIDE_H_
#define _LOHISIDE_H_

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
