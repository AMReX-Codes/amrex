//BL_COPYRIGHT_NOTICE

#ifndef _TagBox_H_
#define _TagBox_H_

//
// $Id: TagBox.H,v 1.2 1997-11-24 18:52:33 lijewski Exp $
//

#include <IntVect.H>
#include <Box.H>
#include <Array.H>
#include <BaseFab.H>
#include <FabArray.H>
#include <BoxArray.H>
#include <BoxDomain.H>
#include <PArray.H>
#include <Geometry.H>

#if !defined(BL_SPACEDIM) || (BL_SPACEDIM != 1 && BL_SPACEDIM !=2 && BL_SPACEDIM != 3)
#error BL_SPACEDIM must be defined
#endif

#if (BL_SPACEDIM==1)
#define TAGBOX	TagBox1d
#elif (BL_SPACEDIM==2)
#define TAGBOX	TagBox2d
#elif (BL_SPACEDIM==3)
#define TAGBOX	TagBox3d
#endif


class TAGBOX
    : public BaseFab<int>
{
  public:
    enum TagVal{CLEAR=0, BUF, SET};

      // administrative
      // construct an invalid FAB with no memory
    TAGBOX();

    TAGBOX(const BOX &bx, int n = 1);

      // destructor
    virtual ~TAGBOX();
    
    void resize(const BOX& b) {BaseFab<int>::resize(b,1);}

    int& operator()(const INTVECT& p)
    {
	return dptr[domain.index(p)];
    }
    const int& operator()(const INTVECT& p) const
    {
	return dptr[domain.index(p)];
    }
 
      // construct and return a new tagbox in which the
      // coarsened cell is tagged of any of the corresponding
      // fine cells are tagged
    TAGBOX* coarsen(const IntVect & ratio);

      // mark neighbors of every tagged cell a distance nbuff away
      // only search interior for initial tagged points where nwid
      // is given as the width of the bndry region.
    void buffer(int nbuff, int nwid);

      // tag cells on intersect with src if corresponding src cell
      // is tagged.
    void merge(const TAGBOX& src);

      // add location of every tagged cell to INTVECT array,
      // starting at given location
    int colate(Array<INTVECT> &ar, int start) const;

#if 0
    void setVal(TagVal val) { BaseFab<int>::setVal(val,0); }
    void setVal(TagVal val, const BOX& bx) {
	BaseFab<int>::setVal(val,bx,0);
    }
#endif

    int numTags(const BOX& bx) const;
    int numTags() const;
  private:
    TAGBOX(const TAGBOX &);
    TAGBOX& operator = (const TAGBOX&);
};


class TagBoxArray : public FabArray<int, TAGBOX>
{
//protected:
    //PArray<TAGBOX> tb;
    int border;  // not to be confused with FabArray::n_grow which defines
		 // the invalid regions of the fabs in the FabArray
public:
    //TagBoxArray(PArrayPolicy _managed=PArrayManage)
	//: tb(_managed) {};
    //TagBoxArray(int len, PArrayPolicy _managed=PArrayManage)
	//: tb(len,_managed) {};
    //TagBoxArray(const BoxArray& bs, int _ngrow=0,
		//PArrayPolicy _managed = PArrayManage);

    TagBoxArray(const BoxArray& bs, int _ngrow=0);

    virtual ~TagBoxArray();

    int borderSize() { return border; }
    void buffer(int nbuf);
    void mergeUnique();

  // map tagged cells through a periodic boundary to other grids in TagBoxArray
  // cells which were outside domain are set to CLEAR
    void mapPeriodic( const Geometry& geom );

    void setVal(BoxDomain& bd, TAGBOX::TagVal val);
    void setVal(BoxArray& ba, TAGBOX::TagVal val);
    void coarsen(const IntVect & ratio);
    //int length() const { return tb.length(); }
    //TAGBOX& operator[](int i) {
	//return tb[i];
    //}
    //const TAGBOX& operator[](int i) const {
	//return tb[i];
    //}

    int numTags() const;

    Array<INTVECT>* colate() const;
private:
    TagBoxArray();
    TagBoxArray(const TagBoxArray &);
    TagBoxArray& operator = (const TagBoxArray &);
};

#endif
