
//
// $Id: BndryData.cpp,v 1.1 1998-03-24 07:05:16 almgren Exp $
//

#include <BndryData.H>
#include <Utility.H>

const char NL = '\n';
const char SP = ' ';

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();
#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

BndryData::BndryData(const BoxArray& _grids, int _ncomp, 
                     const ProxyGeometry& _geom)
    : BndryRegister(), geom(_geom)
{
    define(_grids,_ncomp,_geom);
}

// copy constructor
BndryData::BndryData(const BndryData& src)
{
    (*this) = src;
}

BndryData::BndryData(istream& is)
{
    readFrom(is);
}

BndryData::~BndryData()
{
      // masks was not allocated with PArrayManage, must manually dealloc
    clear_masks();
}

void
BndryData::clear_masks()
{
    for (int i = 0; i <2*BL_SPACEDIM; i++) {
        int len = masks[i].length();
        for (int k = 0; k <len; k++) {
            if (masks[i].defined(k)) {
                Mask *m = masks[i].remove(k);
                delete m;
            }

        }
    }
}

ostream& operator << (ostream& os, const BndryData &mgb)
{
ParallelDescriptor::Abort("BndryData::operator<< not yet implemented in parallel");
    const BoxArray& grds = mgb.boxes();
    int ngrds = grds.length();
    os << "[BndryData with " << ngrds << " grids:\n";
    for (int grd = 0; grd < ngrds; grd++){
        for (OrientationIter face; face; ++face) {
            Orientation f = face();
            os << "::: face " << f << " of grid " << grds[grd] << NL;
            os << "BC = " << mgb.bcond[f][grd]
               << " LOC = " << mgb.bcloc[f][grd] << NL;
            os << mgb.masks[f][grd];
            os << mgb.bndry[f][grd];
        }
        os << "--------------------------------------------------\n";
    }
    return os;
}

void
BndryData::writeOn(ostream &os) const
{
ParallelDescriptor::Abort("BndryData::writeOn not yet implemented in parallel");
    BndryRegister::writeOn(os);
    int len = grids.length();
    for (OrientationIter fi; fi; ++fi) {
        Orientation face = fi();
        for (int grd = 0; grd < len; grd++) {
            os << face << SP << grids[grd] << SP
               << bcond[face][grd] << SP << bcloc[face][grd] << NL;
            masks[face][grd].writeOn(os);
        }
    }
    os << NL;
}

void
BndryData::readFrom(istream& is)
{
ParallelDescriptor::Abort("BndryData::readFrom not yet implemented in parallel");
    BndryRegister::readFrom(is);
    int len = grids.length();
    for (OrientationIter fi; fi; ++fi) {
        Orientation face = fi();
        for (int grd = 0; grd < len; grd++) {
            Orientation face_in;
            int grd_in;
            is >> face_in;
            assert(face == face_in);
            is >> grd_in;
            assert(grd_in);
            int bcond_i;
            is >> bcond_i;
            bcond[face][grd] = bcond_i;
            is >> bcloc[face][grd];
            is.ignore(BL_IGNORE_MAX,'\n');
            masks[face][grd].readFrom(is);
        }
    }
    is.ignore(BL_IGNORE_MAX,'\n');
}


BndryData&
BndryData::operator = (const BndryData& src)
{
      // got to save the geometric info
    geom = src.geom;
   
      // redefine grids and bndry array
    BndryRegister::operator = ( (BndryRegister) src);
    int ngrd = grids.length();
    clear_masks();
    for (int i = 0; i < 2*BL_SPACEDIM; i++) {
        bcond[i].resize(ngrd);
        bcloc[i].resize(ngrd);
        masks[i].resize(ngrd);
        for (int grd = 0; grd < ngrd; grd++) {
            bcond[i][grd] = src.bcond[i][grd];
            bcloc[i][grd] = src.bcloc[i][grd];
            const Mask& src_mask = src.masks[i][grd];
            Mask *m = new Mask(src_mask.box(),src_mask.nComp());
            m->copy(src_mask);
            masks[i].set(grd,m);
        }
    }
    return *this;
}

void
BndryData::define(const BoxArray& _grids, int _ncomp, const ProxyGeometry& _geom)
{
    geom = _geom;
    BndryRegister::setBoxes(_grids);
    int len = grids.length();
    assert( len > 0 );

    for (OrientationIter fi; fi; ++fi) {
        Orientation face = fi();
        int coord_dir = face.coordDir();
        masks[face].resize(len);
        bcloc[face].resize(len);
        bcond[face].resize(len);

        BndryRegister::define(face,IndexType::TheCellType(),0,1,0,_ncomp);
        const FabSet& fabs = bndry[face];


          // alloc mask and set to quad_interp value
        for (int k = 0; k < len; k++) {
            Box face_box = adjCell(grids[k],face,1);
              // extend box in directions orthogonal to face normal
            for (int dir = 0; dir < BL_SPACEDIM; dir++) {
                if (dir == coord_dir) continue;
                face_box.grow(dir,1);
            }
            Mask *m = new Mask(face_box);
            m->setVal(outside_domain,0);
            Box dbox(geom.Domain());
            dbox &= face_box;
            m->setVal(not_covered,dbox,0);
            // now have to set as not_covered the periodic translates as well
            if( geom.isAnyPeriodic() ){
              Box dombox(geom.Domain());
              Array<IntVect> pshifts(27);
              geom.periodicShift( dombox, face_box, pshifts );
              for( int iiv=0; iiv<pshifts.length(); iiv++){
                IntVect iv = pshifts[iiv];
                m->shift(iv);
                Box target(dombox);
                target &= m->box();
                if (target.ok()) m->setVal(not_covered,target,0);
                m->shift(-iv);
              }
            }
            masks[face].set(k,m);
              // turn mask off on intersection with grids at this level
            for (int g = 0; g < len; g++) {
                Box ovlp(grids[g]);
                ovlp &= face_box;
                if (ovlp.ok()) m->setVal(covered,ovlp,0);
            }
            // handle special cases if is periodic
            if( geom.isAnyPeriodic() && 
                !geom.Domain().contains(face_box) ){
              Array<IntVect> pshifts(27);
              geom.periodicShift( geom.Domain(), face_box, pshifts);
              for( int iiv=0; iiv<pshifts.length(); iiv++ ){
                IntVect iv = pshifts[iiv];
                m->shift(iv);
                for( int g=0; g<len; g++){
                  Box ovlp(grids[g]);
                  ovlp &= m->box();
                  if( ovlp.ok() ) m->setVal(covered,ovlp,0);
                }
                m->shift(-iv);
              }
            }
        }
    }
}

