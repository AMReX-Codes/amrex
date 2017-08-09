#ifndef AMREX_PRIMITIVE_TYPES
#define AMREX_PRIMITIVE_TYPES
//Question? email tannguyen@lbl.gov
//Created 07-19-2017
//Last modification 07-24-2017

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <string>
#include <sstream>
#include <cassert>
using namespace std;

namespace amrex{
    //!This class is similar to IntVect, but it supports template
    template <int D=1>
	class PointVect{
	    private:
		int _vect[D];
	    public:
		class shift_hasher{
		    private:
			static const unsigned shift=0;
			static const unsigned shift_stride=0;
			shift_hasher(){
			    shift=shift_stride= 8*sizeof(size_t)/D;
			}
		    public:
			size_t operator()(const PointVect& vec) const
			{
			    size_t ret=_vect[0];
			    for(int i=1; i<D; i++){
				ret ^= (_vect[i] << shift);
				shift+= shift_stride;
			    }
			    return ret;
			}
		};
		PointVect(){bzero(_vect, D*sizeof(int));}
		PointVect(int initVal){for(int i=0; i<D; i++) _vect[i]=initVal;};
		int& operator[] (int i) {
		    return _vect[i];
		}
		PointVect<D>& set(int i, int val){
		    _vect[i]=val;
		    return *this;
		}
		PointVect<D>& operator= (const PointVect<D>& rhs) const{
		    for(int i=0; i<D; i++) _vect[i]=rhs[i];
		    return *this;
		}
		bool operator== (int val) const{
		    for(int i=0; i<D; i++) if(_vect[i]!=val) return false;
		    return true;
		}
		bool operator!= (int val) const{
		    for(int i=0; i<D; i++) if(_vect[i]==val) return false;
		    return true;
		}
		bool operator== (const PointVect& rhs) const{
		    for(int i=0; i<D; i++) if(_vect[i]!=rhs[i]) return false;
		    return true;
		}
		bool operator!= (const PointVect& rhs) const{
		    for(int i=0; i<D; i++) if(_vect[i]==rhs[i]) return false;
		    return true;
		}
	};

    /*
       template <class T, int D=1>
       class TaskTable{
       private:
       std::unorderedMap<PointVect<D>, std::list<T*>, PointVect::shift_hasher>
       public:
       };
     */
}

#endif 
