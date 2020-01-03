
// Used by EB_LSCore only.  Might be removed in the future.  So don't use.

namespace amrex {

    //! General collection of MultiFab utilities
    namespace MFUtil {

        //! Template class defining a copy operation where source and target are
        //! grown by ng. Template parameter T must be either MultiFab or
        //! iMultiFab
        template<typename T>
        struct SymmetricGhost {
            static void copy (T & target, const T & source, int nc, int ng) {
                target.copy(source, 0, 0, nc, ng, ng);
            }
        };


        //! Template class defining a copy operation where only the target is
        //! grown by ng. Template parameter T must be either MultiFab or
        //! iMultiFab
        template<typename T> // T can be either MultiFab or iMultiFab
        struct AsymmetricGhost {
            static void copy (T & target, const T & source, int nc, int ng) {
                target.copy(source, 0, 0, nc, 0, ng);
            }
        };


        //! Duplicate class T (which must be either of type MultiFab or
        //! iMultiFab) using a new BoxArray and DistributionMapping. Data from
        //! the source (mf_in) is copied using a ParallelCopy (ghost cells are
        //! treated based on the Cpy class (which is either SymmetricGhost or
        //! AsymmetricGhost)).
        template<typename T, template<typename> class Cpy>
        T duplicate(const BoxArray & ba, const DistributionMapping & dm, const T & mf_in) {

            int ng = mf_in.nGrow();
            int nc = mf_in.nComp();

            T mf_out(ba, dm, nc, ng);

            Cpy<T>::copy(mf_out, mf_in, nc, ng);

            return mf_out;
        }


        //! Duplicate class T (which must be either of type MultiFab or
        //! iMultiFab) using a new BoxArray and DistributionMapping. Data from
        //! the source (mf_in) is copied using a ParallelCopy (ghost cells are
        //! treated based on the Cpy class (which is either SymmetricGhost or
        //! AsymmetricGhost)). Following the copy operation, a FillBoundary
        //! operation is called.
        template<typename T, template<typename> class Cpy>
        T duplicate(const BoxArray & ba, const DistributionMapping & dm, const T & mf_in,
                    const Periodicity & periodicity) {

            T mf_out = duplicate<T, Cpy>(ba, dm, mf_in);
            mf_out.FillBoundary(periodicity);

            return mf_out;
        }


        //! Regrid (by duplication) class T (which must be either of type
        //! MultiFab or iMultiFab) usign the new BoxArray and
        //! DistributionMapping. The boolean flag `regrid_ghost` switched
        //! between the SymmetricGhost and AsymmetricGhost copy functions.
        template<typename T>
        T regrid(const BoxArray & ba, const DistributionMapping & dm,
                 const T & mf_in, bool regrid_ghost = false ) {

            if (regrid_ghost) {
                T mf_out = duplicate<T, SymmetricGhost>(ba, dm, mf_in);
                return mf_out;
            } else {
                T mf_out = duplicate<T, AsymmetricGhost>(ba, dm, mf_in);
                return mf_out;
            }
        }


        //! Regrid (by duplication) class T (which must be either of type
        //! MultiFab or iMultiFab) using the new BoxArray and
        //! DistributionMapping. The T-object at the end of the pointer `mf_out`
        //! is re-defined. The boolean flag `regrid_ghost` switched between the
        //! SymmetricGhost and AsymmetricGhost copy functions.
        template<typename T>
        void regrid( T & mf_out, const BoxArray & ba, const DistributionMapping & dm,
                     const T & mf_in, bool regrid_ghost = false ) {

            int ng = mf_in.nGrow();
            int nc = mf_in.nComp();

            mf_out.define(ba, dm, nc, ng);
            if (regrid_ghost)  SymmetricGhost<T>::copy(mf_out, mf_in, nc, ng);
            else              AsymmetricGhost<T>::copy(mf_out, mf_in, nc, ng);
        }


        //! Regrid (by duplication) class T. The T-object at the end of the
        //! reference `mf_out` is re-defined using the EBFArrayBoxFactory
        //! provided by `eb_factory`. The boolean flag `regrid_ghost` switched
        //! between the SymmetricGhost and AsymmetricGhost copy functions.
        template<typename T>
        void regrid( T & mf_out, const BoxArray & ba, const DistributionMapping & dm,
                     const FabFactory<FArrayBox> & eb_factory,
                     const T & mf_in, bool regrid_ghost = false ) {

            int ng = mf_in.nGrow();
            int nc = mf_in.nComp();

            mf_out.define(ba, dm, nc, ng, MFInfo(), eb_factory);
            if (regrid_ghost)  SymmetricGhost<T>::copy(mf_out, mf_in, nc, ng);
            else              AsymmetricGhost<T>::copy(mf_out, mf_in, nc, ng);
        }


        //! Convert class T_src (which must either be of type MultiFab or
        //! iMultiFab) to type T_dest (which must also be of type MultiFab or
        //! iMultiFab). This function is useful for convering iMultiFab MultiFab
        //! or vice-versa.
        template<typename T_src, typename T_dest>
            std::unique_ptr<T_dest> convert(const T_src & mf_in) {

            const BoxArray & ba            = mf_in.boxArray();
            const DistributionMapping & dm = mf_in.DistributionMap();
            int ncomp                      = mf_in.nComp();
            int ngrow                      = mf_in.nGrow();

            std::unique_ptr<T_dest> mf_out(new T_dest(ba, dm, ncomp, ngrow));

#ifdef _OPENMP
#pragma omp parallel
#endif
            for(MFIter mfi(mf_in, true); mfi.isValid(); ++ mfi) {
                const IArrayBox & in_tile  =      mf_in[mfi];
                FArrayBox & out_tile = (* mf_out)[mfi];

                for(BoxIterator bit(mfi.growntilebox()); bit.ok(); ++ bit) {

                    for(int i = 0; i < ncomp; i ++)
                        out_tile(bit(), i) = in_tile(bit(), i);

                }
            }

            return mf_out;
        }
    }

}
