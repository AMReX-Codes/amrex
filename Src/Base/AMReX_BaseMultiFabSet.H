#include <AMReX_FabArray.H>

namespace amrex {

template <class MF, std::size_t N>
class BaseFabArraySet {
    /**
    * \brief Return a constant reference to the BoxArray that defines the
    * valid region associated with this FabArray.
    */
    [[nodiscard]] const BoxArray& boxArray () const noexcept { return boxarray; }


}




// template<class T, std::size_t N>
// class BaseMultiFabSet {};

// template <class MF, std::size_t N>
// class BaseMultiFabSet<MF, typename std::enable_if<std::is_base_of<FabArrayBase>::value>::type, N> : public std::array<MF, N> 
// {
//     FabArrayBase x();
// };

// // the partial specialization of A is enabled via a template parameter
// template<class T, class Enable = void, std::size_t N>
// class A {}; // primary template
 
// template<class T, std::size_t N>
// class A<T, typename std::enable_if<std::is_floating_point<T>::value>::type, N>
// {
// }; // specialization for floating point types

template<class T, std::size_t N>
class A
{
    public:
    static void TestFunction ();
};


class B : public A
{
    static void TestFunction2 ();

    static void Test2 ()
    {
        B::TestFunction();
    }    
};

class C {
    static void Test3 ();
};

};