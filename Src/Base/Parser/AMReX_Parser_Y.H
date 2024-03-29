#ifndef AMREX_PARSER_Y_H_
#define AMREX_PARSER_Y_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_Math.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>

void amrex_parsererror (char const *s, ...);

namespace amrex {

enum parser_f1_t {  // Built-in functions with one argument
    PARSER_SQRT,
    PARSER_EXP,
    PARSER_LOG,
    PARSER_LOG10,
    PARSER_SIN,
    PARSER_COS,
    PARSER_TAN,
    PARSER_ASIN,
    PARSER_ACOS,
    PARSER_ATAN,
    PARSER_SINH,
    PARSER_COSH,
    PARSER_TANH,
    PARSER_ASINH,
    PARSER_ACOSH,
    PARSER_ATANH,
    PARSER_ABS,
    PARSER_FLOOR,
    PARSER_CEIL,
    PARSER_COMP_ELLINT_1,
    PARSER_COMP_ELLINT_2
};

static
#if defined(__INTEL_COMPILER) && defined(__EDG__)
    const
#else
    constexpr
#endif
std::string_view parser_f1_s[] =
{
    "sqrt",
    "exp",
    "log",
    "log10",
    "sin",
    "cos",
    "tan",
    "asin",
    "acos",
    "atan",
    "sinh",
    "cosh",
    "tanh",
    "asinh",
    "acosh",
    "atanh",
    "abs",
    "floor",
    "ceil",
    "comp_ellint_1",
    "comp_ellint_2"
};

enum parser_f2_t {  // Built-in functions with two arguments
    PARSER_POW,
    PARSER_ATAN2,
    PARSER_GT,
    PARSER_LT,
    PARSER_GEQ,
    PARSER_LEQ,
    PARSER_EQ,
    PARSER_NEQ,
    PARSER_AND,
    PARSER_OR,
    PARSER_HEAVISIDE,
    PARSER_JN,
    PARSER_MIN,
    PARSER_MAX,
    PARSER_FMOD
};

static
#if defined(__INTEL_COMPILER) && defined(__EDG__)
    const
#else
    constexpr
#endif
std::string_view parser_f2_s[] =
{
    "pow",
    "atan2",
    "gt",
    "lt",
    "geq",
    "leq",
    "eq",
    "neq",
    "and",
    "or",
    "heaviside",
    "jn",
    "min",
    "max",
    "fmod"
};

enum parser_f3_t { // functions with three arguments
    PARSER_IF
};

static
#if defined(__INTEL_COMPILER) && defined(__EDG__)
    const
#else
    constexpr
#endif
std::string_view parser_f3_s[] =
{
    "if"
};

enum parser_node_t {
    PARSER_NUMBER,
    PARSER_SYMBOL,
    PARSER_ADD,
    PARSER_SUB,
    PARSER_MUL,
    PARSER_DIV,
    PARSER_F1,
    PARSER_F2,
    PARSER_F3,
    PARSER_ASSIGN,
    PARSER_LIST
};

static
#if defined(__INTEL_COMPILER) && defined(__EDG__)
    const
#else
    constexpr
#endif
std::string_view parser_node_s[] =
{
    "number",
    "symbol",
    "add",
    "sub",
    "mul",
    "div",
    "f1",
    "f2",
    "f3",
    "assign",
    "list"
};

/* In C, the address of the first member of a struct is the same as
 * the address of the struct itself.  Because of this, all struct parser_*
 * pointers can be passed around as struct parser_node pointer and enum
 * parser_node_t type can be safely checked to determine their real type.
 */

struct parser_node {
    enum parser_node_t type;
    enum parser_node_t padding;
    struct parser_node* l; // NOLINT(misc-confusable-identifiers)
    struct parser_node* r;
    struct parser_node* padding2;
};

struct alignas(parser_node) parser_number {
    enum parser_node_t type;
    double value;
};

struct alignas(parser_node) parser_symbol {
    enum parser_node_t type;
    char* name;
    int ip;
};

struct alignas(parser_node) parser_f1 {  /* Builtin functions with one argument */
    enum parser_node_t type;
    enum parser_f1_t ftype;
    struct parser_node* l; // NOLINT(misc-confusable-identifiers)
    struct parser_node* padding1;
    struct parser_node* padding2;
};

struct alignas(parser_node) parser_f2 {  /* Builtin functions with two arguments */
    enum parser_node_t type;
    enum parser_f2_t ftype;
    struct parser_node* l; // NOLINT(misc-confusable-identifiers)
    struct parser_node* r;
    struct parser_node* padding;
};

struct alignas(parser_node) parser_f3 { /* Builtin functions with three arguments */
    enum parser_node_t type;
    enum parser_f3_t ftype;
    struct parser_node* n1;
    struct parser_node* n2;
    struct parser_node* n3;
};

struct alignas(parser_node) parser_assign {
    enum parser_node_t type;
    struct parser_symbol* s;
    struct parser_node* v;
};

static_assert(sizeof(parser_f3) <= sizeof(parser_node), "amrex parser: sizeof parser_node too small");

/*******************************************************************/

/* These functions are used in bison rules to generate the original AST. */
void parser_defexpr (struct parser_node* body);
struct parser_symbol* parser_makesymbol (char* name);
struct parser_node* parser_newnode (enum parser_node_t type, struct parser_node* l,
                                    struct parser_node* r);
struct parser_node* parser_newneg (struct parser_node* n);
struct parser_node* parser_newnumber (double d);
struct parser_node* parser_newsymbol (struct parser_symbol* sym);
struct parser_node* parser_newf1 (enum parser_f1_t ftype, struct parser_node* l);
struct parser_node* parser_newf2 (enum parser_f2_t ftype, struct parser_node* l,
                                  struct parser_node* r);
struct parser_node* parser_newf3 (enum parser_f3_t ftype, struct parser_node* n1,
                                  struct parser_node* n2, struct parser_node* n3);
struct parser_node* parser_newassign (struct parser_symbol* s, struct parser_node* v);
struct parser_node* parser_newlist (struct parser_node* nl, struct parser_node* nr);

/*******************************************************************/

/* This is our struct for storing AST in a more packed way.  The whole
 * tree is stored in a contiguous chunk of memory starting from void*
 * p_root with a size of sz_mempool.
 */
struct amrex_parser {
    void* p_root;
    void* p_free;
    struct parser_node* ast;
    std::size_t sz_mempool;
};

struct amrex_parser* amrex_parser_new ();
void amrex_parser_delete (struct amrex_parser* parser);

struct amrex_parser* parser_dup (struct amrex_parser* source);
struct parser_node* parser_ast_dup (struct amrex_parser* parser, struct parser_node* node, int move);

void parser_regvar (struct amrex_parser* parser, char const* name, int i);
void parser_setconst (struct amrex_parser* parser, char const* name, double c);
void parser_print (struct amrex_parser* parser);
std::set<std::string> parser_get_symbols (struct amrex_parser* parser);
int parser_depth (struct amrex_parser* parser);

/* We need to walk the tree in these functions */
void parser_ast_optimize (struct parser_node* node);
std::size_t parser_ast_size (struct parser_node* node);
void parser_ast_print (struct parser_node* node, std::string const& space, std::ostream& printer);
void parser_ast_regvar (struct parser_node* node, char const* name, int i);
void parser_ast_setconst (struct parser_node* node, char const* name, double c);
void parser_ast_get_symbols (struct parser_node* node, std::set<std::string>& symbols,
                             std::set<std::string>& local_symbols);
int parser_ast_depth (struct parser_node* node);
void parser_ast_sort (struct parser_node* node);

/*******************************************************************/
double parser_get_number (struct parser_node* node);
void parser_set_number (struct parser_node* node, double v);
bool parser_node_equal (struct parser_node* a, struct parser_node* b);
/*******************************************************************/

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_exp (T a) { return std::exp(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_log (T a) { return std::log(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_log10 (T a) { return std::log10(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_sin (T a) { return std::sin(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_cos (T a) { return std::cos(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_tan (T a) { return std::tan(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_asin (T a) { return std::asin(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_acos (T a) { return std::acos(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_atan (T a) { return std::atan(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_sinh (T a) { return std::sinh(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_cosh (T a) { return std::cosh(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_tanh (T a) { return std::tanh(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_asinh (T a) { return std::asinh(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_acosh (T a) { return std::acosh(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_atanh (T a) { return std::atanh(a); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_comp_ellint_1 (T a)
{
#if defined(__GNUC__) && !defined(__clang__) && !defined(__CUDA_ARCH__) && !defined(__NVCOMPILER)
    return std::comp_ellint_1(a);
#else
    amrex::ignore_unused(a);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "parser: comp_ellint_1 only supported with gcc and cpu");
    return 0.0;
#endif
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_comp_ellint_2 (T a)
{
#if defined(__GNUC__) && !defined(__clang__) && !defined(__CUDA_ARCH__) && !defined(__NVCOMPILER)
    return std::comp_ellint_2(a);
#else
    amrex::ignore_unused(a);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "parser: comp_ellint_2 only supported with gcc and cpu");
    return 0.0;
#endif
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_pow (T a, T b) { return std::pow(a,b); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_atan2 (T a, T b) { return std::atan2(a,b); }

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
T parser_math_jn (int a, T b)
{
#if defined AMREX_USE_SYCL || defined __MINGW32__
    amrex::ignore_unused(a,b);
    // neither jn(f) nor std::cyl_bessel_j work yet
    // https://github.com/oneapi-src/oneAPI-spec/issues/308
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false, "parser: jn in SYCL not supported yet");
    return 0.0;
#else
    return jn(a, b);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
parser_call_f1 (enum parser_f1_t type, double a)
{
    switch (type) {
    case PARSER_SQRT:        return std::sqrt(a);
    case PARSER_EXP:         return parser_math_exp<double>(a);
    case PARSER_LOG:         return parser_math_log<double>(a);
    case PARSER_LOG10:       return parser_math_log10<double>(a);
    case PARSER_SIN:         return parser_math_sin<double>(a);
    case PARSER_COS:         return parser_math_cos<double>(a);
    case PARSER_TAN:         return parser_math_tan<double>(a);
    case PARSER_ASIN:        return parser_math_asin<double>(a);
    case PARSER_ACOS:        return parser_math_acos<double>(a);
    case PARSER_ATAN:        return parser_math_atan<double>(a);
    case PARSER_SINH:        return parser_math_sinh<double>(a);
    case PARSER_COSH:        return parser_math_cosh<double>(a);
    case PARSER_TANH:        return parser_math_tanh<double>(a);
    case PARSER_ASINH:       return parser_math_asinh<double>(a);
    case PARSER_ACOSH:       return parser_math_acosh<double>(a);
    case PARSER_ATANH:       return parser_math_atanh<double>(a);
    case PARSER_ABS:         return std::abs(a);
    case PARSER_FLOOR:       return std::floor(a);
    case PARSER_CEIL:        return std::ceil(a);
    case PARSER_COMP_ELLINT_1: return parser_math_comp_ellint_1<double>(a);
    case PARSER_COMP_ELLINT_2: return parser_math_comp_ellint_2<double>(a);
    default:
        amrex::Abort("parser_call_f1: Unknown function ");
        return 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
parser_call_f2 (enum parser_f2_t type, double a, double b)
{
    switch (type) {
    case PARSER_POW:
        return parser_math_pow<double>(a,b);
    case PARSER_ATAN2:
        return parser_math_atan2<double>(a,b);
    case PARSER_GT:
        return (a > b) ? 1.0 : 0.0;
    case PARSER_LT:
        return (a < b) ? 1.0 : 0.0;
    case PARSER_GEQ:
        return (a >= b) ? 1.0 : 0.0;
    case PARSER_LEQ:
        return (a <= b) ? 1.0 : 0.0;
    case PARSER_EQ:
        return (a == b) ? 1.0 : 0.0;
    case PARSER_NEQ:
        return (a != b) ? 1.0 : 0.0;
    case PARSER_AND:
        return ((a != 0.0) && (b != 0.0)) ? 1.0 : 0.0;
    case PARSER_OR:
        return ((a != 0.0) || (b != 0.0)) ? 1.0 : 0.0;
    case PARSER_HEAVISIDE:
        return (a < 0.0) ? 0.0 : ((a > 0.0) ? 1.0 : b);
    case PARSER_JN:
        return parser_math_jn<double>(int(a),b);
    case PARSER_MIN:
        return (a < b) ? a : b;
    case PARSER_MAX:
        return (a > b) ? a : b;
    case PARSER_FMOD:
        return std::fmod(a,b);
    default:
        amrex::Abort("parser_call_f2: Unknown function");
        return 0.0;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double
parser_call_f3 (enum parser_f3_t /*type*/, double a, double b, double c)
{
    // There is only one type currently
    return (a != 0.0) ? b : c;
}

}

#endif
