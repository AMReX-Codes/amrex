#ifndef AMREX_PARSER_EXE_H_
#define AMREX_PARSER_EXE_H_
#include <AMReX_Config.H>

#include <AMReX_Parser_Y.H>
#include <AMReX_Vector.H>

#ifndef AMREX_PARSER_STACK_SIZE
#define AMREX_PARSER_STACK_SIZE 16
#endif

#define AMREX_PARSER_LOCAL_IDX0 1000
#define AMREX_PARSER_GET_DATA(i) ((i)>=1000) ? pstack[(i)-1000] : x[i]

namespace amrex {

// N: node
// P: pointer offset
// V: value (i.e., double literal)

enum parser_exe_t {
    PARSER_EXE_NULL = 0,
    PARSER_EXE_NUMBER, // 1
    PARSER_EXE_SYMBOL, // 2
    PARSER_EXE_ADD,    // 3
    PARSER_EXE_SUB,    // 4
    PARSER_EXE_MUL,    // 5
    PARSER_EXE_DIV_F,  // 6
    PARSER_EXE_DIV_B,  // 7
    PARSER_EXE_NEG,    // 8
    PARSER_EXE_F1,     // 9
    PARSER_EXE_F2_F,   // 10
    PARSER_EXE_F2_B,   // 11
    PARSER_EXE_ADD_VP, // 12
    PARSER_EXE_SUB_VP, // 13
    PARSER_EXE_MUL_VP, // 14
    PARSER_EXE_DIV_VP, // 15
    PARSER_EXE_ADD_PP, // 16
    PARSER_EXE_SUB_PP, // 17
    PARSER_EXE_MUL_PP, // 18
    PARSER_EXE_DIV_PP, // 19
    PARSER_EXE_NEG_P,  // 20
    PARSER_EXE_ADD_VN, // 21
    PARSER_EXE_SUB_VN, // 22
    PARSER_EXE_MUL_VN, // 23
    PARSER_EXE_DIV_VN, // 24
    PARSER_EXE_ADD_PN, // 25
    PARSER_EXE_SUB_PN, // 26
    PARSER_EXE_MUL_PN, // 27
    PARSER_EXE_DIV_PN, // 28
    PARSER_EXE_IF,     // 29
    PARSER_EXE_JUMP    // 30
};

struct alignas(8) ParserExeNull {
    enum parser_exe_t type = PARSER_EXE_NULL;
};

struct alignas(8) ParserExeNumber {
    enum parser_exe_t type = PARSER_EXE_NUMBER;
    double v;
};

struct alignas(8) ParserExeSymbol {
    enum parser_exe_t type = PARSER_EXE_SYMBOL;
    int i;
};

struct alignas(8) ParserExeADD {
    enum parser_exe_t type = PARSER_EXE_ADD;
};

struct alignas(8) ParserExeSUB {
    enum parser_exe_t type = PARSER_EXE_SUB;
    double sign;
};

struct alignas(8) ParserExeMUL {
    enum parser_exe_t type = PARSER_EXE_MUL;
};

struct alignas(8) ParserExeDIV_F {
    enum parser_exe_t type = PARSER_EXE_DIV_F;
};

struct alignas(8) ParserExeDIV_B {
    enum parser_exe_t type = PARSER_EXE_DIV_B;
};

struct alignas(8) ParserExeNEG {
    enum parser_exe_t type = PARSER_EXE_NEG;
};

struct alignas(8) ParserExeF1 {
    enum parser_exe_t type = PARSER_EXE_F1;
    parser_f1_t ftype;
};

struct alignas(8) ParserExeF2_F {
    enum parser_exe_t type = PARSER_EXE_F2_F;
    parser_f2_t ftype;
};

struct alignas(8) ParserExeF2_B {
    enum parser_exe_t type = PARSER_EXE_F2_B;
    parser_f2_t ftype;
};

struct alignas(8) ParserExeADD_VP {
    enum parser_exe_t type = PARSER_EXE_ADD_VP;
    int i;
    double v;
};

struct alignas(8) ParserExeSUB_VP {
    enum parser_exe_t type = PARSER_EXE_SUB_VP;
    int i;
    double v;
};

struct alignas(8) ParserExeMUL_VP {
    enum parser_exe_t type = PARSER_EXE_MUL_VP;
    int i;
    double v;
};

struct alignas(8) ParserExeDIV_VP {
    enum parser_exe_t type = PARSER_EXE_DIV_VP;
    int i;
    double v;
};

struct alignas(8) ParserExeADD_PP {
    enum parser_exe_t type = PARSER_EXE_ADD_PP;
    int i1;
    int i2;
};

struct alignas(8) ParserExeSUB_PP {
    enum parser_exe_t type = PARSER_EXE_SUB_PP;
    int i1;
    int i2;
};

struct alignas(8) ParserExeMUL_PP {
    enum parser_exe_t type = PARSER_EXE_MUL_PP;
    int i1;
    int i2;
};

struct alignas(8) ParserExeDIV_PP {
    enum parser_exe_t type = PARSER_EXE_DIV_PP;
    int i1;
    int i2;
};

struct alignas(8) ParserExeNEG_P {
    enum parser_exe_t type = PARSER_EXE_NEG_P;
    int i;
};

struct alignas(8) ParserExeADD_VN {
    enum parser_exe_t type = PARSER_EXE_ADD_VN;
    double v;
};

struct alignas(8) ParserExeSUB_VN {
    enum parser_exe_t type = PARSER_EXE_SUB_VN;
    double v;
};

struct alignas(8) ParserExeMUL_VN {
    enum parser_exe_t type = PARSER_EXE_MUL_VN;
    double v;
};

struct alignas(8) ParserExeDIV_VN {
    enum parser_exe_t type = PARSER_EXE_DIV_VN;
    double v;
};

struct alignas(8) ParserExeADD_PN {
    enum parser_exe_t type = PARSER_EXE_ADD_PN;
    int i;
};

struct alignas(8) ParserExeSUB_PN {
    enum parser_exe_t type = PARSER_EXE_SUB_PN;
    int i;
    double sign;
};

struct alignas(8) ParserExeMUL_PN {
    enum parser_exe_t type = PARSER_EXE_MUL_PN;
    int i;
};

struct alignas(8) ParserExeDIV_PN {
    enum parser_exe_t type = PARSER_EXE_DIV_PN;
    int i;
    bool reverse;
};

struct alignas(8) ParserExeIF {
    enum parser_exe_t type = PARSER_EXE_IF;
    int offset;
};

struct alignas(8) ParserExeJUMP {
    enum parser_exe_t type = PARSER_EXE_JUMP;
    int offset;
};

template <int N>
struct ParserStack
{
    double m_data[N];
    int m_size = 0;
    constexpr void push (double v) { m_data[m_size++] = v; }
    constexpr void pop () { --m_size; }
    [[nodiscard]] constexpr double const& top () const { return m_data[m_size-1]; }
    [[nodiscard]] constexpr double      & top ()       { return m_data[m_size-1]; }
    [[nodiscard]] constexpr double operator[] (int i) const { return m_data[i]; }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double parser_exe_eval (const char* p, double const* x)
{
    ParserStack<AMREX_PARSER_STACK_SIZE> pstack;
    while (*((parser_exe_t*)p) != PARSER_EXE_NULL) { // NOLINT
        switch (*((parser_exe_t*)p))
        {
        case PARSER_EXE_NUMBER:
        {
            pstack.push(((ParserExeNumber*)p)->v);
            p   += sizeof(ParserExeNumber);
            break;
        }
        case PARSER_EXE_SYMBOL:
        {
            int i = ((ParserExeSymbol*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(d);
            p     += sizeof(ParserExeSymbol);
            break;
        }
        case PARSER_EXE_ADD:
        {
            double b = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() += b; // NOLINT
            p += sizeof(ParserExeADD);
            break;
        }
        case PARSER_EXE_SUB:
        {
            double b = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() = (pstack.top() - b) * (((ParserExeSUB*)p)->sign); // NOLINT
            p += sizeof(ParserExeSUB);
            break;
        }
        case PARSER_EXE_MUL:
        {
            double b = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() *= b; // NOLINT
            p += sizeof(ParserExeMUL);
            break;
        }
        case PARSER_EXE_DIV_F:
        {
            double v = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() /= v; // NOLINT
            p += sizeof(ParserExeDIV_F);
            break;
        }
        case PARSER_EXE_DIV_B:
        {
            double v = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() = v / pstack.top(); // NOLINT
            p += sizeof(ParserExeDIV_B);
            break;
        }
        case PARSER_EXE_NEG:
        {
            pstack.top() = -pstack.top(); // NOLINT
            p += sizeof(ParserExeNEG);
            break;
        }
        case PARSER_EXE_F1:
        {
            pstack.top() = parser_call_f1(((ParserExeF1*)p)->ftype, pstack.top()); // NOLINT
            p += sizeof(ParserExeF1);
            break;
        }
        case PARSER_EXE_F2_F:
        {
            double v = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() = parser_call_f2(((ParserExeF2_F*)p)->ftype, pstack.top(), v); // NOLINT
            p += sizeof(ParserExeF2_F);
            break;
        }
        case PARSER_EXE_F2_B:
        {
            double v = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() = parser_call_f2(((ParserExeF2_B*)p)->ftype, v, pstack.top()); // NOLINT
            p += sizeof(ParserExeF2_B);
            break;
        }
        case PARSER_EXE_ADD_VP:
        {
            int i = ((ParserExeADD_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(((ParserExeADD_VP*)p)->v + d);
            p     += sizeof(ParserExeADD_VP);
            break;
        }
        case PARSER_EXE_SUB_VP:
        {
            int i = ((ParserExeSUB_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(((ParserExeSUB_VP*)p)->v - d);
            p     += sizeof(ParserExeSUB_VP);
            break;
        }
        case PARSER_EXE_MUL_VP:
        {
            int i = ((ParserExeMUL_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(((ParserExeMUL_VP*)p)->v * d);
            p     += sizeof(ParserExeMUL_VP);
            break;
        }
        case PARSER_EXE_DIV_VP:
        {
            int i = ((ParserExeDIV_VP*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(((ParserExeDIV_VP*)p)->v / d);
            p     += sizeof(ParserExeDIV_VP);
            break;
        }
        case PARSER_EXE_ADD_PP:
        {
            int i = ((ParserExeADD_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeADD_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            pstack.push(d1+d2);
            p     += sizeof(ParserExeADD_PP);
            break;
        }
        case PARSER_EXE_SUB_PP:
        {
            int i = ((ParserExeSUB_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeSUB_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            pstack.push(d1-d2);
            p     += sizeof(ParserExeSUB_PP);
            break;
        }
        case PARSER_EXE_MUL_PP:
        {
            int i = ((ParserExeMUL_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeMUL_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            pstack.push(d1*d2);
            p     += sizeof(ParserExeMUL_PP);
            break;
        }
        case PARSER_EXE_DIV_PP:
        {
            int i = ((ParserExeDIV_PP*)p)->i1;
            double d1 = AMREX_PARSER_GET_DATA(i);
            i     = ((ParserExeDIV_PP*)p)->i2;
            double d2 = AMREX_PARSER_GET_DATA(i);
            pstack.push(d1/d2);
            p      += sizeof(ParserExeDIV_PP);
            break;
        }
        case PARSER_EXE_NEG_P:
        {
            int i = ((ParserExeNEG_P*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(-d);
            p     += sizeof(ParserExeNEG_P);
            break;
        }
        case PARSER_EXE_ADD_VN:
        {
            pstack.top() += ((ParserExeADD_VN*)p)->v; // NOLINT
            p       += sizeof(ParserExeADD_VN);
            break;
        }
        case PARSER_EXE_SUB_VN:
        {
            pstack.top() = ((ParserExeSUB_VN*)p)->v - pstack.top(); // NOLINT
            p      += sizeof(ParserExeSUB_VN);
            break;
        }
        case PARSER_EXE_MUL_VN:
        {
            pstack.top() *= ((ParserExeMUL_VN*)p)->v; // NOLINT
            p       += sizeof(ParserExeMUL_VN);
            break;
        }
        case PARSER_EXE_DIV_VN:
        {
            pstack.top() = ((ParserExeDIV_VN*)p)->v / pstack.top(); // NOLINT
            p      += sizeof(ParserExeDIV_VN);
            break;
        }
        case PARSER_EXE_ADD_PN:
        {
            int i = ((ParserExeADD_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.top() += d; // NOLINT
            p         += sizeof(ParserExeADD_PN);
            break;
        }
        case PARSER_EXE_SUB_PN:
        {
            int i = ((ParserExeSUB_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.top() = (d - pstack.top()) * (((ParserExeSUB_PN*)p)->sign); // NOLINT
            p         += sizeof(ParserExeSUB_PN);
            break;
        }
        case PARSER_EXE_MUL_PN:
        {
            int i = ((ParserExeMUL_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.top() *= d; // NOLINT
            p         += sizeof(ParserExeMUL_PN);
            break;
        }
        case PARSER_EXE_DIV_PN:
        {
            int i = ((ParserExeDIV_PN*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            if (((ParserExeDIV_PN*)p)->reverse) {
                pstack.top() /= d; // NOLINT
            } else {
                pstack.top() = d / pstack.top(); // NOLINT
            }
            p            += sizeof(ParserExeDIV_PN);
            break;
        }
        case PARSER_EXE_IF:
        {
            double cond = pstack.top(); // NOLINT
            pstack.pop();
            if (cond == 0.0) { // false branch
                p += ((ParserExeIF*)p)->offset;
            }
            p += sizeof(ParserExeIF);
            break;
        }
        case PARSER_EXE_JUMP:
        {
            int offset = ((ParserExeJUMP*)p)->offset;
            p += sizeof(ParserExeJUMP) + offset;
            break;
        }
        default:
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,"parser_exe_eval: unknown node type");
        }
    }
    return pstack.top(); // NOLINT
}

void parser_compile_exe_size (struct parser_node* node, char*& p, std::size_t& exe_size,
                              int& max_stack_size, int& stack_size, Vector<char*>& local_variables);

inline std::size_t
parser_exe_size (struct amrex_parser* parser, int& max_stack_size, int& stack_size)
{
    char* p = nullptr;
    std::size_t exe_size = 0;
    max_stack_size = 0;
    stack_size = 0;
    Vector<char*> local_variables;
    parser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables);
    stack_size -= static_cast<int>(local_variables.size())+1;
    return exe_size+sizeof(ParserExeNull);
}

inline void
parser_compile (struct amrex_parser* parser, char* p)
{
    std::size_t exe_size = 0;
    int max_stack_size = 0;
    int stack_size = 0;
    Vector<char*> local_variables;
    parser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables);
    new(p) ParserExeNull;
}

}

#endif
