#ifndef AMREX_IPARSER_EXE_H_
#define AMREX_IPARSER_EXE_H_
#include <AMReX_Config.H>

#include <AMReX_IParser_Y.H>
#include <AMReX_Vector.H>

#ifndef AMREX_IPARSER_STACK_SIZE
#define AMREX_IPARSER_STACK_SIZE 16
#endif

#define AMREX_IPARSER_LOCAL_IDX0 1000
#define AMREX_IPARSER_GET_DATA(i) ((i)>=1000) ? pstack[(i)-1000] : x[i]

namespace amrex {

// N: node
// P: pointer offset
// V: value (i.e., int literal)

enum iparser_exe_t {
    IPARSER_EXE_NULL = 0,
    IPARSER_EXE_NUMBER,
    IPARSER_EXE_SYMBOL,
    IPARSER_EXE_ADD,
    IPARSER_EXE_SUB,
    IPARSER_EXE_MUL,
    IPARSER_EXE_DIV_F,
    IPARSER_EXE_DIV_B,
    IPARSER_EXE_NEG,
    IPARSER_EXE_F1,
    IPARSER_EXE_F2_F,
    IPARSER_EXE_F2_B,
    IPARSER_EXE_ADD_VP,
    IPARSER_EXE_SUB_VP,
    IPARSER_EXE_MUL_VP,
    IPARSER_EXE_DIV_VP,
    IPARSER_EXE_DIV_PV,
    IPARSER_EXE_ADD_PP,
    IPARSER_EXE_SUB_PP,
    IPARSER_EXE_MUL_PP,
    IPARSER_EXE_DIV_PP,
    IPARSER_EXE_NEG_P,
    IPARSER_EXE_ADD_VN,
    IPARSER_EXE_SUB_VN,
    IPARSER_EXE_MUL_VN,
    IPARSER_EXE_DIV_NV,
    IPARSER_EXE_DIV_VN,
    IPARSER_EXE_ADD_PN,
    IPARSER_EXE_SUB_PN,
    IPARSER_EXE_MUL_PN,
    IPARSER_EXE_DIV_PN,
    IPARSER_EXE_IF,
    IPARSER_EXE_JUMP
};

struct alignas(8) IParserExeNull {
    enum iparser_exe_t type = IPARSER_EXE_NULL;
};

struct alignas(8) IParserExeNumber {
    enum iparser_exe_t type = IPARSER_EXE_NUMBER;
    int v;
};

struct alignas(8) IParserExeSymbol {
    enum iparser_exe_t type = IPARSER_EXE_SYMBOL;
    int i;
};

struct alignas(8) IParserExeADD {
    enum iparser_exe_t type = IPARSER_EXE_ADD;
};

struct alignas(8) IParserExeSUB {
    enum iparser_exe_t type = IPARSER_EXE_SUB;
    int sign;
};

struct alignas(8) IParserExeMUL {
    enum iparser_exe_t type = IPARSER_EXE_MUL;
};

struct alignas(8) IParserExeDIV_F {
    enum iparser_exe_t type = IPARSER_EXE_DIV_F;
};

struct alignas(8) IParserExeDIV_B {
    enum iparser_exe_t type = IPARSER_EXE_DIV_B;
};

struct alignas(8) IParserExeNEG {
    enum iparser_exe_t type = IPARSER_EXE_NEG;
};

struct alignas(8) IParserExeF1 {
    enum iparser_exe_t type = IPARSER_EXE_F1;
    iparser_f1_t ftype;
};

struct alignas(8) IParserExeF2_F {
    enum iparser_exe_t type = IPARSER_EXE_F2_F;
    iparser_f2_t ftype;
};

struct alignas(8) IParserExeF2_B {
    enum iparser_exe_t type = IPARSER_EXE_F2_B;
    iparser_f2_t ftype;
};

struct alignas(8) IParserExeADD_VP {
    enum iparser_exe_t type = IPARSER_EXE_ADD_VP;
    int i;
    int v;
};

struct alignas(8) IParserExeSUB_VP {
    enum iparser_exe_t type = IPARSER_EXE_SUB_VP;
    int i;
    int v;
};

struct alignas(8) IParserExeMUL_VP {
    enum iparser_exe_t type = IPARSER_EXE_MUL_VP;
    int i;
    int v;
};

struct alignas(8) IParserExeDIV_VP {
    enum iparser_exe_t type = IPARSER_EXE_DIV_VP;
    int i;
    int v;
};

struct alignas(8) IParserExeDIV_PV {
    enum iparser_exe_t type = IPARSER_EXE_DIV_PV;
    int i;
    int v;
};

struct alignas(8) IParserExeADD_PP {
    enum iparser_exe_t type = IPARSER_EXE_ADD_PP;
    int i1;
    int i2;
};

struct alignas(8) IParserExeSUB_PP {
    enum iparser_exe_t type = IPARSER_EXE_SUB_PP;
    int i1;
    int i2;
};

struct alignas(8) IParserExeMUL_PP {
    enum iparser_exe_t type = IPARSER_EXE_MUL_PP;
    int i1;
    int i2;
};

struct alignas(8) IParserExeDIV_PP {
    enum iparser_exe_t type = IPARSER_EXE_DIV_PP;
    int i1;
    int i2;
};

struct alignas(8) IParserExeNEG_P {
    enum iparser_exe_t type = IPARSER_EXE_NEG_P;
    int i;
};

struct alignas(8) IParserExeADD_VN {
    enum iparser_exe_t type = IPARSER_EXE_ADD_VN;
    int v;
};

struct alignas(8) IParserExeSUB_VN {
    enum iparser_exe_t type = IPARSER_EXE_SUB_VN;
    int v;
};

struct alignas(8) IParserExeMUL_VN {
    enum iparser_exe_t type = IPARSER_EXE_MUL_VN;
    int v;
};

struct alignas(8) IParserExeDIV_VN {
    enum iparser_exe_t type = IPARSER_EXE_DIV_VN;
    int v;
};

struct alignas(8) IParserExeDIV_NV {
    enum iparser_exe_t type = IPARSER_EXE_DIV_NV;
    int v;
};

struct alignas(8) IParserExeADD_PN {
    enum iparser_exe_t type = IPARSER_EXE_ADD_PN;
    int i;
};

struct alignas(8) IParserExeSUB_PN {
    enum iparser_exe_t type = IPARSER_EXE_SUB_PN;
    int i;
    int sign;
};

struct alignas(8) IParserExeMUL_PN {
    enum iparser_exe_t type = IPARSER_EXE_MUL_PN;
    int i;
};

struct alignas(8) IParserExeDIV_PN {
    enum iparser_exe_t type = IPARSER_EXE_DIV_PN;
    int i;
    bool reverse;
};

struct alignas(8) IParserExeIF {
    enum iparser_exe_t type = IPARSER_EXE_IF;
    int offset;
};

struct alignas(8) IParserExeJUMP {
    enum iparser_exe_t type = IPARSER_EXE_JUMP;
    int offset;
};

template <int N>
struct IParserStack
{
    int m_data[N];
    int m_size = 0;
    constexpr void push (int v) { m_data[m_size++] = v; }
    constexpr void pop () { --m_size; }
    [[nodiscard]] constexpr int const& top () const { return m_data[m_size-1]; }
    [[nodiscard]] constexpr int      & top ()       { return m_data[m_size-1]; }
    [[nodiscard]] constexpr int operator[] (int i) const { return m_data[i]; }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int iparser_exe_eval (const char* p, int const* x)
{
    IParserStack<AMREX_IPARSER_STACK_SIZE> pstack;
    while (*((iparser_exe_t*)p) != IPARSER_EXE_NULL) {
        switch (*((iparser_exe_t*)p))
        {
        case IPARSER_EXE_NUMBER:
        {
            pstack.push(((IParserExeNumber*)p)->v);
            p   += sizeof(IParserExeNumber);
            break;
        }
        case IPARSER_EXE_SYMBOL:
        {
            int i = ((IParserExeSymbol*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(d);
            p    += sizeof(IParserExeSymbol);
            break;
        }
        case IPARSER_EXE_ADD:
        {
            int b = pstack.top();
            pstack.pop();
            pstack.top() += b;
            p += sizeof(IParserExeADD);
            break;
        }
        case IPARSER_EXE_SUB:
        {
            int b = pstack.top();
            pstack.pop();
            pstack.top() = (pstack.top() - b) * (((IParserExeSUB*)p)->sign);
            p += sizeof(IParserExeSUB);
            break;
        }
        case IPARSER_EXE_MUL:
        {
            int b = pstack.top();
            pstack.pop();
            pstack.top() *= b;
            p += sizeof(IParserExeMUL);
            break;
        }
        case IPARSER_EXE_DIV_F:
        {
            int v = pstack.top();
            pstack.pop();
            pstack.top() /= v;
            p += sizeof(IParserExeDIV_F);
            break;
        }
        case IPARSER_EXE_DIV_B:
        {
            int v = pstack.top();
            pstack.pop();
            pstack.top() = v / pstack.top();
            p += sizeof(IParserExeDIV_B);
            break;
        }
        case IPARSER_EXE_NEG:
        {
            pstack.top() = -pstack.top();
            p += sizeof(IParserExeNEG);
            break;
        }
        case IPARSER_EXE_F1:
        {
            pstack.top() = iparser_call_f1(((IParserExeF1*)p)->ftype, pstack.top());
            p += sizeof(IParserExeF1);
            break;
        }
        case IPARSER_EXE_F2_F:
        {
            int v = pstack.top();
            pstack.pop();
            pstack.top() = iparser_call_f2(((IParserExeF2_F*)p)->ftype, pstack.top(), v);
            p += sizeof(IParserExeF2_F);
            break;
        }
        case IPARSER_EXE_F2_B:
        {
            int v = pstack.top();
            pstack.pop();
            pstack.top() = iparser_call_f2(((IParserExeF2_B*)p)->ftype, v, pstack.top());
            p += sizeof(IParserExeF2_B);
            break;
        }
        case IPARSER_EXE_ADD_VP:
        {
            int i = ((IParserExeADD_VP*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(((IParserExeADD_VP*)p)->v + d);
            p   += sizeof(IParserExeADD_VP);
            break;
        }
        case IPARSER_EXE_SUB_VP:
        {
            int i = ((IParserExeSUB_VP*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(((IParserExeSUB_VP*)p)->v - d);
            p   += sizeof(IParserExeSUB_VP);
            break;
        }
        case IPARSER_EXE_MUL_VP:
        {
            int i = ((IParserExeMUL_VP*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(((IParserExeMUL_VP*)p)->v * d);
            p   += sizeof(IParserExeMUL_VP);
            break;
        }
        case IPARSER_EXE_DIV_VP:
        {
            int i = ((IParserExeDIV_VP*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(((IParserExeDIV_VP*)p)->v / d);
            p   += sizeof(IParserExeDIV_VP);
            break;
        }
        case IPARSER_EXE_DIV_PV:
        {
            int i = ((IParserExeDIV_PV*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(d / ((IParserExeDIV_PV*)p)->v);
            p       += sizeof(IParserExeDIV_PV);
            break;
        }
        case IPARSER_EXE_ADD_PP:
        {
            int i = ((IParserExeADD_PP*)p)->i1;
            int d1 = AMREX_IPARSER_GET_DATA(i);
            i     = ((IParserExeADD_PP*)p)->i2;
            int d2 = AMREX_IPARSER_GET_DATA(i);
            pstack.push(d1+d2);
            p     += sizeof(IParserExeADD_PP);
            break;
        }
        case IPARSER_EXE_SUB_PP:
        {
            int i = ((IParserExeSUB_PP*)p)->i1;
            int d1 = AMREX_IPARSER_GET_DATA(i);
            i     = ((IParserExeSUB_PP*)p)->i2;
            int d2 = AMREX_IPARSER_GET_DATA(i);
            pstack.push(d1-d2);
            p     += sizeof(IParserExeSUB_PP);
            break;
        }
        case IPARSER_EXE_MUL_PP:
        {
            int i = ((IParserExeMUL_PP*)p)->i1;
            int d1 = AMREX_IPARSER_GET_DATA(i);
            i     = ((IParserExeMUL_PP*)p)->i2;
            int d2 = AMREX_IPARSER_GET_DATA(i);
            pstack.push(d1*d2);
            p     += sizeof(IParserExeMUL_PP);
            break;
        }
        case IPARSER_EXE_DIV_PP:
        {
            int i = ((IParserExeDIV_PP*)p)->i1;
            int d1 = AMREX_IPARSER_GET_DATA(i);
            i     = ((IParserExeDIV_PP*)p)->i2;
            int d2 = AMREX_IPARSER_GET_DATA(i);
            pstack.push(d1/d2);
            p      += sizeof(IParserExeDIV_PP);
            break;
        }
        case IPARSER_EXE_NEG_P:
        {
            int i = ((IParserExeNEG_P*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.push(-d);
            p     += sizeof(IParserExeNEG_P);
            break;
        }
        case IPARSER_EXE_ADD_VN:
        {
            pstack.top() += ((IParserExeADD_VN*)p)->v;
            p       += sizeof(IParserExeADD_VN);
            break;
        }
        case IPARSER_EXE_SUB_VN:
        {
            pstack.top() = ((IParserExeSUB_VN*)p)->v - pstack.top();
            p      += sizeof(IParserExeSUB_VN);
            break;
        }
        case IPARSER_EXE_MUL_VN:
        {
            pstack.top() *= ((IParserExeMUL_VN*)p)->v;
            p       += sizeof(IParserExeMUL_VN);
            break;
        }
        case IPARSER_EXE_DIV_VN:
        {
            pstack.top() = ((IParserExeDIV_VN*)p)->v / pstack.top();
            p      += sizeof(IParserExeDIV_VN);
            break;
        }
        case IPARSER_EXE_DIV_NV:
        {
            pstack.top() /= ((IParserExeDIV_NV*)p)->v;
            p       += sizeof(IParserExeDIV_NV);
            break;
        }
        case IPARSER_EXE_ADD_PN:
        {
            int i = ((IParserExeADD_PN*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.top() += d;
            p         += sizeof(IParserExeADD_PN);
            break;
        }
        case IPARSER_EXE_SUB_PN:
        {
            int i = ((IParserExeSUB_PN*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.top() = (d - pstack.top()) * (((IParserExeSUB_PN*)p)->sign);
            p         += sizeof(IParserExeSUB_PN);
            break;
        }
        case IPARSER_EXE_MUL_PN:
        {
            int i = ((IParserExeMUL_PN*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            pstack.top() *= d;
            p         += sizeof(IParserExeMUL_PN);
            break;
        }
        case IPARSER_EXE_DIV_PN:
        {
            int i = ((IParserExeDIV_PN*)p)->i;
            int d = AMREX_IPARSER_GET_DATA(i);
            if (((IParserExeDIV_PN*)p)->reverse) {
                pstack.top() /= d;
            } else {
                pstack.top() = d / pstack.top();
            }
            p            += sizeof(IParserExeDIV_PN);
            break;
        }
        case IPARSER_EXE_IF:
        {
            int cond = pstack.top();
            pstack.pop();
            if (cond == 0.0) { // false branch
                p += ((IParserExeIF*)p)->offset;
            }
            p += sizeof(IParserExeIF);
            break;
        }
        case IPARSER_EXE_JUMP:
        {
            int offset = ((IParserExeJUMP*)p)->offset;
            p    += sizeof(IParserExeJUMP) + offset;
            break;
        }
        default:
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,"parser_exe_eval: unknown node type");
        }
    }
    return pstack.top();
}

void iparser_compile_exe_size (struct iparser_node* node, char*& p, std::size_t& exe_size,
                               int& max_stack_size, int& stack_size, Vector<char*>& local_variables);

inline std::size_t
iparser_exe_size (struct amrex_iparser* parser, int& max_stack_size, int& stack_size)
{
    char* p = nullptr;
    std::size_t exe_size = 0;
    max_stack_size = 0;
    stack_size = 0;
    Vector<char*> local_variables;
    iparser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables);
    stack_size -= static_cast<int>(local_variables.size())+1;
    return exe_size+sizeof(IParserExeNull);
}

inline void
iparser_compile (struct amrex_iparser* parser, char* p)
{
    std::size_t exe_size = 0;
    int max_stack_size = 0;
    int stack_size = 0;
    Vector<char*> local_variables;
    iparser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables);
    new(p) IParserExeNull;
}

}

#endif
