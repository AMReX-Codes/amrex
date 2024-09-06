#ifndef AMREX_PARSER_EXE_H_
#define AMREX_PARSER_EXE_H_
#include <AMReX_Config.H>

#include <AMReX_Parser_Y.H>
#include <AMReX_Stack.H>
#include <AMReX_Vector.H>

#include <limits>

#ifndef AMREX_PARSER_STACK_SIZE
#define AMREX_PARSER_STACK_SIZE 16
#endif

#define AMREX_PARSER_LOCAL_IDX0 1000
#define AMREX_PARSER_GET_DATA(i) ((i)<1000) ? x[i] : pstack[(i)-1000]

namespace amrex {

// N: node
// P: pointer offset
// V: value (i.e., double literal)

enum parser_exe_t {
    PARSER_EXE_NULL = 0,
    PARSER_EXE_NUMBER,
    PARSER_EXE_SYMBOL,
    PARSER_EXE_ADD,
    PARSER_EXE_SUB_F,
    PARSER_EXE_SUB_B,
    PARSER_EXE_MUL,
    PARSER_EXE_DIV_F,
    PARSER_EXE_DIV_B,
    PARSER_EXE_F1,
    PARSER_EXE_F2_F,
    PARSER_EXE_F2_B,
    PARSER_EXE_ADD_VP,
    PARSER_EXE_SUB_VP,
    PARSER_EXE_MUL_VP,
    PARSER_EXE_DIV_VP,
    PARSER_EXE_ADD_PP,
    PARSER_EXE_SUB_PP,
    PARSER_EXE_MUL_PP,
    PARSER_EXE_DIV_PP,
    PARSER_EXE_ADD_VN,
    PARSER_EXE_SUB_VN,
    PARSER_EXE_MUL_VN,
    PARSER_EXE_DIV_VN,
    PARSER_EXE_ADD_PN,
    PARSER_EXE_SUB_PN,
    PARSER_EXE_MUL_PN,
    PARSER_EXE_DIV_PN,
    PARSER_EXE_SQUARE,
    PARSER_EXE_POWI,
    PARSER_EXE_IF,
    PARSER_EXE_JUMP
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

struct alignas(8) ParserExeSUB_F {
    enum parser_exe_t type = PARSER_EXE_SUB_F;
};

struct alignas(8) ParserExeSUB_B {
    enum parser_exe_t type = PARSER_EXE_SUB_B;
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

struct alignas(8) ParserExeSquare {
    enum parser_exe_t type = PARSER_EXE_SQUARE;
};

struct alignas(8) ParserExePOWI {
    enum parser_exe_t type = PARSER_EXE_POWI;
    int i;
};

struct alignas(8) ParserExeIF {
    enum parser_exe_t type = PARSER_EXE_IF;
    int offset;
};

struct alignas(8) ParserExeJUMP {
    enum parser_exe_t type = PARSER_EXE_JUMP;
    int offset;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double parser_exe_eval (const char* p, double const* x)
{
    if (p == nullptr) { return std::numeric_limits<double>::max(); }

    Stack<double, AMREX_PARSER_STACK_SIZE> pstack;
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
        case PARSER_EXE_SUB_F:
        {
            double b = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() -= b; // NOLINT
            p += sizeof(ParserExeSUB_F);
            break;
        }
        case PARSER_EXE_SUB_B:
        {
            double b = pstack.top(); // NOLINT
            pstack.pop();
            pstack.top() = b - pstack.top(); // NOLINT
            p += sizeof(ParserExeSUB_B);
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
        case PARSER_EXE_SQUARE:
        {
            double& d = pstack.top();
            d *= d;
            p += sizeof(ParserExeSquare);
            break;
        }
        case PARSER_EXE_POWI:
        {
            double& d = pstack.top();
            int n = ((ParserExePOWI*)p)->i;
            if (n != 0) {
                if (n < 0) {
                    d = 1.0/d;
                    n = -n;
                }
                double y = 1.0;
                while (n > 1) {
                    if (n % 2 == 0) {
                        d *= d;
                        n = n/2;
                    } else {
                        y *= d;
                        d *= d;
                        n = (n-1)/2;
                    }
                }
                d *= y;
            } else {
                d = 1.0;
            }
            p += sizeof(ParserExePOWI);
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
                              int& max_stack_size, int& stack_size, Vector<char const*>& local_variables);

inline std::size_t
parser_exe_size (struct amrex_parser* parser, int& max_stack_size, int& stack_size)
{
    parser_ast_sort(parser->ast);
    char* p = nullptr;
    std::size_t exe_size = 0;
    max_stack_size = 0;
    stack_size = 0;
    Vector<char const*> local_variables;
    parser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables);
    stack_size -= static_cast<int>(local_variables.size())+1;
    return exe_size+sizeof(ParserExeNull);
}

inline Vector<char const*>
parser_compile (struct amrex_parser* parser, char* p)
{
    std::size_t exe_size = 0;
    int max_stack_size = 0;
    int stack_size = 0;
    Vector<char const*> local_variables;
    parser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables);
    new(p) ParserExeNull;
    return local_variables;
}

void parser_exe_print(char const* parser, Vector<std::string> const& vars,
                      Vector<char const*> const& locals);

}

#endif
