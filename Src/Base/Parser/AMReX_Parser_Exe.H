#ifndef AMREX_PARSER_EXE_H_
#define AMREX_PARSER_EXE_H_
#include <AMReX_Config.H>

#include <AMReX_Parser_Y.H>
#include <AMReX_Stack.H>
#include <AMReX_Vector.H>

#include <limits>
#include <utility>

#ifndef AMREX_PARSER_STACK_SIZE
#define AMREX_PARSER_STACK_SIZE 16
#endif

#define AMREX_PARSER_LOCAL_IDX0 1000
#define AMREX_PARSER_GET_DATA(i) ((i)<1000) ? x[i] : pstack[(i)-1000]

namespace amrex {

// N: node
// P: pointer offset
// V: value (i.e., double literal)

/*
 TO: Ternary operator, a op0 b op1 c
     L: (a op0 b) op1 c
     R: a op0 (b op1 c)
     op: 0:+, 1:-, 2:*, 3:/
     For example, TO_L03 means (a + b) / c
     Due to soring, some types are not needed.
*/

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
    PARSER_EXE_SQUARE_P,
    PARSER_EXE_CUBIC_P,
    PARSER_EXE_POWI,
    PARSER_EXE_HYPOT2_P,
    PARSER_EXE_HYPOT2,
    PARSER_EXE_HYPOT3_P,
    PARSER_EXE_HYPOT3,
    // PARSER_EXE_TO_L00, similar to R00
    PARSER_EXE_TO_L01,
    // PARSER_EXE_TO_L02, similar to R20
    PARSER_EXE_TO_L03,
    // PARSER_EXE_TO_L10, similar to R01
    PARSER_EXE_TO_L11,
    // PARSER_EXE_TO_L12, similar to R21
    PARSER_EXE_TO_L13,
    // PARSER_EXE_TO_L20, similar to R02
    PARSER_EXE_TO_L21,
    // PARSER_EXE_TO_L22, similar to R22
    PARSER_EXE_TO_L23,
    // PARSER_EXE_TO_L30, similar to R03
    PARSER_EXE_TO_L31,
    // PARSER_EXE_TO_L32, similar to R23
    // PARSER_EXE_TO_L33, (a/b)/c is similar to a/(b*c)
    PARSER_EXE_TO_R00,
    PARSER_EXE_TO_R01,
    PARSER_EXE_TO_R02,
    PARSER_EXE_TO_R03,
    PARSER_EXE_TO_R10,
    // PARSER_EXE_TO_R11, a-(b-c) is similar to a+(c-b)
    PARSER_EXE_TO_R12,
    PARSER_EXE_TO_R13,
    PARSER_EXE_TO_R20,
    PARSER_EXE_TO_R21,
    PARSER_EXE_TO_R22,
    PARSER_EXE_TO_R23,
    PARSER_EXE_TO_R30,
    PARSER_EXE_TO_R31,
    PARSER_EXE_TO_R32,
    // PARSER_EXE_TO_R33, a/(b/c) is similar to (x*c)/y
    PARSER_EXE_IF,
    PARSER_EXE_JUMP,
    PARSER_EXE_USER_FN
};

struct alignas(double) ParserExeNull {
    enum parser_exe_t type = PARSER_EXE_NULL;
};

struct alignas(double) ParserExeNumber {
    enum parser_exe_t type = PARSER_EXE_NUMBER;
    double v;
};

struct alignas(double) ParserExeSymbol {
    enum parser_exe_t type = PARSER_EXE_SYMBOL;
    int i;
};

struct alignas(double) ParserExeADD {
    enum parser_exe_t type = PARSER_EXE_ADD;
};

struct alignas(double) ParserExeSUB_F {
    enum parser_exe_t type = PARSER_EXE_SUB_F;
};

struct alignas(double) ParserExeSUB_B {
    enum parser_exe_t type = PARSER_EXE_SUB_B;
};

struct alignas(double) ParserExeMUL {
    enum parser_exe_t type = PARSER_EXE_MUL;
};

struct alignas(double) ParserExeDIV_F {
    enum parser_exe_t type = PARSER_EXE_DIV_F;
};

struct alignas(double) ParserExeDIV_B {
    enum parser_exe_t type = PARSER_EXE_DIV_B;
};

struct alignas(double) ParserExeF1 {
    enum parser_exe_t type = PARSER_EXE_F1;
    parser_f1_t ftype;
};

struct alignas(double) ParserExeF2_F {
    enum parser_exe_t type = PARSER_EXE_F2_F;
    parser_f2_t ftype;
};

struct alignas(double) ParserExeF2_B {
    enum parser_exe_t type = PARSER_EXE_F2_B;
    parser_f2_t ftype;
};

struct alignas(double) ParserExeADD_VP {
    enum parser_exe_t type = PARSER_EXE_ADD_VP;
    int i;
    double v;
};

struct alignas(double) ParserExeSUB_VP {
    enum parser_exe_t type = PARSER_EXE_SUB_VP;
    int i;
    double v;
};

struct alignas(double) ParserExeMUL_VP {
    enum parser_exe_t type = PARSER_EXE_MUL_VP;
    int i;
    double v;
};

struct alignas(double) ParserExeDIV_VP {
    enum parser_exe_t type = PARSER_EXE_DIV_VP;
    int i;
    double v;
};

struct alignas(double) ParserExeADD_PP {
    enum parser_exe_t type = PARSER_EXE_ADD_PP;
    int i1;
    int i2;
};

struct alignas(double) ParserExeSUB_PP {
    enum parser_exe_t type = PARSER_EXE_SUB_PP;
    int i1;
    int i2;
};

struct alignas(double) ParserExeMUL_PP {
    enum parser_exe_t type = PARSER_EXE_MUL_PP;
    int i1;
    int i2;
};

struct alignas(double) ParserExeDIV_PP {
    enum parser_exe_t type = PARSER_EXE_DIV_PP;
    int i1;
    int i2;
};

struct alignas(double) ParserExeADD_VN {
    enum parser_exe_t type = PARSER_EXE_ADD_VN;
    double v;
};

struct alignas(double) ParserExeSUB_VN {
    enum parser_exe_t type = PARSER_EXE_SUB_VN;
    double v;
};

struct alignas(double) ParserExeMUL_VN {
    enum parser_exe_t type = PARSER_EXE_MUL_VN;
    double v;
};

struct alignas(double) ParserExeDIV_VN {
    enum parser_exe_t type = PARSER_EXE_DIV_VN;
    double v;
};

struct alignas(double) ParserExeADD_PN {
    enum parser_exe_t type = PARSER_EXE_ADD_PN;
    int i;
};

struct alignas(double) ParserExeSUB_PN {
    enum parser_exe_t type = PARSER_EXE_SUB_PN;
    int i;
    double sign;
};

struct alignas(double) ParserExeMUL_PN {
    enum parser_exe_t type = PARSER_EXE_MUL_PN;
    int i;
};

struct alignas(double) ParserExeDIV_PN {
    enum parser_exe_t type = PARSER_EXE_DIV_PN;
    int i;
    bool reverse;
};

struct alignas(double) ParserExeSquare {
    enum parser_exe_t type = PARSER_EXE_SQUARE;
};

struct alignas(double) ParserExeSquare_P {
    enum parser_exe_t type = PARSER_EXE_SQUARE_P;
    int i;
};

struct alignas(double) ParserExeCubic_P {
    enum parser_exe_t type = PARSER_EXE_CUBIC_P;
    int i;
};

struct alignas(double) ParserExePOWI {
    enum parser_exe_t type = PARSER_EXE_POWI;
    int i;
};

struct alignas(double) ParserExeHypot2_P {
    enum parser_exe_t type = PARSER_EXE_HYPOT2_P;
    int i1;
    int i2;
};

struct alignas(double) ParserExeHypot2 {
    enum parser_exe_t type = PARSER_EXE_HYPOT2;
};

struct alignas(double) ParserExeHypot3_P {
    enum parser_exe_t type = PARSER_EXE_HYPOT3_P;
    int i1;
    int i2;
    int i3;
};

struct alignas(double) ParserExeHypot3 {
    enum parser_exe_t type = PARSER_EXE_HYPOT3;
};

struct alignas(double) ParserExeTO_L01 {
    enum parser_exe_t type = PARSER_EXE_TO_L01;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_L03 {
    enum parser_exe_t type = PARSER_EXE_TO_L03;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_L11 {
    enum parser_exe_t type = PARSER_EXE_TO_L11;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_L13 {
    enum parser_exe_t type = PARSER_EXE_TO_L13;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_L21 {
    enum parser_exe_t type = PARSER_EXE_TO_L21;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_L23 {
    enum parser_exe_t type = PARSER_EXE_TO_L23;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_L31 {
    enum parser_exe_t type = PARSER_EXE_TO_L31;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R00 {
    enum parser_exe_t type = PARSER_EXE_TO_R00;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R01 {
    enum parser_exe_t type = PARSER_EXE_TO_R01;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R02 {
    enum parser_exe_t type = PARSER_EXE_TO_R02;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R03 {
    enum parser_exe_t type = PARSER_EXE_TO_R03;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R10 {
    enum parser_exe_t type = PARSER_EXE_TO_R10;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R12 {
    enum parser_exe_t type = PARSER_EXE_TO_R12;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R13 {
    enum parser_exe_t type = PARSER_EXE_TO_R13;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R20 {
    enum parser_exe_t type = PARSER_EXE_TO_R20;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R21 {
    enum parser_exe_t type = PARSER_EXE_TO_R21;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R22 {
    enum parser_exe_t type = PARSER_EXE_TO_R22;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R23 {
    enum parser_exe_t type = PARSER_EXE_TO_R23;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R30 {
    enum parser_exe_t type = PARSER_EXE_TO_R30;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R31 {
    enum parser_exe_t type = PARSER_EXE_TO_R31;
    int ia, ib, ic;
};

struct alignas(double) ParserExeTO_R32 {
    enum parser_exe_t type = PARSER_EXE_TO_R32;
    int ia, ib, ic;
};

struct alignas(double) ParserExeIF {
    enum parser_exe_t type = PARSER_EXE_IF;
    int offset;
};

struct alignas(double) ParserExeJUMP {
    enum parser_exe_t type = PARSER_EXE_JUMP;
    int offset;
};

struct alignas(double) ParserExeUserFn {
    enum parser_exe_t type = PARSER_EXE_USER_FN;
    short argc;
    void* fh;
    void* fd;
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
        case PARSER_EXE_SQUARE_P:
        {
            int i = ((ParserExeSquare_P*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(d*d);
            p += sizeof(ParserExeSquare_P);
            break;
        }
        case PARSER_EXE_CUBIC_P:
        {
            int i = ((ParserExeCubic_P*)p)->i;
            double d = AMREX_PARSER_GET_DATA(i);
            pstack.push(d*d*d);
            p += sizeof(ParserExeCubic_P);
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
        case PARSER_EXE_HYPOT2_P:
        {
            int i1 = ((ParserExeHypot2_P*)p)->i1;
            int i2 = ((ParserExeHypot2_P*)p)->i2;
            double d1 = AMREX_PARSER_GET_DATA(i1);
            double d2 = AMREX_PARSER_GET_DATA(i2);
            pstack.push(std::sqrt(d1*d1+d2*d2));
            p += sizeof(ParserExeHypot2_P);
            break;
        }
        case PARSER_EXE_HYPOT3_P:
        {
            int i1 = ((ParserExeHypot3_P*)p)->i1;
            int i2 = ((ParserExeHypot3_P*)p)->i2;
            int i3 = ((ParserExeHypot3_P*)p)->i3;
            double d1 = AMREX_PARSER_GET_DATA(i1);
            double d2 = AMREX_PARSER_GET_DATA(i2);
            double d3 = AMREX_PARSER_GET_DATA(i3);
            pstack.push(std::sqrt(d1*d1+d2*d2+d3*d3));
            p += sizeof(ParserExeHypot3_P);
            break;
        }
        case PARSER_EXE_HYPOT2:
        {
            double v1 = pstack.top(); // NOLINT;
            pstack.pop();
            double v2 = pstack.top(); // NOLINT;
            pstack.top() = std::sqrt(v1*v1+v2*v2);
            p += sizeof(ParserExeHypot2);
            break;
        }
        case PARSER_EXE_HYPOT3:
        {
            double v1 = pstack.top(); // NOLINT;
            pstack.pop();
            double v2 = pstack.top(); // NOLINT;
            pstack.pop();
            double v3 = pstack.top(); // NOLINT;
            pstack.top() = std::sqrt(v1*v1+v2*v2+v3*v3);
            p += sizeof(ParserExeHypot3);
            break;
        }
        case PARSER_EXE_TO_L01:
        {
            auto* q  = (ParserExeTO_L01*)p;
            p += sizeof(ParserExeTO_L01);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a+b)-c);
            break;
        }
        case PARSER_EXE_TO_L03:
        {
            auto* q  = (ParserExeTO_L03*)p;
            p += sizeof(ParserExeTO_L03);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a+b)/c);
            break;
        }
        case PARSER_EXE_TO_L11:
        {
            auto* q  = (ParserExeTO_L11*)p;
            p += sizeof(ParserExeTO_L11);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a-b)-c);
            break;
        }
        case PARSER_EXE_TO_L13:
        {
            auto* q  = (ParserExeTO_L13*)p;
            p += sizeof(ParserExeTO_L13);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a-b)/c);
            break;
        }
        case PARSER_EXE_TO_L21:
        {
            auto* q  = (ParserExeTO_L21*)p;
            p += sizeof(ParserExeTO_L21);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a*b)-c);
            break;
        }
        case PARSER_EXE_TO_L23:
        {
            auto* q  = (ParserExeTO_L23*)p;
            p += sizeof(ParserExeTO_L23);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a*b)/c);
            break;
        }
        case PARSER_EXE_TO_L31:
        {
            auto* q  = (ParserExeTO_L31*)p;
            p += sizeof(ParserExeTO_L31);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push((a/b)-c);
            break;
        }
        case PARSER_EXE_TO_R00:
        {
            auto* q  = (ParserExeTO_R00*)p;
            p += sizeof(ParserExeTO_R00);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a+(b+c));
            break;
        }
        case PARSER_EXE_TO_R01:
        {
            auto* q  = (ParserExeTO_R01*)p;
            p += sizeof(ParserExeTO_R01);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a+(b-c));
            break;
        }
        case PARSER_EXE_TO_R02:
        {
            auto* q  = (ParserExeTO_R02*)p;
            p += sizeof(ParserExeTO_R02);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a+(b*c));
            break;
        }
        case PARSER_EXE_TO_R03:
        {
            auto* q  = (ParserExeTO_R03*)p;
            p += sizeof(ParserExeTO_R03);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a+(b/c));
            break;
        }
        case PARSER_EXE_TO_R10:
        {
            auto* q  = (ParserExeTO_R10*)p;
            p += sizeof(ParserExeTO_R10);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a-(b+c));
            break;
        }
        case PARSER_EXE_TO_R12:
        {
            auto* q  = (ParserExeTO_R12*)p;
            p += sizeof(ParserExeTO_R12);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a-(b*c));
            break;
        }
        case PARSER_EXE_TO_R13:
        {
            auto* q  = (ParserExeTO_R13*)p;
            p += sizeof(ParserExeTO_R13);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a-(b/c));
            break;
        }
        case PARSER_EXE_TO_R20:
        {
            auto* q  = (ParserExeTO_R20*)p;
            p += sizeof(ParserExeTO_R20);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a*(b+c));
            break;
        }
        case PARSER_EXE_TO_R21:
        {
            auto* q  = (ParserExeTO_R21*)p;
            p += sizeof(ParserExeTO_R21);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a*(b-c));
            break;
        }
        case PARSER_EXE_TO_R22:
        {
            auto* q  = (ParserExeTO_R22*)p;
            p += sizeof(ParserExeTO_R22);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a*(b*c));
            break;
        }
        case PARSER_EXE_TO_R23:
        {
            auto* q  = (ParserExeTO_R23*)p;
            p += sizeof(ParserExeTO_R23);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a*(b/c));
            break;
        }
        case PARSER_EXE_TO_R30:
        {
            auto* q  = (ParserExeTO_R30*)p;
            p += sizeof(ParserExeTO_R30);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a/(b+c));
            break;
        }
        case PARSER_EXE_TO_R31:
        {
            auto* q  = (ParserExeTO_R31*)p;
            p += sizeof(ParserExeTO_R31);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a/(b-c));
            break;
        }
        case PARSER_EXE_TO_R32:
        {
            auto* q  = (ParserExeTO_R32*)p;
            p += sizeof(ParserExeTO_R32);
            double a = AMREX_PARSER_GET_DATA(q->ia);
            double b = AMREX_PARSER_GET_DATA(q->ib);
            double c = AMREX_PARSER_GET_DATA(q->ic);
            pstack.push(a/(b*c));
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
        case PARSER_EXE_USER_FN:
        {
            short argc = ((ParserExeUserFn*)p)->argc;
            if (argc == 1) {
                AMREX_IF_ON_HOST  (( pstack.top() = ParserUserFn1(((ParserExeUserFn*)p)->fh)(pstack.top()); ))
#if !defined(AMREX_USE_SYCL)
                AMREX_IF_ON_DEVICE(( pstack.top() = ParserUserFn1(((ParserExeUserFn*)p)->fd)(pstack.top()); ))
#endif
            } else if (argc == 2) {
                double a1 = pstack.top(); // NOLINT
                pstack.pop();
                AMREX_IF_ON_HOST  (( pstack.top() = ParserUserFn2(((ParserExeUserFn*)p)->fh)(pstack.top(),a1); ))
#if !defined(AMREX_USE_SYCL)
                AMREX_IF_ON_DEVICE(( pstack.top() = ParserUserFn2(((ParserExeUserFn*)p)->fd)(pstack.top(),a1); ))
#endif
            } else if (argc == 3) {
                double a2 = pstack.top(); // NOLINT
                pstack.pop();
                double a1 = pstack.top(); // NOLINT
                pstack.pop();
                AMREX_IF_ON_HOST  (( pstack.top() = ParserUserFn3(((ParserExeUserFn*)p)->fh)(pstack.top(),a1,a2); ))
#if !defined(AMREX_USE_SYCL)
                AMREX_IF_ON_DEVICE(( pstack.top() = ParserUserFn3(((ParserExeUserFn*)p)->fd)(pstack.top(),a1,a2); ))
#endif
            } else if (argc == 4) {
                double a3 = pstack.top(); // NOLINT
                pstack.pop();
                double a2 = pstack.top(); // NOLINT
                pstack.pop();
                double a1 = pstack.top(); // NOLINT
                pstack.pop();
                AMREX_IF_ON_HOST  (( pstack.top() = ParserUserFn4(((ParserExeUserFn*)p)->fh)(pstack.top(),a1,a2,a3); ))
#if !defined(AMREX_USE_SYCL)
                AMREX_IF_ON_DEVICE(( pstack.top() = ParserUserFn4(((ParserExeUserFn*)p)->fd)(pstack.top(),a1,a2,a3); ))
#endif
            }
            p += sizeof(ParserExeUserFn);
            break;
        }
        default:
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,"parser_exe_eval: unknown node type");
        }
    }
    return pstack.top(); // NOLINT
}

void parser_compile_exe_size (struct parser_node* node, char*& p, std::size_t& exe_size,
                              int& max_stack_size, int& stack_size,
                              Vector<char const*>& local_variables,
                              std::map<std::string,std::pair<void*,void*>> const& ufs);

inline std::size_t
parser_exe_size (struct amrex_parser* parser, int& max_stack_size, int& stack_size)
{
    parser_ast_sort(parser->ast);
    char* p = nullptr;
    std::size_t exe_size = 0;
    max_stack_size = 0;
    stack_size = 0;
    Vector<char const*> local_variables;
    std::map<std::string,std::pair<void*,void*>> ufs;
    parser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables, ufs);
    stack_size -= static_cast<int>(local_variables.size())+1;
    return exe_size+sizeof(ParserExeNull);
}

inline Vector<char const*>
parser_compile (struct amrex_parser* parser,
                std::map<std::string,std::pair<void*,void*>> const& ufs, char* p)
{
    std::size_t exe_size = 0;
    int max_stack_size = 0;
    int stack_size = 0;
    Vector<char const*> local_variables;
    parser_compile_exe_size(parser->ast, p, exe_size, max_stack_size, stack_size, local_variables, ufs);
    new(p) ParserExeNull;
    return local_variables;
}

void parser_exe_print(char const* parser, Vector<std::string> const& vars,
                      Vector<char const*> const& locals);

}

#endif
