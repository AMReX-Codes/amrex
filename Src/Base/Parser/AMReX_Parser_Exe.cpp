#include <AMReX_Parser_Exe.H>
#include <limits>
#include <utility>

namespace amrex {

#if AMREX_PARSER_THREADED_DISPATCH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
double parser_exe_eval_threaded (const char* p, double const* x)
{
    if (p == nullptr) { return std::numeric_limits<double>::max(); }
#define AMREX_PARSER_THREADED_BODY 1
#include "AMReX_Parser_Exe_Body.H"
#undef AMREX_PARSER_THREADED_BODY
#undef AMREX_PARSER_CASE
#undef AMREX_PARSER_NEXT
}
#pragma GCC diagnostic pop
#else
double parser_exe_eval_threaded (const char* p, double const* x)
{
    return parser_exe_eval(p, x);
}
#endif

void
parser_compile_exe_size (struct parser_node* node, char*& p, std::size_t& exe_size,
                         int& max_stack_size, int& stack_size,
                         Vector<char const*>& local_variables,
                         std::map<std::string,std::pair<void*,void*>> const& ufs)
{
    auto parser_symbol_idx = [&] (struct parser_node* snode) -> int
    {
        AMREX_ASSERT(snode->type == PARSER_SYMBOL);
        auto* sym = (struct parser_symbol*)snode;
        auto r = std::find_if(local_variables.rbegin(), local_variables.rend(),
                              [=] (char const* i)
                                  { return std::strcmp(sym->name, i) == 0; });
        if (r != local_variables.rend()) {
            return static_cast<int>(std::distance(r, local_variables.rend())) - 1
                + AMREX_PARSER_LOCAL_IDX0;
        } else {
            if (sym->ip < 0) {
                throw std::runtime_error(std::string("Unknown variable ") + sym->name);
            }
            return sym->ip;
        }
    };

    // In parser_exe_eval, we push to the stack for NUMBER, SYMBOL, VP, PP.
    // In parser_exe_eval, we pop the stack for ADD, SUB, MUL, DIV, F2, and IF.

    // Note that for + and * the nodes have been sorted before this function
    // is called. So we don't need to worry about cases like f(x) + x.

    // Note that there is no PARSER_SUB. a-b is actually a+(-b).

    switch (node->type)
    {
    case PARSER_NUMBER:
    {
        if (p) {
            auto *t = new(p) ParserExeNumber;
            p      += sizeof(ParserExeNumber);
            t->v = ((struct parser_number*)node)->value;
        }
        exe_size += sizeof(ParserExeNumber);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_SYMBOL:
    {
        if (p) {
            auto *t = new(p) ParserExeSymbol;
            p      += sizeof(ParserExeSymbol);
            t->i = parser_symbol_idx(node);
        }
        exe_size += sizeof(ParserExeSymbol);
        ++stack_size;
        max_stack_size = std::max(max_stack_size, stack_size);
        break;
    }
    case PARSER_ADD:
    {
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_SYMBOL)
        { // 3 + x
            if (p) {
                auto *t = new(p) ParserExeADD_VP;
                p      += sizeof(ParserExeADD_VP);
                t->i = parser_symbol_idx(node->r);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeADD_VP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 node->r->r->type == PARSER_SYMBOL &&
                 parser_get_number(node->r->l) == -1.0)
        { // 3 + (-1.0)*x => 3 - x
            if (p) {
                auto *t = new(p) ParserExeSUB_VP;
                p      += sizeof(ParserExeSUB_VP);
                t->i = parser_symbol_idx(node->r->r);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeSUB_VP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 node->r->r->type == PARSER_SYMBOL)
        { // b + a*x => fma(a, x, b)
            if (p) {
                auto *t = new(p) ParserExeFMA_VPV;
                p      += sizeof(ParserExeFMA_VPV);
                t->i = parser_symbol_idx(node->r->r);
                t->a = parser_get_number(node->r->l);
                t->b = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeFMA_VPV);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 parser_get_number(node->r->l) == -1.0)
        { // 3 + (-1)*f(x) => 3 - f(x)
            parser_compile_exe_size(node->r->r, p, exe_size, max_stack_size,
                                    stack_size, local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeSUB_VN;
                p      += sizeof(ParserExeSUB_VN);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeSUB_VN);
        }
        else if (node->l->type == PARSER_NUMBER)
        { // 3 + f(x) => 3 + f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeADD_VN;
                p      += sizeof(ParserExeADD_VN);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeADD_VN);
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 node->r->r->type == PARSER_SYMBOL &&
                 parser_get_number(node->r->l) == -1.0)
        { // x + -y => x - y
            if (p) {
                auto *t = new(p) ParserExeSUB_PP;
                p      += sizeof(ParserExeSUB_PP);
                t->i1 = parser_symbol_idx(node->l);
                t->i2 = parser_symbol_idx(node->r->r);
            }
            exe_size += sizeof(ParserExeSUB_PP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        { // x + y
            if (p) {
                auto *t = new(p) ParserExeADD_PP;
                p      += sizeof(ParserExeADD_PP);
                t->i1 = parser_symbol_idx(node->l);
                t->i2 = parser_symbol_idx(node->r);
            }
            exe_size += sizeof(ParserExeADD_PP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 parser_get_number(node->r->l) == -1.0)
        { // x + (-1)*f(x) => x - f(x)
            parser_compile_exe_size(node->r->r, p, exe_size, max_stack_size,
                                    stack_size, local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeSUB_PN;
                p      += sizeof(ParserExeSUB_PN);
                t->i = parser_symbol_idx(node->l);
                t->reverse = false;
            }
            exe_size += sizeof(ParserExeSUB_PN);
        }
        else if (node->l->type == PARSER_SYMBOL)
        { // x + f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeADD_PN;
                p      += sizeof(ParserExeADD_PN);
                t->i = parser_symbol_idx(node->l);
            }
            exe_size += sizeof(ParserExeADD_PN);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->l->r->type == PARSER_SYMBOL &&
                 parser_get_number(node->l->l) == -1.0)
        { // -x + f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeSUB_PN;
                p      += sizeof(ParserExeSUB_PN);
                t->i = parser_symbol_idx(node->l->r);
                t->reverse = true;
            }
            exe_size += sizeof(ParserExeSUB_PN);
        }
        else if (node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 node->r->r->type == PARSER_SYMBOL &&
                 parser_get_number(node->r->l) == -1.0)
        { // f(x) + (-1)*x => -(x-f(x))
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeSUB_PN;
                p      += sizeof(ParserExeSUB_PN);
                t->i = parser_symbol_idx(node->r->r);
                t->reverse = true;
            }
            exe_size += sizeof(ParserExeSUB_PN);
        }
        else if (node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 parser_get_number(node->r->l) == -1.0)
        { // f(x) + (-1)*g(x) => f(x) - g(x)
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeSUB_B;
                    p += sizeof(ParserExeSUB_B);
                }
                exe_size += sizeof(ParserExeSUB_B);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->r->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeSUB_F;
                    p += sizeof(ParserExeSUB_F);
                }
                exe_size += sizeof(ParserExeSUB_F);
            }
            --stack_size;
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 parser_get_number(node->l->l) == -1.0)
        { // (-1)*f(x) + g(x) => g(x) - f(x)
            int d1 = parser_ast_depth(node->l->r);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->l->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeSUB_F;
                    p += sizeof(ParserExeSUB_F);
                }
                exe_size += sizeof(ParserExeSUB_F);
            } else {
                parser_compile_exe_size(node->l->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeSUB_B;
                    p += sizeof(ParserExeSUB_B);
                }
                exe_size += sizeof(ParserExeSUB_B);
            }
            --stack_size;
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
            }
            if (p) {
                new(p)      ParserExeADD;
                p += sizeof(ParserExeADD);
            }
            exe_size += sizeof(ParserExeADD);
            --stack_size;
        }
        break;
    }
    case PARSER_MUL:
    {
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_SYMBOL)
        { // 3 * x
            if (p) {
                auto *t = new(p) ParserExeMUL_VP;
                p      += sizeof(ParserExeMUL_VP);
                t->i = parser_symbol_idx(node->r);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeMUL_VP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_NUMBER)
        { // 3 * f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeMUL_VN;
                p      += sizeof(ParserExeMUL_VN);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeMUL_VN);
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        { // x * y
            if (p) {
                auto *t = new(p) ParserExeMUL_PP;
                p      += sizeof(ParserExeMUL_PP);
                t->i1 = parser_symbol_idx(node->l);
                t->i2 = parser_symbol_idx(node->r);
            }
            exe_size += sizeof(ParserExeMUL_PP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_SYMBOL)
        { // x * f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeMUL_PN;
                p      += sizeof(ParserExeMUL_PN);
                t->i = parser_symbol_idx(node->l);
            }
            exe_size += sizeof(ParserExeMUL_PN);
        }
        else if (parser_node_equal(node->l,node->r))
        { // f(x) * f(x)
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                new(p)      ParserExeSquare;
                p += sizeof(ParserExeSquare);
            }
            exe_size += sizeof(ParserExeSquare);
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
            }
            if (p) {
                new(p)      ParserExeMUL;
                p += sizeof(ParserExeMUL);
            }
            exe_size += sizeof(ParserExeMUL);
            --stack_size;
        }
        break;
    }
    case PARSER_DIV:
    {
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_SYMBOL)
        { // 3 / x
            if (p) {
                auto *t = new(p) ParserExeDIV_VP;
                p      += sizeof(ParserExeDIV_VP);
                t->i = parser_symbol_idx(node->r);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeDIV_VP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_NUMBER)
        { // 3 / f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeDIV_VN;
                p      += sizeof(ParserExeDIV_VN);
                t->v = parser_get_number(node->l);
            }
            exe_size += sizeof(ParserExeDIV_VN);
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        { // x / y
            if (p) {
                auto *t = new(p) ParserExeDIV_PP;
                p      += sizeof(ParserExeDIV_PP);
                t->i1 = parser_symbol_idx(node->l);
                t->i2 = parser_symbol_idx(node->r);
            }
            exe_size += sizeof(ParserExeDIV_PP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (node->l->type == PARSER_SYMBOL)
        { // x / f(x)
            parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeDIV_PN;
                p      += sizeof(ParserExeDIV_PN);
                t->i = parser_symbol_idx(node->l);
                t->reverse = false;
            }
            exe_size += sizeof(ParserExeDIV_PN);
        }
        else if (node->r->type == PARSER_SYMBOL)
        { // f(x) / x
            parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeDIV_PN;
                p      += sizeof(ParserExeDIV_PN);
                t->i = parser_symbol_idx(node->r);
                t->reverse = true;
            }
            exe_size += sizeof(ParserExeDIV_PN);
        }
        else
        {
            int d1 = parser_ast_depth(node->l);
            int d2 = parser_ast_depth(node->r);
            if (d1 < d2) {
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeDIV_B;
                    p += sizeof(ParserExeDIV_B);
                }
                exe_size += sizeof(ParserExeDIV_B);
            } else {
                parser_compile_exe_size(node->l, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                parser_compile_exe_size(node->r, p, exe_size, max_stack_size,
                                        stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeDIV_F;
                    p += sizeof(ParserExeDIV_F);
                }
                exe_size += sizeof(ParserExeDIV_F);
            }
            --stack_size;
        }
        break;
    }
    case PARSER_F1:
    {
        if (((struct parser_f1*)node)->l->type == PARSER_SYMBOL)
        { // f(x)
            if (p) {
                auto *t = new(p) ParserExeF1_P;
                p      += sizeof(ParserExeF1_P);
                t->ftype = ((struct parser_f1*)node)->ftype;
                t->i = parser_symbol_idx(((struct parser_f1*)node)->l);
            }
            exe_size += sizeof(ParserExeF1_P);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else
        {
            parser_compile_exe_size(((struct parser_f1*)node)->l, p, exe_size,
                                    max_stack_size, stack_size, local_variables, ufs);
            if (p) {
                auto *t = new(p) ParserExeF1;
                p      += sizeof(ParserExeF1);
                t->ftype = ((struct parser_f1*)node)->ftype;
            }
            exe_size += sizeof(ParserExeF1);
        }
        break;
    }
    case PARSER_F2:
    {
        auto* f2 = (struct parser_f2*)node;
        bool l_is_symbol = (f2->l->type == PARSER_SYMBOL);
        bool r_is_number = (f2->r->type == PARSER_NUMBER);
        double rv = r_is_number ? parser_get_number(f2->r) : 0.0;
        if (f2->ftype == PARSER_POW && r_is_number && rv == 2.0)
        {
            if (l_is_symbol) { // x^2
                if (p) {
                    auto *t = new(p) ParserExeSquare_P;
                    p      += sizeof(ParserExeSquare_P);
                    t->i = parser_symbol_idx(f2->l);
                }
                exe_size += sizeof(ParserExeSquare_P);
                ++stack_size;
                max_stack_size = std::max(max_stack_size, stack_size);
            } else {
                parser_compile_exe_size(f2->l, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                if (p) {
                    new(p)      ParserExeSquare;
                    p += sizeof(ParserExeSquare);
                }
                exe_size += sizeof(ParserExeSquare);
            }
        }
        else if (f2->ftype == PARSER_POW && r_is_number &&
                 rv == std::floor(rv) &&
                 std::abs(rv) <= double(std::numeric_limits<int>::max()))
        {
            if (l_is_symbol) { // x^n
                if (p) {
                    auto *t = new(p) ParserExePOWI_P;
                    p      += sizeof(ParserExePOWI_P);
                    t->i = parser_symbol_idx(f2->l);
                    t->n = int(std::floor(rv));
                }
                exe_size += sizeof(ParserExePOWI_P);
                ++stack_size;
                max_stack_size = std::max(max_stack_size, stack_size);
            } else {
                parser_compile_exe_size(f2->l, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                if (p) {
                    auto *t = new(p) ParserExePOWI;
                    p      += sizeof(ParserExePOWI);
                    t->i = int(std::floor(rv));
                }
                exe_size += sizeof(ParserExePOWI);
            }
        }
        else if (f2->ftype == PARSER_POW && r_is_number && rv == 0.5)
        {
            if (l_is_symbol) { // sqrt(x)
                if (p) {
                    auto *t = new(p) ParserExeF1_P;
                    p      += sizeof(ParserExeF1_P);
                    t->ftype = PARSER_SQRT;
                    t->i = parser_symbol_idx(f2->l);
                }
                exe_size += sizeof(ParserExeF1_P);
                ++stack_size;
                max_stack_size = std::max(max_stack_size, stack_size);
            } else {
                parser_compile_exe_size(f2->l, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                if (p) {
                    auto *t = new(p) ParserExeF1;
                    p      += sizeof(ParserExeF1);
                    t->ftype = PARSER_SQRT;
                }
                exe_size += sizeof(ParserExeF1);
            }
        }
        else if (l_is_symbol && f2->r->type == PARSER_SYMBOL)
        { // f(x,y)
            if (p) {
                auto *t = new(p) ParserExeF2_PP;
                p      += sizeof(ParserExeF2_PP);
                t->ftype = f2->ftype;
                t->i1 = parser_symbol_idx(f2->l);
                t->i2 = parser_symbol_idx(f2->r);
            }
            exe_size += sizeof(ParserExeF2_PP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (l_is_symbol && r_is_number)
        { // f(x,3)
            if (p) {
                auto *t = new(p) ParserExeF2_PV;
                p      += sizeof(ParserExeF2_PV);
                t->ftype = f2->ftype;
                t->i = parser_symbol_idx(f2->l);
                t->v = rv;
            }
            exe_size += sizeof(ParserExeF2_PV);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else if (f2->l->type == PARSER_NUMBER && f2->r->type == PARSER_SYMBOL)
        { // f(3,x)
            if (p) {
                auto *t = new(p) ParserExeF2_VP;
                p      += sizeof(ParserExeF2_VP);
                t->ftype = f2->ftype;
                t->i = parser_symbol_idx(f2->r);
                t->v = parser_get_number(f2->l);
            }
            exe_size += sizeof(ParserExeF2_VP);
            ++stack_size;
            max_stack_size = std::max(max_stack_size, stack_size);
        }
        else
        {
            int d1 = parser_ast_depth(f2->l);
            int d2 = parser_ast_depth(f2->r);
            if (d1 < d2) {
                parser_compile_exe_size(f2->r, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                parser_compile_exe_size(f2->l, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                if (p) {
                    auto *t = new(p) ParserExeF2_B;
                    p      += sizeof(ParserExeF2_B);
                    t->ftype = f2->ftype;
                }
                exe_size += sizeof(ParserExeF2_B);
            } else {
                parser_compile_exe_size(f2->l, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                parser_compile_exe_size(f2->r, p, exe_size,
                                        max_stack_size, stack_size, local_variables, ufs);
                if (p) {
                    auto *t = new(p) ParserExeF2_F;
                    p      += sizeof(ParserExeF2_F);
                    t->ftype = f2->ftype;
                }
                exe_size += sizeof(ParserExeF2_F);
            }
            --stack_size;
        }
        break;
    }
    case PARSER_F3:
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(((struct parser_f3*)node)->ftype == PARSER_IF,
                                         "parser_compile: unknown f3 type");
        parser_compile_exe_size(((struct parser_f3*)node)->n1, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);

        ParserExeIF* tif = nullptr;
        char* psave = nullptr;
        if (p) {
            tif = new(p) ParserExeIF;
            p  += sizeof(ParserExeIF);
            psave = p;
        }
        exe_size += sizeof(ParserExeIF);
        --stack_size;
        auto stack_size_save = stack_size;

        parser_compile_exe_size(((struct parser_f3*)node)->n2, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);

        ParserExeJUMP* tjump = nullptr;
        if (p) {
            tjump = new(p) ParserExeJUMP;
            p +=    sizeof(ParserExeJUMP);
        }
        exe_size += sizeof(ParserExeJUMP);

        if (psave) {
            tif->offset = static_cast<int>(p-psave);
        }
        stack_size = stack_size_save;

        psave = p;
        parser_compile_exe_size(((struct parser_f3*)node)->n3, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);
        if (tjump) {
            tjump->offset = static_cast<int>(p-psave);
        }

        break;
    }
    case PARSER_ASSIGN:
    {
        auto *asgn = (struct parser_assign*)node;
        parser_compile_exe_size(asgn->v, p, exe_size, max_stack_size, stack_size,
                                local_variables, ufs);
        local_variables.push_back(asgn->s->name);
        break;
    }
    case PARSER_LIST:
    {
        parser_compile_exe_size(node->l, p, exe_size, max_stack_size, stack_size,
                                local_variables, ufs);
        parser_compile_exe_size(node->r, p, exe_size, max_stack_size, stack_size,
                                local_variables, ufs);
        break;
    }
    case PARSER_USRF1:
    {
        parser_compile_exe_size(((struct parser_usrf1*)node)->l, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);
        if (p) {
            auto *t = new(p) ParserExeUserFn;
            p      += sizeof(ParserExeUserFn);
            t->argc = 1;
            auto const& pp = ufs.at(std::string(((struct parser_usrf1*)node)->name));
            t->fh = pp.first;
            t->fd = pp.second;
        }
        exe_size += sizeof(ParserExeUserFn);
        break;
    }
    case PARSER_USRF2:
    {
        parser_compile_exe_size(((struct parser_usrf2*)node)->l, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);
        parser_compile_exe_size(((struct parser_usrf2*)node)->r, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);
        if (p) {
            auto *t = new(p) ParserExeUserFn;
            p      += sizeof(ParserExeUserFn);
            t->argc = 2;
            auto const& pp = ufs.at(std::string(((struct parser_usrf2*)node)->name));
            t->fh = pp.first;
            t->fd = pp.second;
        }
        exe_size += sizeof(ParserExeUserFn);
        --stack_size;
        break;
    }
    case PARSER_USRFN:
    {
        auto argc = ((struct parser_usrfn*)node)->argc;
        parser_compile_exe_size(((struct parser_usrfn*)node)->n1, p, exe_size,
                                max_stack_size, stack_size, local_variables, ufs);
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_compile_exe_size(((struct parser_usrfn*)node)->others[iarg],
                                    p, exe_size, max_stack_size, stack_size,
                                    local_variables, ufs);
        }
        if (p) {
            auto *t = new(p) ParserExeUserFn;
            p      += sizeof(ParserExeUserFn);
            t->argc = argc;
            auto const& pp = ufs.at(std::string(((struct parser_usrfn*)node)->name));
            t->fh = pp.first;
            t->fd = pp.second;
        }
        exe_size += sizeof(ParserExeUserFn);
        stack_size -= (argc-1);
        break;
    }
    default:
        amrex::Abort("parser_compile: unknown node type " + std::to_string(node->type));
    }
}

namespace {
    enum class paren_t {
        plusminus,
        muldiv,
        pow,
        atom
    };

    std::pair<bool,bool> need_parens (paren_t lhs, paren_t op, paren_t rhs)
    {
        std::pair<bool,bool> r;
        if (lhs < op) {
            r.first = true;
        } else if (lhs == op) {
            if (op == paren_t::pow) {
                r.first = true;
            } else {
                r.first = false;
            }
        } else {
            r.first = false;
        }
        if (op < rhs) {
            r.second = false;
        } else if (op == rhs) {
            if (op == paren_t::pow) {
                r.second = false;
            } else {
                r.second = true;
            }
        } else {
            r.second = true;
        }
        return r;
    }

    std::pair<std::string,paren_t> make_op_string (std::pair<std::string,paren_t> const& a,
                                                   std::pair<std::string,paren_t> const& op,
                                                   std::pair<std::string,paren_t> const& b)
    {
        auto pp = need_parens(a.second, op.second, b.second);
        std::string r;
        if (pp.first) {
            r += "(";
        }
        r.append(a.first);
        if (pp.first) {
            r += ")";
        }
        r += op.first;
        if (pp.second) {
            r += "(";
        }
        r.append(b.first);
        if (pp.second) {
            r += ")";
        }
        return {r,op.second};
    }

    std::pair<std::string,paren_t> make_f1_string (std::string_view const& f, std::string const& a)
    {
        std::string r{f};
        r.append("(").append(a).append(")");
        return {r,paren_t::atom};
    }

    std::pair<std::string,paren_t> make_f2_string (std::string_view const& f, std::string const& a,
                                                   std::string const& b)
    {
        std::string r{f};
        r.append("(").append(a).append(",").append(b).append(")");
        return {r,paren_t::atom};
    }
}

void parser_exe_print(char const* p, Vector<std::string> const& vars,
                      Vector<char const*> const& locals)
{
    Vector<std::pair<std::string,paren_t>> pstack;
    int count = 0;
    auto& os = amrex::OutStream();
    os << "  #     Instruction   StackSize   StackTop\n";

    auto get_sym = [&] (int i) -> std::pair<std::string,paren_t>
    {
        if (i >= AMREX_PARSER_LOCAL_IDX0) {
            return {locals[i-AMREX_PARSER_LOCAL_IDX0],paren_t::atom};
        } else {
            return {vars[i],paren_t::atom};
        }
    };

    auto get_val = [&] (double v) -> std::pair<std::string,paren_t>
    {
        return {std::to_string(v), paren_t::atom};
    };

    while (*((parser_exe_t*)p) != PARSER_EXE_NULL) { // NOLINT
        switch (*((parser_exe_t*)p))
        {
        case PARSER_EXE_NUMBER:
        {
            pstack.emplace_back(get_val(((ParserExeNumber*)p)->v));
            os << std::setw(3) << count++
               << std::setw(16) << "push"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeNumber);
            break;
        }
        case PARSER_EXE_SYMBOL:
        {
            int i = ((ParserExeSymbol*)p)->i;
            pstack.emplace_back(get_sym(i));
            os << std::setw(3) << count++
               << std::setw(16) << "push"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first;
            if (i >= AMREX_PARSER_LOCAL_IDX0) {
                os << " = " << pstack[i-AMREX_PARSER_LOCAL_IDX0].first;
            }
            os << "\n";
            p += sizeof(ParserExeSymbol);
            break;
        }
        case PARSER_EXE_ADD:
        {
            auto n = pstack.size();
            pstack[n-2] = make_op_string(pstack[n-2],{"+",paren_t::plusminus}, pstack[n-1]); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << "add"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeADD);
            break;
        }
        case PARSER_EXE_SUB_F:
        {
            auto n = pstack.size();
            pstack[n-2] = make_op_string(pstack[n-2], {"-",paren_t::plusminus}, pstack[n-1]); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << "sub"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSUB_F);
            break;
        }
        case PARSER_EXE_SUB_B:
        {
            auto n = pstack.size();
            pstack[n-2] = make_op_string(pstack[n-1], {"-",paren_t::plusminus}, pstack[n-2]); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << "rsub"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSUB_B);
            break;
        }
        case PARSER_EXE_MUL:
        {
            auto n = pstack.size();
            pstack[n-2] = make_op_string(pstack[n-2], {"*",paren_t::muldiv}, pstack[n-1]); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << "mul"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeMUL);
            break;
        }
        case PARSER_EXE_DIV_F:
        {
            auto n = pstack.size();
            pstack[n-2] = make_op_string(pstack[n-2], {"/",paren_t::muldiv}, pstack[n-1]); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << "div"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeDIV_F);
            break;
        }
        case PARSER_EXE_DIV_B:
        {
            auto n = pstack.size();
            pstack[n-2] = make_op_string(pstack[n-1], {"/",paren_t::muldiv}, pstack[n-2]); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << "rdiv"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeDIV_B);
            break;
        }
        case PARSER_EXE_F1:
        {
            pstack.back() = make_f1_string(parser_f1_s[((ParserExeF1*)p)->ftype],
                                           pstack.back().first);
            os << std::setw(3) << count++
               << std::setw(16) << parser_f1_s[((ParserExeF1*)p)->ftype]
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF1);
            break;
        }
        case PARSER_EXE_F2_F:
        {
            auto n = pstack.size();
            pstack[n-2] = make_f2_string(parser_f2_s[((ParserExeF2_F*)p)->ftype],
                                         pstack[n-2].first, pstack[n-1].first); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << parser_f2_s[((ParserExeF2_F*)p)->ftype]
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF2_F);
            break;
        }
        case PARSER_EXE_F2_B:
        {
            auto n = pstack.size();
            pstack[n-2] = make_f2_string(parser_f2_s[((ParserExeF2_B*)p)->ftype],
                                         pstack[n-1].first, pstack[n-2].first); // NOLINT
            pstack.pop_back();
            os << std::setw(3) << count++
               << std::setw(16) << std::string("r.")
                + std::string(parser_f2_s[((ParserExeF2_B*)p)->ftype])
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF2_B);
            break;
        }
        case PARSER_EXE_ADD_VP:
        {
            int i = ((ParserExeADD_VP*)p)->i;
            auto v = ((ParserExeADD_VP*)p)->v;
            pstack.push_back(make_op_string(get_val(v), {"+",paren_t::plusminus}, get_sym(i)));
            os << std::setw(3) << count++
               << std::setw(16) << "addvp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeADD_VP);
            break;
        }
        case PARSER_EXE_SUB_VP:
        {
            int i = ((ParserExeSUB_VP*)p)->i;
            auto v = ((ParserExeSUB_VP*)p)->v;
            pstack.push_back(make_op_string(get_val(v), {"-",paren_t::plusminus}, get_sym(i)));
            os << std::setw(3) << count++
               << std::setw(16) << "subvp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSUB_VP);
            break;
        }
        case PARSER_EXE_MUL_VP:
        {
            int i = ((ParserExeMUL_VP*)p)->i;
            auto v = ((ParserExeMUL_VP*)p)->v;
            pstack.push_back(make_op_string(get_val(v), {"*",paren_t::muldiv}, get_sym(i)));
            os << std::setw(3) << count++
               << std::setw(16) << "mulvp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeMUL_VP);
            break;
        }
        case PARSER_EXE_DIV_VP:
        {
            int i = ((ParserExeDIV_VP*)p)->i;
            auto v = ((ParserExeDIV_VP*)p)->v;
            pstack.push_back(make_op_string(get_val(v), {"/",paren_t::muldiv}, get_sym(i)));
            os << std::setw(3) << count++
               << std::setw(16) << "divvp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeDIV_VP);
            break;
        }
        case PARSER_EXE_ADD_PP:
        {
            int i = ((ParserExeADD_PP*)p)->i1;
            int j = ((ParserExeADD_PP*)p)->i2;
            pstack.push_back(make_op_string(get_sym(i), {"+",paren_t::plusminus}, get_sym(j)));
            os << std::setw(3) << count++
               << std::setw(16) << "addpp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeADD_PP);
            break;
        }
        case PARSER_EXE_SUB_PP:
        {
            int i = ((ParserExeSUB_PP*)p)->i1;
            int j = ((ParserExeSUB_PP*)p)->i2;
            pstack.push_back(make_op_string(get_sym(i), {"-",paren_t::plusminus}, get_sym(j)));
            os << std::setw(3) << count++
               << std::setw(16) << "subpp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSUB_PP);
            break;
        }
        case PARSER_EXE_MUL_PP:
        {
            int i = ((ParserExeMUL_PP*)p)->i1;
            int j = ((ParserExeMUL_PP*)p)->i2;
            pstack.push_back(make_op_string(get_sym(i), {"*",paren_t::muldiv}, get_sym(j)));
            os << std::setw(3) << count++
               << std::setw(16) << "mulpp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeMUL_PP);
            break;
        }
        case PARSER_EXE_DIV_PP:
        {
            int i = ((ParserExeDIV_PP*)p)->i1;
            int j = ((ParserExeDIV_PP*)p)->i2;
            pstack.push_back(make_op_string(get_sym(i), {"/",paren_t::muldiv}, get_sym(j)));
            os << std::setw(3) << count++
               << std::setw(16) << "divpp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeDIV_PP);
            break;
        }
        case PARSER_EXE_ADD_VN:
        {
            auto v = ((ParserExeADD_VN*)p)->v;
            pstack.back() = make_op_string(get_val(v), {"+",paren_t::plusminus}, pstack.back());
            os << std::setw(3) << count++
               << std::setw(16) << "addvn"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeADD_VN);
            break;
        }
        case PARSER_EXE_SUB_VN:
        {
            auto v = ((ParserExeSUB_VN*)p)->v;
            pstack.back() = make_op_string(get_val(v), {"-",paren_t::plusminus}, pstack.back());
            os << std::setw(3) << count++
               << std::setw(16) << "subvn"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSUB_VN);
            break;
        }
        case PARSER_EXE_MUL_VN:
        {
            auto v = ((ParserExeMUL_VN*)p)->v;
            pstack.back() = make_op_string(get_val(v), {"*",paren_t::muldiv}, pstack.back());
            os << std::setw(3) << count++
               << std::setw(16) << "mulvn"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeMUL_VN);
            break;
        }
        case PARSER_EXE_DIV_VN:
        {
            auto v = ((ParserExeDIV_VN*)p)->v;
            pstack.back() = make_op_string(get_val(v), {"/",paren_t::muldiv}, pstack.back());
            os << std::setw(3) << count++
               << std::setw(16) << "divvn"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeDIV_VN);
            break;
        }
        case PARSER_EXE_ADD_PN:
        {
            int i = ((ParserExeADD_PN*)p)->i;
            pstack.back() = make_op_string(get_sym(i), {"+",paren_t::plusminus}, pstack.back());
            os << std::setw(3) << count++
               << std::setw(16) << "addpn"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeADD_PN);
            break;
        }
        case PARSER_EXE_SUB_PN:
        {
            int i = ((ParserExeSUB_PN*)p)->i;
            std::string op;
            if (!((ParserExeSUB_PN*)p)->reverse) {
                pstack.back() = make_op_string(get_sym(i), {"-",paren_t::plusminus}, pstack.back());
            } else {
                pstack.back() = make_op_string(pstack.back(), {"-",paren_t::plusminus}, get_sym(i));
                op = "r";
            }
            op.append("subpn");
            os << std::setw(3) << count++
               << std::setw(16) << op
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p  += sizeof(ParserExeSUB_PN);
            break;
        }
        case PARSER_EXE_MUL_PN:
        {
            int i = ((ParserExeMUL_PN*)p)->i;
            pstack.back() = make_op_string(get_sym(i), {"*",paren_t::muldiv}, pstack.back());
            os << std::setw(3) << count++
               << std::setw(16) << "mulpn"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeMUL_PN);
            break;
        }
        case PARSER_EXE_DIV_PN:
        {
            int i = ((ParserExeDIV_PN*)p)->i;
            std::string op;
            if (((ParserExeDIV_PN*)p)->reverse) {
                pstack.back() = make_op_string(pstack.back(), {"/",paren_t::muldiv}, get_sym(i));
                op = "r";
            } else {
                pstack.back() = make_op_string(get_sym(i), {"/",paren_t::muldiv}, pstack.back());
            }
            op.append("divpn");
            os << std::setw(3) << count++
               << std::setw(16) << op
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeDIV_PN);
            break;
        }
        case PARSER_EXE_SQUARE:
        {
            pstack.back() = make_op_string(pstack.back(), {"^",paren_t::pow}, {"2",paren_t::atom});
            os << std::setw(3) << count++
               << std::setw(16) << "square"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSquare);
            break;
        }
        case PARSER_EXE_POWI:
        {
            int n = ((ParserExePOWI*)p)->i;
            pstack.back() = make_op_string(pstack.back(), {"^",paren_t::pow},
                                           {std::to_string(n),paren_t::atom});
            os << std::setw(3) << count++
               << std::setw(16) << "powi"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExePOWI);
            break;
        }
        case PARSER_EXE_FMA_VPV:
        {
            int i = ((ParserExeFMA_VPV*)p)->i;
            auto a = ((ParserExeFMA_VPV*)p)->a;
            auto b = ((ParserExeFMA_VPV*)p)->b;
            auto tmp = make_op_string(get_val(a), {"*",paren_t::muldiv}, get_sym(i));
            pstack.push_back(make_op_string(tmp, {"+",paren_t::plusminus}, get_val(b)));
            os << std::setw(3) << count++
               << std::setw(16) << "fmavpv"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeFMA_VPV);
            break;
        }
        case PARSER_EXE_F1_P:
        {
            auto* t = (ParserExeF1_P*)p;
            pstack.emplace_back(make_f1_string(parser_f1_s[t->ftype], get_sym(t->i).first));
            os << std::setw(3) << count++
               << std::setw(16) << std::string(parser_f1_s[t->ftype]) + "p"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF1_P);
            break;
        }
        case PARSER_EXE_SQUARE_P:
        {
            auto* t = (ParserExeSquare_P*)p;
            pstack.emplace_back(make_op_string(get_sym(t->i), {"^",paren_t::pow}, get_val(2.0)));
            os << std::setw(3) << count++
               << std::setw(16) << "squarep"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeSquare_P);
            break;
        }
        case PARSER_EXE_POWI_P:
        {
            auto* t = (ParserExePOWI_P*)p;
            pstack.emplace_back(make_op_string(get_sym(t->i), {"^",paren_t::pow},
                                               {std::to_string(t->n),paren_t::atom}));
            os << std::setw(3) << count++
               << std::setw(16) << "powip"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExePOWI_P);
            break;
        }
        case PARSER_EXE_F2_PP:
        {
            auto* t = (ParserExeF2_PP*)p;
            pstack.emplace_back(make_f2_string(parser_f2_s[t->ftype], get_sym(t->i1).first,
                                               get_sym(t->i2).first));
            os << std::setw(3) << count++
               << std::setw(16) << std::string(parser_f2_s[t->ftype]) + "pp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF2_PP);
            break;
        }
        case PARSER_EXE_F2_PV:
        {
            auto* t = (ParserExeF2_PV*)p;
            pstack.emplace_back(make_f2_string(parser_f2_s[t->ftype], get_sym(t->i).first,
                                               get_val(t->v).first));
            os << std::setw(3) << count++
               << std::setw(16) << std::string(parser_f2_s[t->ftype]) + "pv"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF2_PV);
            break;
        }
        case PARSER_EXE_F2_VP:
        {
            auto* t = (ParserExeF2_VP*)p;
            pstack.emplace_back(make_f2_string(parser_f2_s[t->ftype], get_val(t->v).first,
                                               get_sym(t->i).first));
            os << std::setw(3) << count++
               << std::setw(16) << std::string(parser_f2_s[t->ftype]) + "vp"
               << std::setw(12) << pstack.size()
               << "   "
               << pstack.back().first << "\n";
            p += sizeof(ParserExeF2_VP);
            break;
        }
        case PARSER_EXE_IF:
        {
            os << "parser_exe_print: cannot handle IF yet\n";
            return;
        }
        case PARSER_EXE_JUMP:
        {
            os << "parser_exe_print: cannot handle JUMP yet\n";
            return;
        }
        case PARSER_EXE_USER_FN:
        {
            os << "parser_exe_print: cannot handle user function yet\n";
            return;
        }
        default:
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,"parser_print: unknown node type");
        }
    }
}

}
