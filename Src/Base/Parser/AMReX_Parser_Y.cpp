#include <AMReX.H>
#include <AMReX_Parser_Y.H>
#include <amrex_parser.tab.h>

#include <algorithm>
#include <cstdarg>
#include <vector>

void
amrex_parsererror (char const *s, ...)
{
    char print_buff[512];
    std::va_list vl;
    va_start(vl, s);
    std::vsnprintf(print_buff, 512, s, vl);
    va_end(vl);
    throw std::runtime_error(print_buff);
}

namespace amrex {

namespace {
    struct parser_node* parser_root = nullptr;
    std::vector<void*>  parser_ptrs;
}

// This is called by a bison rule to store the original AST in a static variable.
void
parser_defexpr (struct parser_node* body)
{
    parser_root = body;
}

struct parser_symbol*
parser_makesymbol (char* name)
{
    // We allocate more than enough space so that late we can turn parser_symbol
    // into into parser_node if necessary.
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *symbol = (struct parser_symbol*) parser_ptrs.back(); // NOLINT
    symbol->type = PARSER_SYMBOL;
    symbol->name = strdup(name);
    parser_ptrs.push_back(symbol->name);
    symbol->ip = -1;
    return symbol;
}

struct parser_node*
parser_newnode (enum parser_node_t type, struct parser_node* l, struct parser_node* r)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *tmp = (struct parser_node*) parser_ptrs.back();
    if (type == PARSER_SUB) {
        tmp->type = PARSER_ADD;
        tmp->l = l;
        tmp->r = parser_newnode(PARSER_MUL, parser_newnumber(-1.0), r);
    } else {
        tmp->type = type;
        tmp->l = l;
        tmp->r = r;
    }
    return tmp;
}

struct parser_node*
parser_newneg (struct parser_node* n)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *tmp = (struct parser_node*) parser_ptrs.back();
    tmp->type = PARSER_MUL;
    tmp->l = parser_newnumber(-1.0);
    tmp->r = n;
    return tmp;
}

struct parser_node*
parser_newnumber (double d)
{
    // We allocate more than enough space so that late we can turn parser_number
    // into into parser_node if necessary.
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *r = (struct parser_number*) parser_ptrs.back(); // NOLINT
    r->type = PARSER_NUMBER;
    r->value = d;
    return (struct parser_node*) r;
}

struct parser_node*
parser_newsymbol (struct parser_symbol* symbol)
{
    return (struct parser_node*) symbol;
}

struct parser_node*
parser_newf1 (enum parser_f1_t ftype, struct parser_node* l)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *tmp = (struct parser_f1*) parser_ptrs.back(); // NOLINT
    tmp->type = PARSER_F1;
    tmp->l = l;
    tmp->ftype = ftype;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newf2 (enum parser_f2_t ftype, struct parser_node* l, struct parser_node* r)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *tmp = (struct parser_f2*) parser_ptrs.back(); // NOLINT
    tmp->type = PARSER_F2;
    tmp->l = l;
    tmp->r = r;
    tmp->ftype = ftype;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newf3 (enum parser_f3_t ftype, struct parser_node* n1, struct parser_node* n2,
              struct parser_node* n3)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *tmp = (struct parser_f3*) parser_ptrs.back(); // NOLINT
    tmp->type = PARSER_F3;
    tmp->n1 = n1;
    tmp->n2 = n2;
    tmp->n3 = n3;
    tmp->ftype = ftype;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newusrf1 (struct parser_symbol* fname, struct parser_node* l)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto* tmp = (struct parser_usrf1*) parser_ptrs.back(); // NOLINT
    tmp->type = PARSER_USRF1;
    tmp->name = strdup(fname->name);
    parser_ptrs.push_back(tmp->name);
    tmp->l = l;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newusrf2 (struct parser_symbol* fname, struct parser_node* l,
                 struct parser_node* r)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto* tmp = (struct parser_usrf2*) parser_ptrs.back(); // NOLINT
    tmp->type = PARSER_USRF2;
    tmp->name = strdup(fname->name);
    parser_ptrs.push_back(tmp->name);
    tmp->l = l;
    tmp->r = r;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newusrfn (struct parser_symbol* fname, std::vector<struct parser_node*> const& nv)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto* tmp = (struct parser_usrfn*) parser_ptrs.back(); // NOLINT
    tmp->type = PARSER_USRFN;
    tmp->argc = short(nv.size());
    tmp->name = strdup(fname->name);
    parser_ptrs.push_back(tmp->name);
    tmp->n1 = nv[0];
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node*)*(tmp->argc-1)));
    tmp->others = (struct parser_node**) parser_ptrs.back(); // NOLINT
    for (short iarg = 0; iarg < tmp->argc-1; ++iarg) {
        tmp->others[iarg] = nv[iarg+1];
    }
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newassign (struct parser_symbol* sym, struct parser_node* v)
{
    parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
    auto *r = (struct parser_assign*) parser_ptrs.back(); // NOLINT
    r->type = PARSER_ASSIGN;
    r->s = sym;
    r->v = v;
    return (struct parser_node*) r;
}

struct parser_node*
parser_newlist (struct parser_node* nl, struct parser_node* nr)
{
    if (nr == nullptr) {
        return nl;
    } else {
        parser_ptrs.push_back(std::malloc(sizeof(struct parser_node)));
        auto *r = (struct parser_node*) parser_ptrs.back();
        r->type = PARSER_LIST;
        r->l = nl;
        r->r = nr;
        return r;
    }
}

namespace {
struct parser_node* parser_get_rightmost_operand (struct parser_node* node)
{
    if (node && node->type == PARSER_F2) {
        auto ftype = ((struct parser_f2*)node)->ftype;
        if (ftype == PARSER_CMP_CHAIN) {
            return parser_get_rightmost_operand(node->r);
        } else if (ftype == PARSER_LT || ftype == PARSER_GT ||
                   ftype == PARSER_LEQ || ftype == PARSER_GEQ ||
                   ftype == PARSER_EQ || ftype == PARSER_NEQ) {
            return node->r;
        }
    }
    return nullptr;
}

bool parser_is_comparison (struct parser_node* node)
{
    if (node && node->type == PARSER_F2) {
        auto ftype = ((struct parser_f2*)node)->ftype;
        return (ftype == PARSER_LT || ftype == PARSER_GT ||
                ftype == PARSER_LEQ || ftype == PARSER_GEQ ||
                ftype == PARSER_EQ || ftype == PARSER_NEQ ||
                (ftype == PARSER_CMP_CHAIN && parser_is_comparison(node->r)));
    } else {
        return false;
    }
}
}

struct parser_node* parser_newcmpchain (struct parser_node* nl, enum parser_f2_t cmp,
                                        struct parser_node* nr)
{
    /* If left side is already a comparison, this extends the chain */
    if (amrex::parser_is_comparison(nl)) {
       return amrex::parser_newf2(amrex::PARSER_CMP_CHAIN, nl,
                                  amrex::parser_newf2(cmp,
                                                      amrex::parser_get_rightmost_operand(nl),
                                                      nr));
    } else {
        return amrex::parser_newf2(cmp, nl, nr); // Initial comparison
    }
}

/*******************************************************************/

struct amrex_parser*
amrex_parser_new ()
{
    auto *my_parser = (struct amrex_parser*) std::malloc(sizeof(struct amrex_parser));

    my_parser->sz_mempool = parser_ast_size(parser_root);
    my_parser->p_root = std::malloc(my_parser->sz_mempool);
    my_parser->p_free = my_parser->p_root;

    my_parser->ast = parser_ast_dup(my_parser, parser_root);

    amrex_parser_delete_ptrs();

    if ((char*)my_parser->p_root + my_parser->sz_mempool != (char*)my_parser->p_free) {
        amrex::Abort("amrex_parser_new: error in memory size");
    }

    std::map<std::string,double> local_consts;
    parser_ast_optimize(my_parser->ast, local_consts);
    parser_ast_sort(my_parser->ast);

    return my_parser;
}

void
amrex_parser_delete (struct amrex_parser* parser)
{
    std::free(parser->p_root);
    std::free(parser);
}

void
amrex_parser_delete_ptrs ()
{
    for (auto* p : parser_ptrs) {
        std::free(p);
    }
    parser_ptrs.clear();
}

namespace {

std::size_t
parser_aligned_size (std::size_t N)
{
    const unsigned int align_size = 16;
    std::size_t x = N + (align_size-1);
    x -= x & (align_size-1);
    return x;
}

void*
parser_allocate (struct amrex_parser* my_parser, std::size_t N)
{
    void* r = my_parser->p_free;
    my_parser->p_free = (char*)r + parser_aligned_size(N);
    return r;
}

}

std::size_t
parser_ast_size (struct parser_node* node)
{
    std::size_t result = 0;

    switch (node->type)
    {
    case PARSER_NUMBER:
        result = parser_aligned_size(sizeof(struct parser_node));
        break;
    case PARSER_SYMBOL:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_aligned_size(std::strlen(((struct parser_symbol*)node)->name)+1);
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_ast_size(node->l) + parser_ast_size(node->r);
        break;
    case PARSER_F1:
        result = parser_aligned_size(sizeof(struct parser_node))
            +             parser_ast_size(((struct parser_f1*)node)->l);
        break;
    case PARSER_F2:
        result = parser_aligned_size(sizeof(struct parser_node))
            +             parser_ast_size(((struct parser_f2*)node)->l)
            +             parser_ast_size(((struct parser_f2*)node)->r);
        break;
    case PARSER_F3:
        result = parser_aligned_size(sizeof(struct parser_node))
            +             parser_ast_size(((struct parser_f3*)node)->n1)
            +             parser_ast_size(((struct parser_f3*)node)->n2)
            +             parser_ast_size(((struct parser_f3*)node)->n3);
        break;
    case PARSER_USRF1:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_aligned_size(std::strlen(((struct parser_usrf1*)node)->name)+1)
            + parser_ast_size(((struct parser_usrf1*)node)->l);
        break;
    case PARSER_USRF2:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_aligned_size(std::strlen(((struct parser_usrf2*)node)->name)+1)
            + parser_ast_size(((struct parser_usrf2*)node)->l)
            + parser_ast_size(((struct parser_usrf2*)node)->r);
        break;
    case PARSER_USRFN:
    {
        short argc = ((struct parser_usrfn*)node)->argc;
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_aligned_size(std::strlen(((struct parser_usrfn*)node)->name)+1)
            + parser_aligned_size(sizeof(struct parser_node*)*(argc-1))
            + parser_ast_size(((struct parser_usrfn*)node)->n1);
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            result += parser_ast_size(((struct parser_usrfn*)node)->others[iarg]);
        }
        break;
    }
    case PARSER_ASSIGN:
        result += parser_aligned_size(sizeof(struct parser_node))
            + parser_ast_size((struct parser_node*)(((struct parser_assign*)node)->s))
            + parser_ast_size(((struct parser_assign*)node)->v);
        break;
    default:
        amrex::Abort("parser_ast_size: unknown node type " + std::to_string(node->type));
    }

    return result;
}

struct parser_node*
parser_ast_dup (struct amrex_parser* my_parser, struct parser_node* node)
{
    void* result = nullptr;

    switch (node->type)
    {
    case PARSER_NUMBER:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_number));
        break;
    case PARSER_SYMBOL:
    {
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_symbol));
        const auto len = std::strlen(((struct parser_symbol*)node)->name)+1;
        ((struct parser_symbol*)result)->name = (char*) parser_allocate
            (my_parser, len);
        std::strncpy(((struct parser_symbol*)result)->name,
                     ((struct parser_symbol*)node  )->name, len);
        break;
    }
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_node));
        ((struct parser_node*)result)->l = parser_ast_dup(my_parser, node->l);
        ((struct parser_node*)result)->r = parser_ast_dup(my_parser, node->r);
        break;
    case PARSER_F1:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_f1));
        ((struct parser_f1*)result)->l = parser_ast_dup(my_parser,
                                                 ((struct parser_f1*)node)->l);
        break;
    case PARSER_F2:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_f2));
        ((struct parser_f2*)result)->l = parser_ast_dup(my_parser,
                                                 ((struct parser_f2*)node)->l);
        ((struct parser_f2*)result)->r = parser_ast_dup(my_parser,
                                                 ((struct parser_f2*)node)->r);
        break;
    case PARSER_F3:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_f3));
        ((struct parser_f3*)result)->n1 = parser_ast_dup(my_parser,
                                                 ((struct parser_f3*)node)->n1);
        ((struct parser_f3*)result)->n2 = parser_ast_dup(my_parser,
                                                 ((struct parser_f3*)node)->n2);
        ((struct parser_f3*)result)->n3 = parser_ast_dup(my_parser,
                                                 ((struct parser_f3*)node)->n3);
        break;
    case PARSER_USRF1:
    {
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_usrf1));
        const auto len = std::strlen(((struct parser_usrf1*)node)->name)+1;
        ((struct parser_usrf1*)result)->name = (char*) parser_allocate(my_parser, len);
        std::strncpy(((struct parser_usrf1*)result)->name,
                     ((struct parser_usrf1*)node  )->name, len);
        ((struct parser_usrf1*)result)->l = parser_ast_dup
            (my_parser, ((struct parser_usrf1*)node)->l);
        break;
    }
    case PARSER_USRF2:
    {
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_usrf2));
        const auto len = std::strlen(((struct parser_usrf2*)node)->name)+1;
        ((struct parser_usrf2*)result)->name = (char*) parser_allocate(my_parser, len);
        std::strncpy(((struct parser_usrf2*)result)->name,
                     ((struct parser_usrf2*)node  )->name, len);
        ((struct parser_usrf2*)result)->l = parser_ast_dup
            (my_parser, ((struct parser_usrf2*)node)->l);
        ((struct parser_usrf2*)result)->r = parser_ast_dup
            (my_parser, ((struct parser_usrf2*)node)->r);
        break;
    }
    case PARSER_USRFN:
    {
        short argc = ((struct parser_usrfn*)node)->argc;
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_usrfn));
        const auto len = std::strlen(((struct parser_usrfn*)node)->name)+1;
        ((struct parser_usrfn*)result)->name = (char*) parser_allocate(my_parser, len);
        std::strncpy(((struct parser_usrfn*)result)->name,
                     ((struct parser_usrfn*)node  )->name, len);
        ((struct parser_usrfn*)result)->others = (struct parser_node**)
            parser_allocate(my_parser,     sizeof(struct parser_node*)*(argc-1));
        ((struct parser_usrfn*)result)->n1 = parser_ast_dup
            (my_parser, ((struct parser_usrfn*)node)->n1);
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            ((struct parser_usrfn*)result)->others[iarg] = parser_ast_dup
                (my_parser, ((struct parser_usrfn*)node)->others[iarg]);
        }
        break;
    }
    case PARSER_ASSIGN:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_assign));
        ((struct parser_assign*)result)->s = (struct parser_symbol*)
            parser_ast_dup(my_parser, (struct parser_node*)
                                                (((struct parser_assign*)node)->s));
        ((struct parser_assign*)result)->v = parser_ast_dup(my_parser,
                                                 ((struct parser_assign*)node)->v);
        break;
    default:
        amrex::Abort("parser_ast_dup: unknown node type " + std::to_string(node->type));
    }
    return (struct parser_node*)result;
}

namespace {
    char* parser_get_name (struct parser_node* node)
    {
        AMREX_ASSERT(node->type == PARSER_SYMBOL);
        return ((struct parser_symbol*)node)->name;
    }

    bool is_add_combinable (struct parser_node* a, struct parser_node*b)
    {
        if ((a->type == PARSER_NUMBER) &&
            (b->type == PARSER_NUMBER))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if (parser_node_equal(a, b))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_MUL) &&
                 (a->l->type == PARSER_NUMBER) &&
                 parser_node_equal(a->r, b))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((b->type == PARSER_MUL) &&
                 (b->l->type == PARSER_NUMBER) &&
                 parser_node_equal(a, b->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_MUL) &&
                 (b->type == PARSER_MUL) &&
                 (a->l->type == PARSER_NUMBER) &&
                 (b->l->type == PARSER_NUMBER) &&
                 parser_node_equal(a->r, b->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_DIV) &&
                 (b->type == PARSER_DIV) &&
                 (a->l->type == PARSER_NUMBER) &&
                 (b->l->type == PARSER_NUMBER) &&
                 parser_node_equal(a->r, b->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else
        {
            return false;
        }
    }

    bool is_mul_combinable (struct parser_node*a, struct parser_node*b)
    {
        if ((a->type == PARSER_NUMBER) &&
            (b->type == PARSER_NUMBER))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_NUMBER) &&
                 (b->type == PARSER_MUL) &&
                 (b->l->type == PARSER_NUMBER))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_NUMBER) &&
                 (b->type == PARSER_DIV) &&
                 (b->l->type == PARSER_NUMBER))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((b->type == PARSER_NUMBER) &&
                 (a->type == PARSER_MUL) &&
                 (a->l->type == PARSER_NUMBER))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_MUL) &&
                 (b->type == PARSER_MUL) &&
                 (a->l->type == PARSER_NUMBER) &&
                 (b->l->type == PARSER_NUMBER))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((b->type == PARSER_DIV) &&
                 parser_node_equal(a, b->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_MUL) &&
                 (b->type == PARSER_DIV) &&
                 parser_node_equal(a->l, b->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_MUL) &&
                 (b->type == PARSER_DIV) &&
                 parser_node_equal(a->r, b->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((a->type == PARSER_DIV) &&
                 parser_node_equal(a->r, b))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((b->type == PARSER_MUL) &&
                 (a->type == PARSER_DIV) &&
                 parser_node_equal(b->l, a->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if ((b->type == PARSER_MUL) &&
                 (a->type == PARSER_DIV) &&
                 parser_node_equal(b->r, a->r))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if (b->type == PARSER_F2 &&
                 ((struct parser_f2*)b)->ftype == PARSER_POW &&
                 parser_node_equal(((struct parser_f2*)b)->l, a))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if (a->type == PARSER_F2 &&
                 ((struct parser_f2*)a)->ftype == PARSER_POW &&
                 parser_node_equal(((struct parser_f2*)a)->l, b))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else
        {
            return false;
        }
    }

    bool parser_node_compare (struct parser_node* a, struct parser_node* b)
    {
        if ((a->type) < (b->type)) {
            return true;
        } else if ((a->type) == (b->type)) {
            switch (a->type)
            {
            case PARSER_NUMBER:
                return parser_get_number(a) < parser_get_number(b);
            case PARSER_SYMBOL:
                return std::strcmp(parser_get_name(a),
                                   parser_get_name(b)) < 0;
            case PARSER_ADD:
            case PARSER_SUB:
            case PARSER_MUL:
            case PARSER_DIV:
                return parser_node_compare(a->r, b->r) ||
                    (parser_node_equal(a->r, b->r) &&
                     parser_node_compare(a->l, b->l));
            case PARSER_F1:
                return (((struct parser_f1*)(a))->ftype <
                        ((struct parser_f1*)(b))->ftype) ||
                    ((((struct parser_f1*)(a))->ftype ==
                      ((struct parser_f1*)(b))->ftype) &&
                     parser_node_compare(a->l,b->l));
            case PARSER_F2:
                if (((struct parser_f2*)(a))->ftype <
                    ((struct parser_f2*)(b))->ftype) {
                    return true;
                } else if (((struct parser_f2*)(a))->ftype ==
                           ((struct parser_f2*)(b))->ftype) {
                    return parser_node_compare(a->r, b->r) ||
                        (parser_node_equal(a->r, b->r) &&
                         parser_node_compare(a->l, b->l));
                } else {
                    return false;
                }
            default:
                return false;
            }
        }
        return false;
    }

    template <typename F>
    bool group_combinables (struct parser_node*& a, struct parser_node*& b,
                            F const& f, parser_node_t type)
    {
        if (a->type == type && f(a->l, b))
        {
            std::swap(a->r,b);
            return true;
        }
        else if (a->type == type && f(a->r, b))
        {
            std::swap(a->l,b);
            return true;
        }
        else if (b->type == type && f(a, b->l))
        {
            std::swap(a, b->r);
            return true;
        }
        else if (b->type == type && f(a, b->r))
        {
            std::swap(a, b->l);
            return true;
        }
        else if (a->type == type && group_combinables(a->l, b, f, type))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if (a->type == type && group_combinables(a->r, b, f, type))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if (b->type == type && group_combinables(a, b->l, f, type))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        else if (b->type == type && group_combinables(a, b->r, f, type))
        { // NOLINT(bugprone-branch-clone)
            return true;
        }
        return false;
    }

    // Try to find x**idx in y, where idx is either 1 or -1. If found, we
    // can cancel them by setting x and the x**idx node in y to one.
    bool try_cancel (struct parser_node* x, int idx, struct parser_node* y)
    {
        if (idx == 1 && parser_node_equal(x,y)) {
            parser_set_number(x, 1.0);
            parser_set_number(y, 1.0);
            return true;
        } else if (y->type == PARSER_MUL) {
            if (try_cancel(x, idx, y->l)) {
                return true;
            } else if (try_cancel(x, idx, y->r)) {
                return true;
            }
        } else if (y->type == PARSER_DIV) {
            if (try_cancel(x, idx, y->l)) {
                return true;
            } else if (try_cancel(x, -idx, y->r)) {
                return true;
            }
        } else if (x->type == PARSER_MUL) {
            if (try_cancel(x->l, idx, y)) {
                return true;
            } else if (try_cancel(x->r, idx, y)) {
                return true;
            }
        } else if (x->type == PARSER_DIV) {
            if (try_cancel(x->l, idx, y)) {
                return true;
            } else if (try_cancel(x->r, -idx, y)) {
                return true;
            }
        }
        return false;
    }
}

bool parser_same_symbol (struct parser_node* a, struct parser_node* b)
{
    return (a->type == PARSER_SYMBOL)
        && (b->type == PARSER_SYMBOL)
        && (std::strcmp(((struct parser_symbol*)a)->name,
                        ((struct parser_symbol*)b)->name) == 0);
}

bool parser_node_equal (struct parser_node* a, struct parser_node* b)
{
    if (a->type != b->type) { return false; }
    switch (a->type)
    {
    case PARSER_NUMBER:
        return parser_get_number(a) == parser_get_number(b);
    case PARSER_SYMBOL:
        return parser_same_symbol(a,b);
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
        return parser_node_equal(a->l,b->l) && parser_node_equal(a->r,b->r);
    case PARSER_F1:
        return (((struct parser_f1*)a)->ftype == ((struct parser_f1*)b)->ftype)
            && parser_node_equal(((struct parser_f1*)a)->l,
                                 ((struct parser_f1*)b)->l);
    case PARSER_F2:
        return (((struct parser_f2*)a)->ftype == ((struct parser_f2*)b)->ftype)
            && parser_node_equal(((struct parser_f2*)a)->l,
                                 ((struct parser_f2*)b)->l)
            && parser_node_equal(((struct parser_f2*)a)->r,
                                 ((struct parser_f2*)b)->r);
    case PARSER_F3:
        return (((struct parser_f3*)a)->ftype == ((struct parser_f3*)b)->ftype)
            && parser_node_equal(((struct parser_f3*)a)->n1,
                                 ((struct parser_f3*)b)->n1)
            && parser_node_equal(((struct parser_f3*)a)->n2,
                                 ((struct parser_f3*)b)->n2)
            && parser_node_equal(((struct parser_f3*)a)->n3,
                                 ((struct parser_f3*)b)->n3);
    case PARSER_USRF1:
        return (std::strcmp(((struct parser_usrf1*)a)->name,
                            ((struct parser_usrf1*)b)->name) == 0)
            && parser_node_equal(((struct parser_usrf1*)a)->l,
                                 ((struct parser_usrf1*)b)->l);
    case PARSER_USRF2:
        return (std::strcmp(((struct parser_usrf2*)a)->name,
                            ((struct parser_usrf2*)b)->name) == 0)
            && parser_node_equal(((struct parser_usrf2*)a)->l,
                                 ((struct parser_usrf2*)b)->l)
            && parser_node_equal(((struct parser_usrf2*)a)->r,
                                 ((struct parser_usrf2*)b)->r);
    case PARSER_USRFN:
    {
        short argca = ((struct parser_usrfn*)a)->argc;
        short argcb = ((struct parser_usrfn*)b)->argc;
        bool r = (argca == argcb)
            && (std::strcmp(((struct parser_usrfn*)a)->name,
                            ((struct parser_usrfn*)b)->name) == 0)
            && parser_node_equal(((struct parser_usrfn*)a)->n1,
                                 ((struct parser_usrfn*)b)->n1);
        for (short iarg = 0; (iarg < argca-1) && r; ++iarg) {
            r = r && parser_node_equal(((struct parser_usrfn*)a)->others[iarg],
                                       ((struct parser_usrfn*)b)->others[iarg]);
        }
        return r;
    }
    case PARSER_LIST:
    case PARSER_ASSIGN:
        return false;
    default:
        amrex::Abort("parser_node_equal: unknown node type " + std::to_string(a->type));
        return false;
    }
}

void
parser_ast_optimize (struct parser_node*& node, std::map<std::string,double>& local_consts)
{
    // No need to free memory because we only call this on ASTs in
    // amrex_parser that are allocated from the memory pool.

    if (node == nullptr) { return; };

    switch (node->type)
    {
    case PARSER_NUMBER:
        break;
    case PARSER_SYMBOL:
        if (auto it = local_consts.find(((struct parser_symbol*)node)->name);
            it != local_consts.end())
        {
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = it->second;
        }
        break;
    case PARSER_ADD:
        parser_ast_optimize(node->l,local_consts);
        parser_ast_optimize(node->r,local_consts);
        parser_ast_sort(node);
        if (node->l->type == PARSER_NUMBER && parser_get_number(node->l) == 0.0)
        { // 0 + ?
            std::memcpy(node, node->r, sizeof(struct parser_node));
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_NUMBER)
        { // 3 + 4 => 7
            double a = parser_get_number(node->l);
            double b = parser_get_number(node->r);
            parser_set_number(node, a+b);
        }
        else if (parser_node_equal(node->l, node->r))
        { // x + x = 2*x
            parser_set_number(node->l, 2.0);
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 parser_node_equal(node->l->r, node->r))
        { // (3 * x) + x => 4 * x
            parser_set_number(node->l, parser_get_number(node->l->l)+1.0);
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 parser_node_equal(node->r->r, node->l))
        { // x + (3 * x) => 4 * x
            parser_set_number(node->r, parser_get_number(node->r->l)+1.0);
            std::swap(node->l, node->r);
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER &&
                 parser_node_equal(node->l->r, node->r->r))
        { // (3*x) + (4*x) => 7 * x
            double c = parser_get_number(node->l->l) + parser_get_number(node->r->l);
            parser_set_number(node->l, c);
            node->r = node->r->r;
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_DIV &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV &&
                 node->r->l->type == PARSER_NUMBER &&
                 parser_node_equal(node->l->r, node->r->r))
        { // (3/x) + (4/x) => 7 / x
            double c = parser_get_number(node->l->l) + parser_get_number(node->r->l);
            parser_set_number(node->l, c);
            node->r = node->r->r;
            node->type = PARSER_DIV;
            parser_ast_optimize(node,local_consts);
        }
        // At this point, we have handled all directly combinable cases.
        else if (group_combinables(node->l, node->r, is_add_combinable, PARSER_ADD))
        {
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type != PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 parser_node_equal(node->l, node->r->l))
        { // x + x*y => x*(1.+y)
            parser_set_number(node->r->l, 1.0);
            node->r->type = PARSER_ADD;
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type != PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 parser_node_equal(node->l, node->r->r))
        { // x + y*x => x*(1+y)
            parser_set_number(node->r->r, 1.0);
            std::swap(node->r->l, node->r->r);
            node->r->type = PARSER_ADD;
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type != PARSER_NUMBER &&
                 node->l->type == PARSER_MUL &&
                 parser_node_equal(node->l->l, node->r))
        { // x*y + x => x*(1.+y)
            std::swap(node->l, node->r);
            node->type = PARSER_MUL;
            parser_set_number(node->r->l, 1.0);
            node->r->type = PARSER_ADD;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type != PARSER_NUMBER &&
                 node->l->type == PARSER_MUL &&
                 parser_node_equal(node->r, node->l->r))
        { // y*x + x => x*(1+y)
            std::swap(node->l, node->r);
            node->type = PARSER_MUL;
            std::swap(node->r->l, node->r->r);
            parser_set_number(node->r->l, 1.0);
            node->r->type = PARSER_ADD;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->r->type == PARSER_MUL &&
                 node->l->l->type != PARSER_NUMBER &&
                 parser_node_equal(node->l->l, node->r->l))
        { // a*x + a*y = a*(x+y)
            std::swap(node->l->r, node->r->l); // (a*a) + (x*y)
            node->type = PARSER_MUL;           // (a*a) * (x*y)
            node->r->type = PARSER_ADD;        // (a*a) * (x+y)
            node->l = node->l->l;              // a * (x+y)
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->r->type == PARSER_MUL &&
                 parser_node_equal(node->l->l, node->r->r))
        { // a*x + y*a
            std::swap(node->l->r, node->r->r); // (a*a) + (y*x)
            node->type = PARSER_MUL;           // (a*a) * (y*x)
            node->r->type = PARSER_ADD;        // (a*x) * (y+x)
            node->l = node->l->l;              // a * (y+x)
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->r->type == PARSER_MUL &&
                 parser_node_equal(node->l->r, node->r->l))
        { // x*a + a*y
            std::swap(node->l->l, node->r->l); // (a*a) + (x*y)
            node->type = PARSER_MUL;           // (a*a) * (x*y)
            node->r->type = PARSER_ADD;        // (a*a) * (x+y)
            node->l = node->l->l;              // a * (x+y)
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->r->type == PARSER_MUL &&
                 parser_node_equal(node->l->r, node->r->r))
        { // x*a + y*a
            std::swap(node->l->l, node->r->r); // (a*a) + (y*x)
            node->type = PARSER_MUL;           // (a*x) * (y*x)
            node->r->type = PARSER_ADD;        // (a*x) * (y+x)
            node->l = node->l->l;              // a * (y+x)
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_ADD &&
                 node->l->type != PARSER_NUMBER &&
                 node->r->l->type == PARSER_NUMBER)
        { // L + (# + RR)
            std::swap(node->l, node->r->l);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_ADD &&
                 node->r->type != PARSER_NUMBER &&
                 node->l->l->type == PARSER_NUMBER)
        { // (# + LR) + R
            std::swap(node->l->l,node->r);
            std::swap(node->l, node->r);
            parser_ast_optimize(node,local_consts);
        }
        break;
    case PARSER_MUL:
        parser_ast_optimize(node->l,local_consts);
        parser_ast_optimize(node->r,local_consts);
        parser_ast_sort(node);
        if (node->l->type == PARSER_NUMBER && parser_get_number(node->l) == 0.0)
        {
            parser_set_number(node, 0.0);
        }
        else if (node->l->type == PARSER_NUMBER && parser_get_number(node->l) == 1.0)
        {
            std::memcpy(node, node->r, sizeof(struct parser_node));
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_NUMBER)
        {
            parser_set_number(node, parser_get_number(node->l)
                              *     parser_get_number(node->r));
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER)
        { // 3 * (4*x)
            parser_set_number(node->l, parser_get_number(node->l)
                              *        parser_get_number(node->r->l));
            node->r = node->r->r;
            parser_ast_optimize(node,local_consts); // The number might be 1.
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV &&
                 node->r->l->type == PARSER_NUMBER)
        { // 3 * (4/x)
            parser_set_number(node->l, parser_get_number(node->l)
                              *        parser_get_number(node->r->l));
            node->type = PARSER_DIV;
            node->r = node->r->r;
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER)
        { // (3*x) * (4*y)
            std::swap(node->l->r, node->r->l); // (3*4) * (x*y)
            parser_set_number(node->l, parser_get_number(node->l->l) *
                              parser_get_number(node->l->r));
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_F2 &&
                 ((struct parser_f2*)(node->r))->ftype == PARSER_POW &&
                 parser_node_equal(((struct parser_f2*)(node->r))->l, node->l))
        { // x * pow(x,n)
            auto* xtmp = node->l;
            auto* ptmp = node->r;
            std::memcpy(node, node->r, sizeof(struct parser_node));
            ptmp->type = PARSER_ADD;
            ptmp->l = ((struct parser_f2*)node)->r;
            ptmp->r = xtmp;
            parser_set_number(ptmp->r, 1.0);
            ((struct parser_f2*)node)->r = ptmp;
            parser_ast_optimize(((struct parser_f2*)node)->r,local_consts);
        }
        else if (node->l->type == PARSER_F2 &&
                 ((struct parser_f2*)(node->l))->ftype == PARSER_POW &&
                 parser_node_equal(((struct parser_f2*)(node->l))->l, node->r))
        { // pow(x,n) * x
            auto* xtmp = node->r;
            auto* ptmp = node->l;
            std::memcpy(node, node->l, sizeof(struct parser_node));
            ptmp->type = PARSER_ADD;
            ptmp->l = ((struct parser_f2*)node)->r;
            ptmp->r = xtmp;
            parser_set_number(ptmp->r, 1.0);
            ((struct parser_f2*)node)->r = ptmp;
            parser_ast_optimize(((struct parser_f2*)node)->r,local_consts);
        }
        else if (node->r->type == PARSER_F2 &&
                 ((struct parser_f2*)(node->r))->ftype == PARSER_POW &&
                 node->l->type == PARSER_DIV &&
                 parser_node_equal(((struct parser_f2*)(node->r))->l, node->l->r))
        { // (a/x) * pow(x,n)
            std::swap(node->l, ((struct parser_f2*)(node->r))->r);
            std::swap(node->l, ((struct parser_f2*)(node->r))->r->l);
            ((struct parser_f2*)(node->r))->r->type = PARSER_ADD;
            parser_set_number(((struct parser_f2*)(node->r))->r->r, -1.0);
            parser_ast_optimize(node,local_consts);
        }
        // At this point, we have handled all directdly combinable cases.
        else if (group_combinables(node->l, node->r, is_mul_combinable, PARSER_MUL))
        {
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type != PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER)
        { // x * (3*y) = 3 * (x*y)  // NOLINT(bugprone-branch-clone)
            std::swap(node->l, node->r->l);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type != PARSER_NUMBER &&
                 node->r->type == PARSER_DIV &&
                 node->r->l->type == PARSER_NUMBER)
        { // x * (3/y) = 3 * (x/y)  // NOLINT(bugprone-branch-clone)
            std::swap(node->l, node->r->l);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD &&
                 node->r->l->type == PARSER_NUMBER)
        { // 3 * (4 + x) => 12 + 3*x
            std::swap(node->l, node->r->l);
            parser_set_number(node->l, parser_get_number(node->l) *
                              parser_get_number(node->r->l));
            node->type = PARSER_ADD;
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD &&
                 node->r->l->type == PARSER_MUL &&
                 node->r->l->l->type == PARSER_NUMBER)
        { // 3 * (4*x + y) => 12*x + 3*y
            std::swap(node->l, node->r->l); // (4*x) * (3 + y)
            parser_set_number(node->l->l, parser_get_number(node->l->l) *
                              parser_get_number(node->r->l)); // (12*x) * (3+y)
            node->type = PARSER_ADD;
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD &&
                 node->r->r->type == PARSER_MUL &&
                 node->r->r->l->type == PARSER_NUMBER)
        { // 3 * (x + 4*y) => 3*x + 12*y
            std::swap(node->l, node->r->r); // (4*y) * (x+3)
            parser_set_number(node->l->l, parser_get_number(node->l->l) *
                              parser_get_number(node->r->r)); // (12*y) * (x+3)
            node->type = PARSER_ADD;
            node->r->type = PARSER_MUL;
            std::swap(node->r->l, node->r->r);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER)
        { // (4*x) * y => 4*(x*y)  // NOLINT(bugprone-branch-clone)
            std::swap(node->l->l, node->r);
            std::swap(node->l, node->r);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_DIV &&
                 node->l->l->type == PARSER_NUMBER)
        { // (4/x) * y => 4*(y/x)  // NOLINT(bugprone-branch-clone)
            std::swap(node->l->l, node->r);
            std::swap(node->l, node->r);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_DIV &&
                 node->r->type == PARSER_DIV)
        { // (x/y) * (a/b) => (x*a)/(y*b)
            std::swap(node->l->r, node->r->l);
            node->type = PARSER_DIV;
            node->l->type = PARSER_MUL;
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_F2 &&
                 node->r->type == PARSER_F2 &&
                 ((struct parser_f2*)node->l)->ftype == PARSER_POW &&
                 ((struct parser_f2*)node->r)->ftype == PARSER_POW &&
                 parser_node_equal(((struct parser_f2*)(node->l))->l,
                                   ((struct parser_f2*)(node->r))->l))
        { // pow(x^m) * pow(x^n) => pow(x^(m+n))
            auto* l = (struct parser_f2*)(node->l);
            auto* r = (struct parser_f2*)(node->r);
            std::swap(l->r, r->l);
            std::swap(l->r, node->r);
            node->l->r->type = PARSER_ADD;
            std::memcpy(node, node->l, sizeof(struct parser_node));
            parser_ast_optimize(node,local_consts);
        }
        else if (try_cancel(node->l, -1, node->r))
        { // x * (.../(...x...))
            parser_ast_optimize(node,local_consts);
        }
        else if (try_cancel(node->r, -1, node->l))
        { // (.../(...x...)) * x
            parser_ast_optimize(node,local_consts);
        }
        break;
    case PARSER_DIV:
        parser_ast_optimize(node->l,local_consts);
        parser_ast_optimize(node->r,local_consts);
        if (node->l->type == PARSER_NUMBER &&
            parser_get_number(node->l) == 0.0)
        {
            parser_set_number(node, 0.0);
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_NUMBER)
        {
            parser_set_number(node, parser_get_number(node->l) /
                              parser_get_number(node->r));
        }
        else if (parser_node_equal(node->l, node->r))
        {
            parser_set_number(node, 1.0);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_NUMBER)
        { // (4*x)/3 => (4/3) * x
            std::swap(node->l, node->r);
            parser_set_number(node->l, parser_get_number(node->r->l) /
                              parser_get_number(node->l));
            node->r = node->r->r;
            node->type = PARSER_MUL;
        }
        else if (node->l->type == PARSER_DIV &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_NUMBER)
        { // (4/x)/3 => (4/3) / x
            std::swap(node->l->r, node->r);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER)
        { // (4*x)/(3*y) => (4/3) * (x/y)
            std::swap(node->l->r, node->r->l);
            node->type = PARSER_MUL;
            node->l->type = PARSER_DIV;
            node->r->type = PARSER_DIV;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_DIV &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER)
        { // (4/x)/(3*y) => (4/3) / (x*y)
            std::swap(node->l->r, node->r->l);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_MUL &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV &&
                 node->r->l->type == PARSER_NUMBER)
        { // (4*x)/(3/y) => (4/3)*(x*y)
            std::swap(node->l->r, node->r->l);
            node->l->type = PARSER_DIV;
            node->type = PARSER_MUL;
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_DIV &&
                 node->l->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV &&
                 node->r->l->type == PARSER_NUMBER)
        { // (4/x)/(3/y) => (4/3) * (y/x)
            std::swap(node->l->r, node->r->l);
            node->type = PARSER_MUL;
            std::swap(node->r->l, node->r->r);
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_NUMBER)
        { // x / 3 => (1/3) * x
            std::swap(node->l, node->r);
            parser_set_number(node->l, 1./parser_get_number(node->l));
            node->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_MUL &&
                 node->r->l->type == PARSER_NUMBER)
        { // x / (3*y) => (1./3)*(x/y)
            std::swap(node->l, node->r->l);
            parser_set_number(node->l, 1.0/parser_get_number(node->l));
            node->type = PARSER_MUL;
            node->r->type = PARSER_DIV;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_DIV &&
                 node->r->l->type == PARSER_NUMBER)
        { // x / (3/y) => (1./3.)*(x*y)
            std::swap(node->l, node->r->l);
            parser_set_number(node->l, 1.0/parser_get_number(node->l));
            node->type = PARSER_MUL;
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_DIV &&
                 node->l->type == PARSER_DIV)
        { // (x/y)/(a/b) => (b*x) / (a*y)
            std::swap(node->l->r, node->r->r); // (x/b) / (a/y)
            node->l->type = PARSER_MUL;
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_DIV)
        { // x / (y/z) => (x*z)/y     // NOLINT(bugprone-branch-clone)
            std::swap(node->l, node->r->l); // y/(x/z)
            std::swap(node->l, node->r);
            node->l->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_DIV)
        { // (x/y) / z => x/(y*z)     // NOLINT(bugprone-branch-clone)
            std::swap(node->l, node->r); // z / (x/y)
            std::swap(node->l, node->r->l); // x / (z/y)
            node->r->type = PARSER_MUL;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->l->type == PARSER_F2 &&
                 node->r->type == PARSER_F2 &&
                 ((struct parser_f2*)(node->l))->ftype == PARSER_POW &&
                 ((struct parser_f2*)(node->r))->ftype == PARSER_POW &&
                 parser_node_equal(((struct parser_f2*)(node->l))->l,
                                   ((struct parser_f2*)(node->r))->l))
        { // pow(x,m) / pow(x,n) => pow(x,(m-n))
            auto* l = node->l;
            auto* r = node->r;
            std::memcpy(node, node->l, sizeof(struct parser_node));
            std::swap(l->l, l->r);
            l->type = PARSER_ADD;
            l->r = r;
            l->r->type = PARSER_MUL;
            parser_set_number(l->r->l, -1.0);
            node->r = l;
            parser_ast_optimize(node,local_consts);
        }
        else if (node->r->type == PARSER_F2 &&
                 ((struct parser_f2*)(node->r))->ftype == PARSER_POW &&
                 ((struct parser_f2*)(node->r))->r->type == PARSER_NUMBER)
        { // f(.) / pow(x,n) => f(.) * pow(x,-n)
            node->type = PARSER_MUL;
            parser_set_number(((struct parser_f2*)(node->r))->r,
                              -parser_get_number(((struct parser_f2*)(node->r))->r));
            parser_ast_optimize(node,local_consts);
        }
        else if (try_cancel(node->l, 1, node->r))
        { // x / (...x...)
            parser_ast_optimize(node,local_consts);
        }
        else if (try_cancel(node->r, 1, node->l))
        { // (...x...) / x
            parser_ast_optimize(node,local_consts);
        }
        break;
    case PARSER_F1:
        parser_ast_optimize(((struct parser_f1*)node)->l,local_consts);
        if (((struct parser_f1*)node)->l->type == PARSER_NUMBER)
        {
            double v = parser_call_f1
                (((struct parser_f1*)node)->ftype,
                 ((struct parser_number*)(((struct parser_f1*)node)->l))->value);
            parser_set_number(node, v);
        }
        break;
    case PARSER_F2:
        parser_ast_optimize(((struct parser_f2*)node)->l,local_consts);
        parser_ast_optimize(((struct parser_f2*)node)->r,local_consts);
        if (((struct parser_f2*)node)->ftype == PARSER_CMP_CHAIN) {
            ((struct parser_f2*)node)->ftype = PARSER_AND;
        }
        if (((struct parser_f2*)node)->l->type == PARSER_NUMBER &&
            ((struct parser_f2*)node)->r->type == PARSER_NUMBER)
        {
            double v = parser_call_f2
                (((struct parser_f2*)node)->ftype,
                 ((struct parser_number*)(((struct parser_f2*)node)->l))->value,
                 ((struct parser_number*)(((struct parser_f2*)node)->r))->value);
            parser_set_number(node, v);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_AND &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) == 0.0)
        { // ? and false => false
            parser_set_number(node, 0.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_AND &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) != 0.0)
        { // ? and true => ?
            std::memcpy(node, node->l, sizeof(struct parser_node));
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_AND &&
                 ((struct parser_f2*)node)->l->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->l) == 0.0)
        { // false and ? => false
            parser_set_number(node, 0.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_AND &&
                 ((struct parser_f2*)node)->l->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->l) != 0.0)
        { // true and ? => ?
            std::memcpy(node, node->r, sizeof(struct parser_node));
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_OR &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) != 0.0)
        { // ? or true => true
            parser_set_number(node, 1.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_OR &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) == 0.0)
        { // ? or false => ?
            std::memcpy(node, node->l, sizeof(struct parser_node));
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_OR &&
                 ((struct parser_f2*)node)->l->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->l) != 0.0)
        { // true or ? => true
            parser_set_number(node, 1.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_OR &&
                 ((struct parser_f2*)node)->l->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->l) == 0.0)
        { // false or ? => ?
            std::memcpy(node, node->r, sizeof(struct parser_node));
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_POW &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) == 0.0)
        {
            parser_set_number(node, 1.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_POW &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) == 1.0)
        {
            std::memcpy(node, ((struct parser_f2*)node)->l,
                        sizeof(struct parser_node));
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_POW &&
                 ((struct parser_f2*)node)->l->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->l) == 0.0)
        {
            parser_set_number(node, 0.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_POW &&
                 ((struct parser_f2*)node)->r->type == PARSER_NUMBER &&
                 parser_get_number(((struct parser_f2*)node)->r) == -1.0)
        {
            std::swap(node->l, node->r);
            node->type = PARSER_DIV;
            parser_set_number(node->l, 1.0);
        }
        else if (((struct parser_f2*)node)->ftype == PARSER_POW &&
                 ((struct parser_f2*)node)->l->type == PARSER_F2 &&
                 ((struct parser_f2*)((struct parser_f2*)node)->l)->ftype == PARSER_POW)
        { // pow(pow(,),)
            std::swap(node->l, node->r);
            std::swap(node->l, node->r->l);
            node->r->type = PARSER_MUL;
            parser_ast_sort(node);
            parser_ast_optimize(node, local_consts);
        }
        break;
    case PARSER_F3:
        parser_ast_optimize(((struct parser_f3*)node)->n1,local_consts);
        parser_ast_optimize(((struct parser_f3*)node)->n2,local_consts);
        parser_ast_optimize(((struct parser_f3*)node)->n3,local_consts);
        if (((struct parser_f3*)node)->n1->type == PARSER_NUMBER &&
            ((struct parser_f3*)node)->n2->type == PARSER_NUMBER &&
            ((struct parser_f3*)node)->n3->type == PARSER_NUMBER)
        {
            double v = parser_call_f3
                (((struct parser_f3*)node)->ftype,
                 ((struct parser_number*)(((struct parser_f3*)node)->n1))->value,
                 ((struct parser_number*)(((struct parser_f3*)node)->n2))->value,
                 ((struct parser_number*)(((struct parser_f3*)node)->n3))->value);
            parser_set_number(node, v);
        }
        else if (((struct parser_f3*)node)->n1->type == PARSER_NUMBER &&
                 ((struct parser_f3*)node)->ftype == PARSER_IF)
        {
            if (parser_get_number(((struct parser_f3*)node)->n1) == 0.0) {
                std::memcpy(node, ((struct parser_f3*)node)->n3,
                            sizeof(struct parser_node));
            } else {
                std::memcpy(node, ((struct parser_f3*)node)->n2,
                            sizeof(struct parser_node));
            }
        }
        break;
    case PARSER_USRF1:
        parser_ast_optimize(((struct parser_usrf1*)node)->l,local_consts);
        break;
    case PARSER_USRF2:
        parser_ast_optimize(((struct parser_usrf2*)node)->l,local_consts);
        parser_ast_optimize(((struct parser_usrf2*)node)->r,local_consts);
        break;
    case PARSER_USRFN:
    {
        parser_ast_optimize(((struct parser_usrfn*)node)->n1,local_consts);
        short argc = ((struct parser_usrfn*)node)->argc;
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_ast_optimize(((struct parser_usrfn*)node)->others[iarg],local_consts);
        }
        break;
    }
    case PARSER_ASSIGN:
    {   // lhs = rhs
        parser_ast_optimize(((struct parser_assign*)node)->v,local_consts);
        std::string lhs(((struct parser_assign*)node)->s->name);
        if (auto it = local_consts.find(lhs); it != local_consts.end()) {
            // lhs exists as a local constant
            if (((struct parser_assign*)node)->v->type == PARSER_NUMBER) {
                // Update its value. It's okay if it's the same value.
                // Then delete this node.
                it->second = parser_get_number(((struct parser_assign*)node)->v);
                node = nullptr;
            } else {
                // No longer a constant.
                local_consts.erase(it);
            }
        } else {
            // If lhs is not an existing constant and rhs is a number, we
            // add it to the local constant map and delete this node.
            if (((struct parser_assign*)node)->v->type == PARSER_NUMBER) {

                local_consts[lhs] = parser_get_number(((struct parser_assign*)node)->v);
                node = nullptr;
            }
        }
        break;
    }
    case PARSER_LIST:
        parser_ast_optimize(node->l,local_consts);
        parser_ast_optimize(node->r,local_consts);
        if (node->l == nullptr && node->r == nullptr) {
            node = nullptr;
        } else if (node->l && node->r == nullptr) {
            std::memcpy(node, node->l, sizeof(struct parser_node));
        } else if (node->r && node->l == nullptr) {
            std::memcpy(node, node->r, sizeof(struct parser_node));
        }
        break;
    case PARSER_SUB:
        amrex::Abort("parser_ast_optimize: should not have PARSER_SUB");
        break;
    default:
        amrex::Abort("parser_ast_optimize: unknown node type " +
                     std::to_string(node->type));
    }
}

namespace {

void
parser_ast_print_f1 (struct parser_f1* f1, std::string const& space, std::ostream& printer)
{
    printer << space << parser_f1_s[f1->ftype] << "\n";
    parser_ast_print(f1->l, space+"  ", printer);
}

void
parser_ast_print_f2 (struct parser_f2* f2, std::string const& space, std::ostream& printer)
{
    printer << space << parser_f2_s[f2->ftype] << "\n";
    parser_ast_print(f2->l, space+"  ", printer);
    parser_ast_print(f2->r, space+"  ", printer);
}

void
parser_ast_print_f3 (struct parser_f3* f3, std::string const& space, std::ostream& printer)
{
    printer << space << parser_f3_s[f3->ftype] << "\n";
    std::string const& more_space = space + "  ";
    parser_ast_print(f3->n1, more_space, printer);
    parser_ast_print(f3->n2, more_space, printer);
    parser_ast_print(f3->n3, more_space, printer);
}

void
parser_ast_print_usrf1 (struct parser_usrf1* f1, std::string const& space, std::ostream& printer)
{
    printer << space << f1->name << "\n";
    parser_ast_print(f1->l, space+"  ", printer);
}

void
parser_ast_print_usrf2 (struct parser_usrf2* f2, std::string const& space, std::ostream& printer)
{
    printer << space << f2->name << "\n";
    parser_ast_print(f2->l, space+"  ", printer);
    parser_ast_print(f2->r, space+"  ", printer);
}

void
parser_ast_print_usrfn (struct parser_usrfn* fn, std::string const& space, std::ostream& printer)
{
    printer << space << fn->name << "\n";
    std::string const& more_space = space + "  ";
    parser_ast_print(fn->n1, more_space, printer);
    for (short iarg = 0; iarg < fn->argc-1; ++iarg) {
        parser_ast_print(fn->others[iarg], more_space, printer);
    }
}

}

void
parser_ast_print (struct parser_node* node, std::string const& space, std::ostream& printer)
{
    std::string const& more_space = space + "  ";

    switch (node->type)
    {
    case PARSER_NUMBER:
        printer << space << parser_node_s[node->type] << ": "
                << ((struct parser_number*)node)->value << "\n";
        break;
    case PARSER_SYMBOL:
        printer << space << parser_node_s[node->type] << ": "
                << ((struct parser_symbol*)node)->name << "\n";
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        printer << space << parser_node_s[node->type] << "\n";
        parser_ast_print(node->l, more_space, printer);
        parser_ast_print(node->r, more_space, printer);
        break;
    case PARSER_F1:
        parser_ast_print_f1((struct parser_f1*)node, space, printer);
        break;
    case PARSER_F2:
        parser_ast_print_f2((struct parser_f2*)node, space, printer);
        break;
    case PARSER_F3:
        parser_ast_print_f3((struct parser_f3*)node, space, printer);
        break;
    case PARSER_USRF1:
        parser_ast_print_usrf1((struct parser_usrf1*)node, space, printer);
        break;
    case PARSER_USRF2:
        parser_ast_print_usrf2((struct parser_usrf2*)node, space, printer);
        break;
    case PARSER_USRFN:
        parser_ast_print_usrfn((struct parser_usrfn*)node, space, printer);
        break;
    case PARSER_ASSIGN:
        printer << space << "=: " << ((struct parser_assign*)node)->s->name << " =\n";
        parser_ast_print(((struct parser_assign*)node)->v, more_space, printer);
        break;
    default:
        amrex::Abort("parser_ast_print: unknown node type " + std::to_string(node->type));
    }
}

int
parser_ast_depth (struct parser_node* node)
{
    if (node == nullptr) { return 0; }

    switch (node->type)
    {
    case PARSER_NUMBER:
    case PARSER_SYMBOL:
        return 1;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
    {
        int d1 = parser_ast_depth(node->l);
        int d2 = parser_ast_depth(node->r);
        return std::max(d1,d2)+1;
    }
    case PARSER_F1:
        return parser_ast_depth(((struct parser_f1*)node)->l) + 1;
    case PARSER_F2:
    {
        int d1 = parser_ast_depth(((struct parser_f2*)node)->l);
        int d2 = parser_ast_depth(((struct parser_f2*)node)->r);
        return std::max(d1,d2)+1;
    }
    case PARSER_F3:
    {
        int d1 = parser_ast_depth(((struct parser_f3*)node)->n1);
        int d2 = parser_ast_depth(((struct parser_f3*)node)->n2);
        int d3 = parser_ast_depth(((struct parser_f3*)node)->n3);
        return std::max({d1,d2,d3})+1;
    }
    case PARSER_USRF1:
        return parser_ast_depth(((struct parser_usrf1*)node)->l) + 1;
    case PARSER_USRF2:
    {
        int d1 = parser_ast_depth(((struct parser_usrf2*)node)->l);
        int d2 = parser_ast_depth(((struct parser_usrf2*)node)->r);
        return std::max(d1,d2)+1;
    }
    case PARSER_USRFN:
    {
        int d1 = parser_ast_depth(((struct parser_usrfn*)node)->n1);
        short argc = ((struct parser_usrfn*)node)->argc;
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            int dn = parser_ast_depth(((struct parser_usrfn*)node)->others[iarg]);
            d1 = std::max(d1,dn);
        }
        return d1+1;
    }
    case PARSER_ASSIGN:
    {
        int d = parser_ast_depth(((struct parser_assign*)node)->v);
        return d+1;
    }
    default:
        amrex::Abort("parser_ast_depth: unknown node type " + std::to_string(node->type));
        return 0;
    }
}

void parser_ast_sort (struct parser_node* node)
{
    if (node == nullptr) { return; }

    switch (node->type)
    {
    case PARSER_NUMBER:
    case PARSER_SYMBOL:
        break;
    case PARSER_ADD:
    case PARSER_MUL:
    {
        parser_ast_sort(node->l);
        parser_ast_sort(node->r);
        if (parser_node_compare(node->r, node->l)) {
            std::swap(node->l, node->r);
        }
        break;
    }
    case PARSER_SUB:
    case PARSER_DIV:
    case PARSER_LIST:
        parser_ast_sort(node->l);
        parser_ast_sort(node->r);
        break;
    case PARSER_F1:
        parser_ast_sort(((struct parser_f1*)node)->l);
        break;
    case PARSER_F2:
        parser_ast_sort(((struct parser_f2*)node)->l);
        parser_ast_sort(((struct parser_f2*)node)->r);
        break;
    case PARSER_F3:
        parser_ast_sort(((struct parser_f3*)node)->n1);
        parser_ast_sort(((struct parser_f3*)node)->n2);
        parser_ast_sort(((struct parser_f3*)node)->n3);
        break;
    case PARSER_USRF1:
        parser_ast_sort(((struct parser_usrf1*)node)->l);
        break;
    case PARSER_USRF2:
        parser_ast_sort(((struct parser_usrf2*)node)->l);
        parser_ast_sort(((struct parser_usrf2*)node)->r);
        break;
    case PARSER_USRFN:
    {
        parser_ast_sort(((struct parser_usrfn*)node)->n1);
        short argc = ((struct parser_usrfn*)node)->argc;
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_ast_sort(((struct parser_usrfn*)node)->others[iarg]);
        }
        break;
    }
    case PARSER_ASSIGN:
        parser_ast_sort(((struct parser_assign*)node)->v);
        break;
    default:
        amrex::Abort("parser_ast_sort: unknown node type " + std::to_string(node->type));
    }
}

void
parser_ast_regvar (struct parser_node* node, char const* name, int i)
{
    switch (node->type)
    {
    case PARSER_NUMBER:
        break;
    case PARSER_SYMBOL:
        if (std::strcmp(name, ((struct parser_symbol*)node)->name) == 0) {
            ((struct parser_symbol*)node)->ip = i;
        }
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        parser_ast_regvar(node->l, name, i);
        parser_ast_regvar(node->r, name, i);
        break;
    case PARSER_F1:
        parser_ast_regvar(((struct parser_f1*)node)->l, name, i);
        break;
    case PARSER_F2:
        parser_ast_regvar(((struct parser_f2*)node)->l, name, i);
        parser_ast_regvar(((struct parser_f2*)node)->r, name, i);
        break;
    case PARSER_F3:
        parser_ast_regvar(((struct parser_f3*)node)->n1, name, i);
        parser_ast_regvar(((struct parser_f3*)node)->n2, name, i);
        parser_ast_regvar(((struct parser_f3*)node)->n3, name, i);
        break;
    case PARSER_USRF1:
        parser_ast_regvar(((struct parser_usrf1*)node)->l, name, i);
        break;
    case PARSER_USRF2:
        parser_ast_regvar(((struct parser_usrf2*)node)->l, name, i);
        parser_ast_regvar(((struct parser_usrf2*)node)->r, name, i);
        break;
    case PARSER_USRFN:
    {
        parser_ast_regvar(((struct parser_usrfn*)node)->n1, name, i);
        short argc = ((struct parser_usrfn*)node)->argc;
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_ast_regvar(((struct parser_usrfn*)node)->others[iarg], name, i);
        }
        break;
    }
    case PARSER_ASSIGN:
        parser_ast_regvar(((struct parser_assign*)node)->v, name, i);
        break;
    default:
        amrex::Abort("parser_ast_regvar: unknown node type "+std::to_string(node->type));
    }
}

void parser_ast_setconst (struct parser_node* node, char const* name, double c)
{
    switch (node->type)
    {
    case PARSER_NUMBER:
        break;
    case PARSER_SYMBOL:
        if (std::strcmp(name, ((struct parser_symbol*)node)->name) == 0) {
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = c;
        }
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        parser_ast_setconst(node->l, name, c);
        parser_ast_setconst(node->r, name, c);
        break;
    case PARSER_F1:
        parser_ast_setconst(((struct parser_f1*)node)->l, name, c);
        break;
    case PARSER_F2:
        parser_ast_setconst(((struct parser_f2*)node)->l, name, c);
        parser_ast_setconst(((struct parser_f2*)node)->r, name, c);
        break;
    case PARSER_F3:
        parser_ast_setconst(((struct parser_f3*)node)->n1, name, c);
        parser_ast_setconst(((struct parser_f3*)node)->n2, name, c);
        parser_ast_setconst(((struct parser_f3*)node)->n3, name, c);
        break;
    case PARSER_USRF1:
        parser_ast_setconst(((struct parser_usrf1*)node)->l, name, c);
        break;
    case PARSER_USRF2:
        parser_ast_setconst(((struct parser_usrf2*)node)->l, name, c);
        parser_ast_setconst(((struct parser_usrf2*)node)->r, name, c);
        break;
    case PARSER_USRFN:
    {
        parser_ast_setconst(((struct parser_usrfn*)node)->n1, name, c);
        short argc = ((struct parser_usrfn*)node)->argc;
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_ast_setconst(((struct parser_usrfn*)node)->others[iarg], name, c);
        }
        break;
    }
    case PARSER_ASSIGN:
        parser_ast_setconst(((struct parser_assign*)node)->v, name, c);
        break;
    default:
        amrex::Abort("parser_ast_setconst: unknown node type " + std::to_string(node->type));
    }
}

void parser_ast_get_symbols (struct parser_node* node, std::set<std::string>& symbols,
                             std::set<std::string>& local_symbols)
{
    switch (node->type)
    {
    case PARSER_NUMBER:
        break;
    case PARSER_SYMBOL:
        symbols.emplace(((struct parser_symbol*)node)->name);
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        parser_ast_get_symbols(node->l, symbols, local_symbols);
        parser_ast_get_symbols(node->r, symbols, local_symbols);
        break;
    case PARSER_F1:
        parser_ast_get_symbols(((struct parser_f1*)node)->l, symbols, local_symbols);
        break;
    case PARSER_F2:
        parser_ast_get_symbols(((struct parser_f2*)node)->l, symbols, local_symbols);
        parser_ast_get_symbols(((struct parser_f2*)node)->r, symbols, local_symbols);
        break;
    case PARSER_F3:
        parser_ast_get_symbols(((struct parser_f3*)node)->n1, symbols, local_symbols);
        parser_ast_get_symbols(((struct parser_f3*)node)->n2, symbols, local_symbols);
        parser_ast_get_symbols(((struct parser_f3*)node)->n3, symbols, local_symbols);
        break;
    case PARSER_USRF1:
        parser_ast_get_symbols(((struct parser_usrf1*)node)->l, symbols, local_symbols);
        break;
    case PARSER_USRF2:
        parser_ast_get_symbols(((struct parser_usrf2*)node)->l, symbols, local_symbols);
        parser_ast_get_symbols(((struct parser_usrf2*)node)->r, symbols, local_symbols);
        break;
    case PARSER_USRFN:
    {
        parser_ast_get_symbols(((struct parser_usrfn*)node)->n1, symbols, local_symbols);
        short argc = ((struct parser_usrfn*)node)->argc;
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_ast_get_symbols(((struct parser_usrfn*)node)->others[iarg], symbols, local_symbols);
        }
        break;
    }
    case PARSER_ASSIGN:
        local_symbols.emplace(((struct parser_assign*)node)->s->name);
        parser_ast_get_symbols(((struct parser_assign*)node)->v, symbols, local_symbols);
        break;
    default:
        amrex::Abort("parser_ast_get_symbols: unknown node type " + std::to_string(node->type));
    }
}

void parser_ast_get_user_functions (struct parser_node* node,
                                    std::map<std::string,int>& user_functions)
{
    switch (node->type)
    {
    case PARSER_NUMBER:
        break;
    case PARSER_SYMBOL:
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_LIST:
        parser_ast_get_user_functions(node->l, user_functions);
        parser_ast_get_user_functions(node->r, user_functions);
        break;
    case PARSER_F1:
        parser_ast_get_user_functions(((struct parser_f1*)node)->l, user_functions);
        break;
    case PARSER_F2:
        parser_ast_get_user_functions(((struct parser_f2*)node)->l, user_functions);
        parser_ast_get_user_functions(((struct parser_f2*)node)->r, user_functions);
        break;
    case PARSER_F3:
        parser_ast_get_user_functions(((struct parser_f3*)node)->n1, user_functions);
        parser_ast_get_user_functions(((struct parser_f3*)node)->n2, user_functions);
        parser_ast_get_user_functions(((struct parser_f3*)node)->n3, user_functions);
        break;
    case PARSER_USRF1:
    {
        std::string fname(((struct parser_usrf1*)node)->name);
        auto it = user_functions.find(fname);
        if (it != user_functions.end()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it->second == 1,
                "amrex::Parser: user function overloading not supported");
        } else {
            user_functions[fname] = 1;
        }
        parser_ast_get_user_functions(((struct parser_usrf1*)node)->l, user_functions);
        break;
    }
    case PARSER_USRF2:
    {
        std::string fname(((struct parser_usrf2*)node)->name);
        auto it = user_functions.find(fname);
        if (it != user_functions.end()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it->second == 2,
                "amrex::Parser: user function overloading not supported");
        } else {
            user_functions[fname] = 2;
        }
        parser_ast_get_user_functions(((struct parser_usrf2*)node)->l, user_functions);
        parser_ast_get_user_functions(((struct parser_usrf2*)node)->r, user_functions);
        break;
    }
    case PARSER_USRFN:
    {
        short argc = ((struct parser_usrfn*)node)->argc;
        std::string fname(((struct parser_usrfn*)node)->name);
        auto it = user_functions.find(fname);
        if (it != user_functions.end()) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it->second == argc,
                "amrex::Parser: user function overloading not supported");
        } else {
            user_functions[fname] = argc;
        }
        parser_ast_get_user_functions(((struct parser_usrfn*)node)->n1, user_functions);
        for (short iarg = 0; iarg < argc-1; ++iarg) {
            parser_ast_get_user_functions(((struct parser_usrfn*)node)->others[iarg], user_functions);
        }
        break;
    }
    case PARSER_ASSIGN:
        parser_ast_get_user_functions(((struct parser_assign*)node)->v, user_functions);
        break;
    default:
        amrex::Abort("parser_ast_get_user_functions: unknown node type " + std::to_string(node->type));
    }
}

void
parser_regvar (struct amrex_parser* parser, char const* name, int i)
{
    parser_ast_regvar(parser->ast, name, i);
}

void
parser_setconst (struct amrex_parser* parser, char const* name, double c)
{
    parser_ast_setconst(parser->ast, name, c);
    std::map<std::string,double> local_consts;
    parser_ast_optimize(parser->ast, local_consts);
    parser_ast_sort(parser->ast);
}

void
parser_print (struct amrex_parser* parser)
{
    auto& printer = amrex::OutStream();
    auto oldprec = printer.precision(17);
    parser_ast_print(parser->ast, std::string("  "), printer);
    printer.precision(oldprec);
}

std::set<std::string>
parser_get_symbols (struct amrex_parser* parser)
{
    std::set<std::string> symbols;
    std::set<std::string> local_symbols;
    parser_ast_get_symbols(parser->ast, symbols, local_symbols);
    for (auto const& ls : local_symbols) {
        symbols.erase(ls);
    }
    return symbols;
}

std::map<std::string,int>
parser_get_user_functions (struct amrex_parser* parser)
{
    std::map<std::string,int> user_functions;
    parser_ast_get_user_functions(parser->ast, user_functions);
    return user_functions;
}

int
parser_depth (struct amrex_parser* parser)
{
    return parser_ast_depth(parser->ast);
}

double parser_get_number (struct parser_node* node)
{
    AMREX_ASSERT(node->type == PARSER_NUMBER);
    return ((struct parser_number*)node)->value;
}

void parser_set_number (struct parser_node* node, double v)
{
    ((struct parser_number*)node)->value = v;
    node->type = PARSER_NUMBER;
}

}
