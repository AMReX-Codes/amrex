#include <AMReX.H>
#include <AMReX_Parser_Y.H>
#include <amrex_parser.tab.h>

#include <algorithm>
#include <cstdarg>
#include <string>

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

static struct parser_node* parser_root = nullptr;

// This is called by a bison rule to store the original AST in a static variable.
void
parser_defexpr (struct parser_node* body)
{
    parser_root = body;
}

struct parser_symbol*
parser_makesymbol (char* name)
{
    auto *symbol = (struct parser_symbol*) std::malloc(sizeof(struct parser_symbol));
    symbol->type = PARSER_SYMBOL;
    symbol->name = strdup(name);
    symbol->ip = -1;
    return symbol;
}

struct parser_node*
parser_newnode (enum parser_node_t type, struct parser_node* l, struct parser_node* r)
{
    auto *tmp = (struct parser_node*) std::malloc(sizeof(struct parser_node));
    tmp->type = type;
    tmp->l = l;
    tmp->r = r;
    return tmp;
}

struct parser_node*
parser_newnumber (double d)
{
    auto *r = (struct parser_number*) std::malloc(sizeof(struct parser_number));
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
    auto *tmp = (struct parser_f1*) std::malloc(sizeof(struct parser_f1));
    tmp->type = PARSER_F1;
    tmp->l = l;
    tmp->ftype = ftype;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newf2 (enum parser_f2_t ftype, struct parser_node* l, struct parser_node* r)
{
    auto *tmp = (struct parser_f2*) std::malloc(sizeof(struct parser_f2));
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
    auto *tmp = (struct parser_f3*) std::malloc(sizeof(struct parser_f3));
    tmp->type = PARSER_F3;
    tmp->n1 = n1;
    tmp->n2 = n2;
    tmp->n3 = n3;
    tmp->ftype = ftype;
    return (struct parser_node*) tmp;
}

struct parser_node*
parser_newassign (struct parser_symbol* sym, struct parser_node* v)
{
    auto *r = (struct parser_assign*) std::malloc(sizeof(struct parser_assign));
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
        auto *r = (struct parser_node*) std::malloc(sizeof(struct parser_node));
        r->type = PARSER_LIST;
        r->l = nl;
        r->r = nr;
        return r;
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

    my_parser->ast = parser_ast_dup(my_parser, parser_root, 1); /* 1: free the source parser_root */

    if ((char*)my_parser->p_root + my_parser->sz_mempool != (char*)my_parser->p_free) {
        amrex::Abort("amrex_parser_new: error in memory size");
    }

    parser_ast_optimize(my_parser->ast);

    return my_parser;
}

void
amrex_parser_delete (struct amrex_parser* parser)
{
    std::free(parser->p_root);
    std::free(parser);
}

static
std::size_t
parser_aligned_size (std::size_t N)
{
    const unsigned int align_size = 16;
    std::size_t x = N + (align_size-1);
    x -= x & (align_size-1);
    return x;
}

static
void*
parser_allocate (struct amrex_parser* my_parser, std::size_t N)
{
    void* r = my_parser->p_free;
    my_parser->p_free = (char*)r + parser_aligned_size(N);
    return r;
}

struct amrex_parser*
parser_dup (struct amrex_parser* source)
{
    auto *dest = (struct amrex_parser*) std::malloc(sizeof(struct amrex_parser));
    dest->sz_mempool = source->sz_mempool;
    dest->p_root = std::malloc(dest->sz_mempool);
    dest->p_free = dest->p_root;

    dest->ast = parser_ast_dup(dest, source->ast, 0); /* 0: don't free the source */

    return dest;
}

std::size_t
parser_ast_size (struct parser_node* node)
{
    std::size_t result = 0;

    switch (node->type)
    {
    case PARSER_NUMBER:
        result = parser_aligned_size(sizeof(struct parser_number));
        break;
    case PARSER_SYMBOL:
        result = parser_aligned_size(    sizeof(struct parser_symbol))
            + parser_aligned_size(std::strlen(((struct parser_symbol*)node)->name)+1);
        break;
    case PARSER_ADD:
    case PARSER_SUB:
    case PARSER_MUL:
    case PARSER_DIV:
    case PARSER_ADD_PP:
    case PARSER_SUB_PP:
    case PARSER_MUL_PP:
    case PARSER_DIV_PP:
    case PARSER_LIST:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_ast_size(node->l) + parser_ast_size(node->r);
        break;
    case PARSER_NEG:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_ast_size(node->l);
        break;
    case PARSER_F1:
        result = parser_aligned_size(sizeof(struct parser_f1))
            +             parser_ast_size(((struct parser_f1*)node)->l);
        break;
    case PARSER_F2:
        result = parser_aligned_size(sizeof(struct parser_f2))
            +             parser_ast_size(((struct parser_f2*)node)->l)
            +             parser_ast_size(((struct parser_f2*)node)->r);
        break;
    case PARSER_F3:
        result = parser_aligned_size(sizeof(struct parser_f3))
            +             parser_ast_size(((struct parser_f3*)node)->n1)
            +             parser_ast_size(((struct parser_f3*)node)->n2)
            +             parser_ast_size(((struct parser_f3*)node)->n3);
        break;
    case PARSER_ASSIGN:
        result += parser_aligned_size(sizeof(struct parser_assign))
            + parser_ast_size((struct parser_node*)(((struct parser_assign*)node)->s))
            + parser_ast_size(((struct parser_assign*)node)->v);
        break;
    case PARSER_ADD_VP:
    case PARSER_SUB_VP:
    case PARSER_MUL_VP:
    case PARSER_DIV_VP:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_ast_size(node->r);
        break;
    case PARSER_NEG_P:
        result = parser_aligned_size(sizeof(struct parser_node))
            + parser_ast_size(node->l);
        break;
    default:
        amrex::Abort("parser_ast_size: unknown node type " + std::to_string(node->type));
    }

    return result;
}

struct parser_node*
parser_ast_dup (struct amrex_parser* my_parser, struct parser_node* node, int move)
{
    void* result = nullptr;

    switch (node->type)
    {
    case PARSER_NUMBER:
        result = parser_allocate(my_parser, sizeof(struct parser_number));
        std::memcpy(result, node          , sizeof(struct parser_number));
        break;
    case PARSER_SYMBOL:
    {
        result = parser_allocate(my_parser, sizeof(struct parser_symbol));
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
    case PARSER_ADD_PP:
    case PARSER_SUB_PP:
    case PARSER_MUL_PP:
    case PARSER_DIV_PP:
    case PARSER_LIST:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_node));
        ((struct parser_node*)result)->l = parser_ast_dup(my_parser, node->l, move);
        ((struct parser_node*)result)->r = parser_ast_dup(my_parser, node->r, move);
        break;
    case PARSER_NEG:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_node));
        ((struct parser_node*)result)->l = parser_ast_dup(my_parser, node->l, move);
        ((struct parser_node*)result)->r = nullptr;
        break;
    case PARSER_F1:
        result = parser_allocate(my_parser, sizeof(struct parser_f1));
        std::memcpy(result, node          , sizeof(struct parser_f1));
        ((struct parser_f1*)result)->l = parser_ast_dup(my_parser,
                                                 ((struct parser_f1*)node)->l, move);
        break;
    case PARSER_F2:
        result = parser_allocate(my_parser, sizeof(struct parser_f2));
        std::memcpy(result, node          , sizeof(struct parser_f2));
        ((struct parser_f2*)result)->l = parser_ast_dup(my_parser,
                                                 ((struct parser_f2*)node)->l, move);
        ((struct parser_f2*)result)->r = parser_ast_dup(my_parser,
                                                 ((struct parser_f2*)node)->r, move);
        break;
    case PARSER_F3:
        result = parser_allocate(my_parser, sizeof(struct parser_f3));
        std::memcpy(result, node          , sizeof(struct parser_f3));
        ((struct parser_f3*)result)->n1 = parser_ast_dup(my_parser,
                                                 ((struct parser_f3*)node)->n1, move);
        ((struct parser_f3*)result)->n2 = parser_ast_dup(my_parser,
                                                 ((struct parser_f3*)node)->n2, move);
        ((struct parser_f3*)result)->n3 = parser_ast_dup(my_parser,
                                                 ((struct parser_f3*)node)->n3, move);
        break;
    case PARSER_ASSIGN:
        result = parser_allocate(my_parser, sizeof(struct parser_assign));
        std::memcpy(result, node          , sizeof(struct parser_assign));
        ((struct parser_assign*)result)->s = (struct parser_symbol*)
            parser_ast_dup(my_parser, (struct parser_node*)
                                                (((struct parser_assign*)node)->s), move);
        ((struct parser_assign*)result)->v = parser_ast_dup(my_parser,
                                                 ((struct parser_assign*)node)->v, move);
        break;
    case PARSER_ADD_VP:
    case PARSER_SUB_VP:
    case PARSER_MUL_VP:
    case PARSER_DIV_VP:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_node));
        ((struct parser_node*)result)->r = parser_ast_dup(my_parser, node->r, move);
        break;
    case PARSER_NEG_P:
        result = parser_allocate(my_parser, sizeof(struct parser_node));
        std::memcpy(result, node          , sizeof(struct parser_node));
        ((struct parser_node*)result)->l = parser_ast_dup(my_parser, node->l, move);
        break;
    default:
        amrex::Abort("parser_ast_dup: unknown node type " + std::to_string(node->type));
    }
    if (move) {
        /* Note that we only do this on the original AST.  We do not
         * need to call free for AST stored in amrex_parser because the
         * memory is not allocated with std::malloc directly.
         */
        if (node->type == PARSER_SYMBOL) {
            std::free(((struct parser_symbol*)node)->name);
        }
        std::free((void*)node);
    }
    return (struct parser_node*)result;
}

#define PARSER_MOVEUP_R(node, v) \
    struct parser_node* n = (node)->r->r; \
    int ip = (node)->r->rip; \
    (node)->r = n; \
    (node)->lvp.v = v; \
    (node)->rip   = ip;
#define PARSER_MOVEUP_L(node, v) \
    struct parser_node* n = (node)->l->r; \
    int ip = (node)->l->rip; \
    (node)->r = n; \
    (node)->lvp.v = v; \
    (node)->rip   = ip;
#define PARSER_EVAL_R(node) (node)->r->lvp.v
#define PARSER_EVAL_L(node) (node)->l->lvp.v

#define PARSER_NEG_MOVEUP(node) \
    (node)->r = (node)->l->r; \
    (node)->lvp.v = -(node)->l->lvp.v; \
    (node)->rip = (node)->l->rip;

void
parser_ast_optimize (struct parser_node* node)
{
    /* No need to free memory because we only call this on ASTs in
     * amrex_parser that are allocated from the memory pool.
     */
    switch (node->type)
    {
    case PARSER_NUMBER:
    case PARSER_SYMBOL:
        break;
    case PARSER_ADD:
    case PARSER_ADD_PP:
        parser_ast_optimize(node->l);
        parser_ast_optimize(node->r);
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_NUMBER)
        {
            double v = ((struct parser_number*)(node->l))->value
                +    ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.v = ((struct parser_number*)(node->l))->value;
            node->rip   = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_NUMBER)
        {
            node->lvp.v = ((struct parser_number*)(node->r))->value;
            node->rip   = ((struct parser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
            node->rip    = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_ADD_PP; // For *_PP, the names are stored in the l and r nodes.
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD_VP)
        {
            double v = ((struct parser_number*)(node->l))->value + PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SUB_VP)
        {
            double v = ((struct parser_number*)(node->l))->value + PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_SUB_VP;
        }
        else if (node->l->type == PARSER_ADD_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) + ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_SUB_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) + ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_SUB_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l + (#rl + node_rr) -> (#l + #rl) + node_rr, same type
                double v = ((struct parser_number*)(node->l))->value
                    + ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l + (node_rl + #rr) -> (#l + #rr) + node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    + ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SUB)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l + (#rl - node_rr) -> (#l + #rl) - node_rr, type change
                double v = ((struct parser_number*)(node->l))->value
                    + ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
                node->type = PARSER_SUB;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l + (node_rl - #rr) -> (#l - #rr) + node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    - ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_ADD &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll + node_lr) + #r -> nodel_lr + (#ll + #r), same type
                double v = ((struct parser_number*)(node->l->l))->value
                    + ((struct parser_number*)(node->r))->value;
                node->l = node->l->r;
                ((struct parser_number*)(node->r))->value = v;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll + #lr) + #r -> node_ll + (#lr + #r), same type
                double v = ((struct parser_number*)(node->l->r))->value
                    + ((struct parser_number*)(node->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
            }
        }
        else if (node->l->type == PARSER_SUB &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll - node_lr) + #r -> (#ll + #r) - node_lr, type change
                double v = ((struct parser_number*)(node->l->l))->value
                    + ((struct parser_number*)(node->r))->value;
                node->r = node->l->r;
                ((struct parser_number*)(node->l))->type = PARSER_NUMBER;
                ((struct parser_number*)(node->l))->value = v;
                node->type = PARSER_SUB;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll - #lr) + #r -> node_ll + (#r - #lr), same type
                double v = ((struct parser_number*)(node->r))->value
                    - ((struct parser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
            }
        }
        break;
    case PARSER_SUB:
    case PARSER_SUB_PP:
        parser_ast_optimize(node->l);
        parser_ast_optimize(node->r);
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_NUMBER)
        {
            double v = ((struct parser_number*)(node->l))->value
                -    ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.v = ((struct parser_number*)(node->l))->value;
            node->rip   = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_SUB_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_NUMBER)
        {
            node->lvp.v = -((struct parser_number*)(node->r))->value;
            node->rip   =  ((struct parser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
            node->rip    = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_SUB_PP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD_VP)
        {
            double v = ((struct parser_number*)(node->l))->value - PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_SUB_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SUB_VP)
        {
            double v = ((struct parser_number*)(node->l))->value - PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_ADD_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) - ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_SUB_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) - ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_SUB_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_ADD)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l - (#rl + node_rr) -> (#l - #rl) - node_rr, same type
                double v = ((struct parser_number*)(node->l))->value
                    - ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l - (node_rl + #rr) -> (#l - #rr) - node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    - ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SUB)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l - (#rl - node_rr) -> (#l - #rl) + node_rr, type change
                double v = ((struct parser_number*)(node->l))->value
                    - ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
                node->type = PARSER_ADD;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l - (node_rl - #rr) -> (#l + #rr) - node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    + ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_ADD &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll + node_lr) - #r -> node_lr - (#r - #ll), same type
                double v = ((struct parser_number*)(node->r))->value
                    - ((struct parser_number*)(node->l->l))->value;
                node->l = node->l->r;
                ((struct parser_number*)(node->r))->value = v;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll + #lr) - #r -> node_ll - (#r - #lr), same type
                double v = ((struct parser_number*)(node->r))->value
                    - ((struct parser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
            }
        }
        else if (node->l->type == PARSER_SUB &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll - node_lr) - #r -> (#ll - #r) - node_lr, type change
                double v = ((struct parser_number*)(node->l->l))->value
                    - ((struct parser_number*)(node->r))->value;
                node->r = node->l->r;
                node->l->type = PARSER_NUMBER;
                ((struct parser_number*)(node->l))->value = v;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll - #lr) - #r -> node_ll - (#r + #lr), same type
                double v = ((struct parser_number*)(node->r))->value
                    + ((struct parser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
            }
        }
        break;
    case PARSER_MUL:
    case PARSER_MUL_PP:
        parser_ast_optimize(node->l);
        parser_ast_optimize(node->r);
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_NUMBER)
        {
            double v = ((struct parser_number*)(node->l))->value
                *    ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.v = ((struct parser_number*)(node->l))->value;
            node->rip   = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_NUMBER)
        {
            node->lvp.v = ((struct parser_number*)(node->r))->value;
            node->rip   = ((struct parser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
            node->rip    = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_MUL_PP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL_VP)
        {
            double v = ((struct parser_number*)(node->l))->value * PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV_VP)
        {
            double v = ((struct parser_number*)(node->l))->value * PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_DIV_VP;
        }
        else if (node->l->type == PARSER_MUL_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) * ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_DIV_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) * ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_DIV_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l * (#rl * node_rr) -> (#l * #rl) * node_rr, same type
                double v = ((struct parser_number*)(node->l))->value
                    * ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l * (node_rl * #rr) -> (#l * #rr) * node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    * ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l * (#rl / node_rr) -> (#l * #rl) / node_rr, type change
                double v = ((struct parser_number*)(node->l))->value
                    * ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
                node->type = PARSER_DIV;

            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l * (node_rl / #rr) -> (#l / #rr) * node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    / ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_MUL &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll * node_lr) * #r -> nodel_lr * (#ll * #r), same type
                double v = ((struct parser_number*)(node->l->l))->value
                    * ((struct parser_number*)(node->r))->value;
                node->l = node->l->r;
                ((struct parser_number*)(node->r))->value = v;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll * #lr) * #r -> node_ll + (#lr * #r), same type
                double v = ((struct parser_number*)(node->l->r))->value
                    * ((struct parser_number*)(node->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
            }
        }
        else if (node->l->type == PARSER_DIV &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll / node_lr) * #r -> (#ll * #r) / node_lr, type change
                double v = ((struct parser_number*)(node->l->l))->value
                    * ((struct parser_number*)(node->r))->value;
                node->r = node->l->r;
                ((struct parser_number*)(node->l))->type = PARSER_NUMBER;
                ((struct parser_number*)(node->l))->value = v;
                node->type = PARSER_DIV;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll / #lr) * #r -> node_ll * (#r / #lr), same type
                double v = ((struct parser_number*)(node->r))->value
                    / ((struct parser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
            }
        }
        break;
    case PARSER_DIV:
    case PARSER_DIV_PP:
        parser_ast_optimize(node->l);
        parser_ast_optimize(node->r);
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_NUMBER)
        {
            double v = ((struct parser_number*)(node->l))->value
                /    ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.v = ((struct parser_number*)(node->l))->value;
            node->rip   = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_DIV_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_NUMBER)
        {
            node->lvp.v = double(1.)/((struct parser_number*)(node->r))->value;
            node->rip   =          ((struct parser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_SYMBOL &&
                 node->r->type == PARSER_SYMBOL)
        {
            node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
            node->rip    = ((struct parser_symbol*)(node->r))->ip;
            node->type = PARSER_DIV_PP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL_VP)
        {
            double v = ((struct parser_number*)(node->l))->value / PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_DIV_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV_VP)
        {
            double v = ((struct parser_number*)(node->l))->value / PARSER_EVAL_R(node);
            PARSER_MOVEUP_R(node, v);
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_MUL_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) / ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_DIV_VP &&
                 node->r->type == PARSER_NUMBER)
        {
            double v = PARSER_EVAL_L(node) / ((struct parser_number*)(node->r))->value;
            PARSER_MOVEUP_L(node, v);
            node->type = PARSER_DIV_VP;
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_MUL)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l / (#rl * node_rr) -> (#l / #rl) / node_rr, same type
                double v = ((struct parser_number*)(node->l))->value
                    / ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l / (node_rl * #rr) -> (#l / #rr) / node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    / ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_NUMBER &&
                 node->r->type == PARSER_DIV)
        {
            if (node->r->l->type == PARSER_NUMBER)
            { // #l / (#rl / node_rr) -> (#l / #rl) * node_rr, type change
                double v = ((struct parser_number*)(node->l))->value
                    / ((struct parser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct parser_number*)(node->l))->value = v;
                node->type = PARSER_MUL;
            }
            else if (node->r->r->type == PARSER_NUMBER)
            { // #l / (node_rl / #rr) -> (#l * #rr) / node_rl, same type
                double v = ((struct parser_number*)(node->l))->value
                    * ((struct parser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct parser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == PARSER_MUL &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll * node_lr) / #r -> node_lr * (#ll / #r), type change
                double v = ((struct parser_number*)(node->l->l))->value
                    / ((struct parser_number*)(node->r))->value;
                node->l = node->l->r;
                ((struct parser_number*)(node->r))->value = v;
                node->type = PARSER_MUL;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll * #lr) / #r -> node_ll * (#lr / #r), type change
                double v = ((struct parser_number*)(node->l->r))->value
                    / ((struct parser_number*)(node->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = v;
                node->type = PARSER_MUL;
            }
        }
        else if (node->l->type == PARSER_DIV &&
                 node->r->type == PARSER_NUMBER)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // (#ll / node_lr) / #r -> (#ll / #r) / node_lr, type change
                double v = ((struct parser_number*)(node->l->l))->value
                    / ((struct parser_number*)(node->r))->value;
                node->r = node->l->r;
                node->l->type = PARSER_NUMBER;
                ((struct parser_number*)(node->l))->value = v;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // (node_ll / #lr) / #r -> node_ll * 1./(#r * #lr), type change
                double v = ((struct parser_number*)(node->r))->value
                    * ((struct parser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct parser_number*)(node->r))->value = double(1.)/v;
                node->type = PARSER_MUL;
            }
        }
        break;
    case PARSER_NEG:
        parser_ast_optimize(node->l);
        if (node->l->type == PARSER_NUMBER)
        {
            double v = -((struct parser_number*)(node->l))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        else if (node->l->type == PARSER_SYMBOL)
        {
            node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
            node->type = PARSER_NEG_P;
        }
        else if (node->l->type == PARSER_ADD_VP)
        {
            PARSER_NEG_MOVEUP(node);
            node->type = PARSER_SUB_VP;
        }
        else if (node->l->type == PARSER_SUB_VP)
        {
            PARSER_NEG_MOVEUP(node);
            node->type = PARSER_ADD_VP;
        }
        else if (node->l->type == PARSER_MUL_VP)
        {
            PARSER_NEG_MOVEUP(node);
            node->type = PARSER_MUL_VP;
        }
        else if (node->l->type == PARSER_DIV_VP)
        {
            PARSER_NEG_MOVEUP(node);
            node->type = PARSER_DIV_VP;
        }
        else if (node->l->type == PARSER_ADD)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // -(#ll + node_lr) -> -#ll - node_lr
                node->r = node->l->r;
                ((struct parser_number*)(node->l))->value =
                    -((struct parser_number*)(node->l->l))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_SUB;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // -(node_ll + #lr) -> -#lr - node_ll
                node->r = node->l->l;
                ((struct parser_number*)(node->l))->value =
                    -((struct parser_number*)(node->l->r))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_SUB;
            }
        }
        else if (node->l->type == PARSER_SUB)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // -(#ll - node_lr) -> -#ll + node_lr
                node->r = node->l->r;
                ((struct parser_number*)(node->l))->value =
                    -((struct parser_number*)(node->l->l))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_ADD;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // -(node_ll - #lr) -> #lr - node_ll
                node->r = node->l->l;
                ((struct parser_number*)(node->l))->value =
                    ((struct parser_number*)(node->l->r))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_SUB;
            }
        }
        else if (node->l->type == PARSER_MUL)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // -(#ll * node_lr) -> -#ll * node_lr
                node->r = node->l->r;
                ((struct parser_number*)(node->l))->value =
                    -((struct parser_number*)(node->l->l))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_MUL;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // -(node_ll * #lr) -> -#lr * node_ll
                node->r = node->l->l;
                ((struct parser_number*)(node->l))->value =
                    -((struct parser_number*)(node->l->r))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_MUL;
            }
        }
        else if (node->l->type == PARSER_DIV)
        {
            if (node->l->l->type == PARSER_NUMBER)
            { // -(#ll / node_lr) -> -#ll / node_lr
                node->r = node->l->r;
                ((struct parser_number*)(node->l))->value =
                    -((struct parser_number*)(node->l->l))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_DIV;
            }
            else if (node->l->r->type == PARSER_NUMBER)
            { // -(node_ll / #lr) -> (-1/#lr) * node_ll
                node->r = node->l->l;
                ((struct parser_number*)(node->l))->value =
                    double(-1.0) / ((struct parser_number*)(node->l->r))->value;
                node->l->type = PARSER_NUMBER;
                node->type = PARSER_MUL;
            }
        }
        break;
    case PARSER_F1:
        parser_ast_optimize(node->l);
        if (node->l->type == PARSER_NUMBER)
        {
            double v = parser_call_f1
                (((struct parser_f1*)node)->ftype,
                 ((struct parser_number*)(((struct parser_f1*)node)->l))->value);
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_F2:
        parser_ast_optimize(node->l);
        parser_ast_optimize(node->r);
        if (node->l->type == PARSER_NUMBER &&
            node->r->type == PARSER_NUMBER)
        {
            double v = parser_call_f2
                (((struct parser_f2*)node)->ftype,
                 ((struct parser_number*)(((struct parser_f2*)node)->l))->value,
                 ((struct parser_number*)(((struct parser_f2*)node)->r))->value);
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        else if (node->r->type == PARSER_NUMBER && ((struct parser_f2*)node)->ftype == PARSER_POW)
        {
            struct parser_node* n = node->l;
            double v = ((struct parser_number*)(node->r))->value;
            if (-3.0 == v) {
                ((struct parser_f1*)node)->type = PARSER_F1;
                ((struct parser_f1*)node)->l = n;
                ((struct parser_f1*)node)->ftype = PARSER_POW_M3;
            } else if (-2.0 == v) {
                ((struct parser_f1*)node)->type = PARSER_F1;
                ((struct parser_f1*)node)->l = n;
                ((struct parser_f1*)node)->ftype = PARSER_POW_M2;
            } else if (-1.0 == v) {
                ((struct parser_f1*)node)->type = PARSER_F1;
                ((struct parser_f1*)node)->l = n;
                ((struct parser_f1*)node)->ftype = PARSER_POW_M1;
            } else if (0.0 == v) {
                ((struct parser_number*)node)->type = PARSER_NUMBER;
                ((struct parser_number*)node)->value = 1.0;
            } else if (1.0 == v) {
                ((struct parser_f1*)node)->type = PARSER_F1;
                ((struct parser_f1*)node)->l = n;
                ((struct parser_f1*)node)->ftype = PARSER_POW_P1;
            } else if (2.0 == v) {
                ((struct parser_f1*)node)->type = PARSER_F1;
                ((struct parser_f1*)node)->l = n;
                ((struct parser_f1*)node)->ftype = PARSER_POW_P2;
            } else if (3.0 == v) {
                ((struct parser_f1*)node)->type = PARSER_F1;
                ((struct parser_f1*)node)->l = n;
                ((struct parser_f1*)node)->ftype = PARSER_POW_P3;
            }
        }
        break;
    case PARSER_F3:
        parser_ast_optimize(((struct parser_f3*)node)->n1);
        parser_ast_optimize(((struct parser_f3*)node)->n2);
        parser_ast_optimize(((struct parser_f3*)node)->n3);
        if (((struct parser_f3*)node)->n1->type == PARSER_NUMBER &&
            ((struct parser_f3*)node)->n2->type == PARSER_NUMBER &&
            ((struct parser_f3*)node)->n3->type == PARSER_NUMBER)
        {
            double v = parser_call_f3
                (((struct parser_f3*)node)->ftype,
                 ((struct parser_number*)(((struct parser_f3*)node)->n1))->value,
                 ((struct parser_number*)(((struct parser_f3*)node)->n2))->value,
                 ((struct parser_number*)(((struct parser_f3*)node)->n3))->value);
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_ADD_VP:
        parser_ast_optimize(node->r);
        if (node->r->type == PARSER_NUMBER)
        {
            double v = node->lvp.v + ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_SUB_VP:
        parser_ast_optimize(node->r);
        if (node->r->type == PARSER_NUMBER)
        {
            double v = node->lvp.v - ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_MUL_VP:
        parser_ast_optimize(node->r);
        if (node->r->type == PARSER_NUMBER)
        {
            double v = node->lvp.v * ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_DIV_VP:
        parser_ast_optimize(node->r);
        if (node->r->type == PARSER_NUMBER)
        {
            double v = node->lvp.v / ((struct parser_number*)(node->r))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_NEG_P:
        parser_ast_optimize(node->l);
        if (node->l->type == PARSER_NUMBER)
        {
            double v = -((struct parser_number*)(node->l))->value;
            ((struct parser_number*)node)->type = PARSER_NUMBER;
            ((struct parser_number*)node)->value = v;
        }
        break;
    case PARSER_ASSIGN:
        parser_ast_optimize(((struct parser_assign*)node)->v);
        break;
    case PARSER_LIST:
        parser_ast_optimize(node->l);
        parser_ast_optimize(node->r);
        break;
    default:
        amrex::Abort("parser_ast_optimize: unknown node type " + std::to_string(node->type));
    }
}

static
void
parser_ast_print_f1 (struct parser_f1* f1, std::string const& space, AllPrint& printer)
{
    printer << space;
    switch (f1->ftype) {
    case PARSER_SQRT:        printer << "SQRT\n";        break;
    case PARSER_EXP:         printer << "EXP\n";         break;
    case PARSER_LOG:         printer << "LOG\n";         break;
    case PARSER_LOG10:       printer << "LOG10\n";       break;
    case PARSER_SIN:         printer << "SIN\n";         break;
    case PARSER_COS:         printer << "COS\n";         break;
    case PARSER_TAN:         printer << "TAN\n";         break;
    case PARSER_ASIN:        printer << "ASIN\n";        break;
    case PARSER_ACOS:        printer << "ACOS\n";        break;
    case PARSER_ATAN:        printer << "ATAN\n";        break;
    case PARSER_SINH:        printer << "SINH\n";        break;
    case PARSER_COSH:        printer << "COSH\n";        break;
    case PARSER_TANH:        printer << "TANH\n";        break;
    case PARSER_ABS:         printer << "ABS\n";         break;
    case PARSER_FLOOR:       printer << "FLOOR\n";       break;
    case PARSER_CEIL:        printer << "CEIL\n";        break;
    case PARSER_POW_M3:      printer << "POW(,-3)\n";    break;
    case PARSER_POW_M2:      printer << "POW(,-2)\n";    break;
    case PARSER_POW_M1:      printer << "POW(,-1)\n";    break;
    case PARSER_POW_P1:      printer << "POW(,1)\n";     break;
    case PARSER_POW_P2:      printer << "POW(,2)\n";     break;
    case PARSER_POW_P3:      printer << "POW(,3)\n";     break;
    case PARSER_COMP_ELLINT_1: printer << "COMP_ELLINT_1\n"; break;
    case PARSER_COMP_ELLINT_2: printer << "COMP_ELLINT_2\n"; break;
    default:
        amrex::AllPrint() << "parser_ast_print_f1: Unknown function " << f1->ftype << "\n";
    }
    parser_ast_print(f1->l, space+"  ", printer);
}

static
void
parser_ast_print_f2 (struct parser_f2* f2, std::string const& space, AllPrint& printer)
{
    printer << space;
    switch (f2->ftype) {
    case PARSER_POW:
        printer << "POW\n";
        break;
    case PARSER_GT:
        printer << "GT\n";
        break;
    case PARSER_LT:
        printer << "LT\n";
        break;
    case PARSER_GEQ:
        printer << "GEQ\n";
        break;
    case PARSER_LEQ:
        printer << "LEQ\n";
        break;
    case PARSER_EQ:
        printer << "EQ\n";
        break;
    case PARSER_NEQ:
        printer << "NEQ\n";
        break;
    case PARSER_AND:
        printer << "AND\n";
        break;
    case PARSER_OR:
        printer << "OR\n";
        break;
    case PARSER_HEAVISIDE:
        printer << "HEAVISIDE\n";
        break;
    case PARSER_JN:
        printer << "JN\n";
        break;
    case PARSER_MIN:
        printer << "MIN\n";
        break;
    case PARSER_MAX:
        printer << "MAX\n";
        break;
    case PARSER_FMOD:
        printer << "FMOD\n";
        break;
    default:
        amrex::AllPrint() << "parser_ast_print_f2: Unknown function " << f2->ftype << "\n";
    }
    parser_ast_print(f2->l, space+"  ", printer);
    parser_ast_print(f2->r, space+"  ", printer);
}

static
void
parser_ast_print_f3 (struct parser_f3* f3, std::string const& space, AllPrint& printer)
{
    std::string const& more_space = space + "  ";
    switch (f3->ftype) {
    case PARSER_IF:
        printer << space << "IF\n";
        break;
    default:
        amrex::AllPrint() << "parser_ast_print_f3: Unknown function " << f3->ftype << "\n";
    }
    parser_ast_print(f3->n1, more_space, printer);
    parser_ast_print(f3->n2, more_space, printer);
    parser_ast_print(f3->n3, more_space, printer);
}

void
parser_ast_print (struct parser_node* node, std::string const& space, AllPrint& printer)
{
    std::string const& more_space = space + "  ";
    switch (node->type)
    {
    case PARSER_NUMBER:
        printer << space << "NUMBER: " << ((struct parser_number*)node)->value << "\n";
        break;
    case PARSER_SYMBOL:
        printer << space << "VARIABLE: " << ((struct parser_symbol*)node)->name << "\n";
        break;
    case PARSER_ADD:
        printer << space << "ADD\n";
        parser_ast_print(node->l, more_space, printer);
        parser_ast_print(node->r, more_space, printer);
        break;
    case PARSER_SUB:
        printer << space << "SUB\n";
        parser_ast_print(node->l, more_space, printer);
        parser_ast_print(node->r, more_space, printer);
        break;
    case PARSER_MUL:
        printer << space << "MUL\n";
        parser_ast_print(node->l, more_space, printer);
        parser_ast_print(node->r, more_space, printer);
        break;
    case PARSER_DIV:
        printer << space << "DIV\n";
        parser_ast_print(node->l, more_space, printer);
        parser_ast_print(node->r, more_space, printer);
        break;
    case PARSER_NEG:
        printer << space << "NEG\n";
        parser_ast_print(node->l, more_space, printer);
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
    case PARSER_ASSIGN:
        printer << space << "=: " << ((struct parser_assign*)node)->s->name << " =\n";
        parser_ast_print(((struct parser_assign*)node)->v, more_space, printer);
        break;
    case PARSER_LIST:
        printer << space <<"LIST\n";
        parser_ast_print(node->l, more_space, printer);
        parser_ast_print(node->r, more_space, printer);
        break;
    case PARSER_ADD_VP:
        printer << space << "ADD: " << node->lvp.v << " "
                << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_SUB_VP:
        printer << space << "SUB: " << node->lvp.v << " "
                << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_MUL_VP:
        printer << space << "MUL: " << node->lvp.v << " "
                << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_DIV_VP:
        printer << space << "DIV: " << node->lvp.v << " "
                << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_NEG_P:
        printer << space << "NEG: " << ((struct parser_symbol*)(node->l))->name << "\n";
        break;
    case PARSER_ADD_PP:
        printer << space << "ADD: " << ((struct parser_symbol*)(node->l))->name
                << "  "             << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_SUB_PP:
        printer << space << "SUB: " << ((struct parser_symbol*)(node->l))->name
                << "  "             << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_MUL_PP:
        printer << space << "MUL: " << ((struct parser_symbol*)(node->l))->name
                << "  "             << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    case PARSER_DIV_PP:
        printer << space << "DIV: " << ((struct parser_symbol*)(node->l))->name
                << "  "             << ((struct parser_symbol*)(node->r))->name << "\n";
        break;
    default:
        amrex::Abort("parser_ast_print: unknown node type " + std::to_string(node->type));
    }
}

int
parser_ast_depth (struct parser_node* node)
{
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
    case PARSER_NEG:
        return parser_ast_depth(node->l)+1;
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
    case PARSER_ASSIGN:
    {
        int d = parser_ast_depth(((struct parser_assign*)node)->v);
        return d+1;
    }
    case PARSER_ADD_VP:
    case PARSER_SUB_VP:
    case PARSER_MUL_VP:
    case PARSER_DIV_VP:
    case PARSER_NEG_P:
    case PARSER_ADD_PP:
    case PARSER_SUB_PP:
    case PARSER_MUL_PP:
    case PARSER_DIV_PP:
        return 1;
    default:
        amrex::Abort("parser_ast_print: unknown node type " + std::to_string(node->type));
        return 0;
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
    case PARSER_NEG:
        parser_ast_regvar(node->l, name, i);
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
    case PARSER_ASSIGN:
        parser_ast_regvar(((struct parser_assign*)node)->v, name, i);
        break;
    case PARSER_ADD_VP:
    case PARSER_SUB_VP:
    case PARSER_MUL_VP:
    case PARSER_DIV_VP:
        parser_ast_regvar(node->r, name, i);
        node->rip = ((struct parser_symbol*)(node->r))->ip;
        break;
    case PARSER_NEG_P:
        parser_ast_regvar(node->l, name, i);
        node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
        break;
    case PARSER_ADD_PP:
    case PARSER_SUB_PP:
    case PARSER_MUL_PP:
    case PARSER_DIV_PP:
        parser_ast_regvar(node->l, name, i);
        parser_ast_regvar(node->r, name, i);
        node->lvp.ip = ((struct parser_symbol*)(node->l))->ip;
        node->rip = ((struct parser_symbol*)(node->r))->ip;
        break;
    default:
        amrex::AllPrint() << "parser_ast_regvar: unknown node type " << node->type << "\n";
        amrex::Abort();
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
    case PARSER_ADD_PP:
    case PARSER_SUB_PP:
    case PARSER_MUL_PP:
    case PARSER_DIV_PP:
    case PARSER_LIST:
        parser_ast_setconst(node->l, name, c);
        parser_ast_setconst(node->r, name, c);
        break;
    case PARSER_NEG:
    case PARSER_NEG_P:
        parser_ast_setconst(node->l, name, c);
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
    case PARSER_ASSIGN:
        parser_ast_setconst(((struct parser_assign*)node)->v, name, c);
        break;
    case PARSER_ADD_VP:
    case PARSER_SUB_VP:
    case PARSER_MUL_VP:
    case PARSER_DIV_VP:
        parser_ast_setconst(node->r, name, c);
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
    case PARSER_ADD_PP:
    case PARSER_SUB_PP:
    case PARSER_MUL_PP:
    case PARSER_DIV_PP:
    case PARSER_LIST:
        parser_ast_get_symbols(node->l, symbols, local_symbols);
        parser_ast_get_symbols(node->r, symbols, local_symbols);
        break;
    case PARSER_NEG:
    case PARSER_NEG_P:
        parser_ast_get_symbols(node->l, symbols, local_symbols);
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
    case PARSER_ASSIGN:
        local_symbols.emplace(((struct parser_assign*)node)->s->name);
        parser_ast_get_symbols(((struct parser_assign*)node)->v, symbols, local_symbols);
        break;
    case PARSER_ADD_VP:
    case PARSER_SUB_VP:
    case PARSER_MUL_VP:
    case PARSER_DIV_VP:
        parser_ast_get_symbols(node->r, symbols, local_symbols);
        break;
    default:
        amrex::Abort("parser_ast_get_symbols: unknown node type " + std::to_string(node->type));
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
    parser_ast_optimize(parser->ast);
}

void
parser_print (struct amrex_parser* parser)
{
    amrex::AllPrint printer{};
    printer.SetPrecision(17);
    parser_ast_print(parser->ast, std::string("  "), printer);
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

int
parser_depth (struct amrex_parser* parser)
{
    return parser_ast_depth(parser->ast);
}

}
