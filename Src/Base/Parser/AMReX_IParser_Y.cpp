#include <AMReX.H>
#include <AMReX_IParser_Y.H>
#include <amrex_iparser.tab.h>

#include <algorithm>
#include <cstdarg>
#include <string>

void
amrex_iparsererror (char const *s, ...)
{
    char print_buff[512];
    std::va_list vl;
    va_start(vl, s);
    std::vsnprintf(print_buff, 512, s, vl);
    va_end(vl);
    throw std::runtime_error(print_buff);
}

namespace amrex {

static struct iparser_node* iparser_root = nullptr;

// This is called by a bison rule to store the original AST in a static variable.
void
iparser_defexpr (struct iparser_node* body)
{
    iparser_root = body;
}

struct iparser_symbol*
iparser_makesymbol (char* name)
{
    auto *symbol = (struct iparser_symbol*) std::malloc(sizeof(struct iparser_symbol));
    symbol->type = IPARSER_SYMBOL;
    symbol->name = strdup(name);
    symbol->ip = -1;
    return symbol;
}

struct iparser_node*
iparser_newnode (enum iparser_node_t type, struct iparser_node* l, struct iparser_node* r)
{
    auto *tmp = (struct iparser_node*) std::malloc(sizeof(struct iparser_node));
    tmp->type = type;
    tmp->l = l;
    tmp->r = r;
    return tmp;
}

struct iparser_node*
iparser_newnumber (int d)
{
    auto *r = (struct iparser_number*) std::malloc(sizeof(struct iparser_number));
    r->type = IPARSER_NUMBER;
    r->value = d;
    return (struct iparser_node*) r;
}

struct iparser_node*
iparser_newsymbol (struct iparser_symbol* symbol)
{
    return (struct iparser_node*) symbol;
}

struct iparser_node*
iparser_newf1 (enum iparser_f1_t ftype, struct iparser_node* l)
{
    auto *tmp = (struct iparser_f1*) std::malloc(sizeof(struct iparser_f1));
    tmp->type = IPARSER_F1;
    tmp->l = l;
    tmp->ftype = ftype;
    return (struct iparser_node*) tmp;
}

struct iparser_node*
iparser_newf2 (enum iparser_f2_t ftype, struct iparser_node* l, struct iparser_node* r)
{
    auto *tmp = (struct iparser_f2*) std::malloc(sizeof(struct iparser_f2));
    tmp->type = IPARSER_F2;
    tmp->l = l;
    tmp->r = r;
    tmp->ftype = ftype;
    return (struct iparser_node*) tmp;
}

struct iparser_node*
iparser_newf3 (enum iparser_f3_t ftype, struct iparser_node* n1, struct iparser_node* n2,
               struct iparser_node* n3)
{
    auto *tmp = (struct iparser_f3*) std::malloc(sizeof(struct iparser_f3));
    tmp->type = IPARSER_F3;
    tmp->n1 = n1;
    tmp->n2 = n2;
    tmp->n3 = n3;
    tmp->ftype = ftype;
    return (struct iparser_node*) tmp;
}

struct iparser_node*
iparser_newassign (struct iparser_symbol* sym, struct iparser_node* v)
{
    auto *r = (struct iparser_assign*) std::malloc(sizeof(struct iparser_assign));
    r->type = IPARSER_ASSIGN;
    r->s = sym;
    r->v = v;
    return (struct iparser_node*) r;
}

struct iparser_node*
iparser_newlist (struct iparser_node* nl, struct iparser_node* nr)
{
    if (nr == nullptr) {
        return nl;
    } else {
        auto *r = (struct iparser_node*) std::malloc(sizeof(struct iparser_node));
        r->type = IPARSER_LIST;
        r->l = nl;
        r->r = nr;
        return r;
    }
}

/*******************************************************************/

struct amrex_iparser*
amrex_iparser_new ()
{
    auto *my_iparser = (struct amrex_iparser*) std::malloc(sizeof(struct amrex_iparser));

    my_iparser->sz_mempool = iparser_ast_size(iparser_root);
    my_iparser->p_root = std::malloc(my_iparser->sz_mempool);
    my_iparser->p_free = my_iparser->p_root;

    my_iparser->ast = iparser_ast_dup(my_iparser, iparser_root, 1); /* 1: free the source iparser_root */

    if ((char*)my_iparser->p_root + my_iparser->sz_mempool != (char*)my_iparser->p_free) {
        amrex::Abort("amrex_iparser_new: error in memory size");
    }

    iparser_ast_optimize(my_iparser->ast);

    return my_iparser;
}

void
amrex_iparser_delete (struct amrex_iparser* iparser)
{
    std::free(iparser->p_root);
    std::free(iparser);
}

static
std::size_t
iparser_aligned_size (std::size_t N)
{
    const unsigned int align_size = 16;
    std::size_t x = N + (align_size-1);
    x -= x & (align_size-1);
    return x;
}

static
void*
iparser_allocate (struct amrex_iparser* my_iparser, std::size_t N)
{
    void* r = my_iparser->p_free;
    my_iparser->p_free = (char*)r + iparser_aligned_size(N);
    return r;
}

struct amrex_iparser*
iparser_dup (struct amrex_iparser* source)
{
    auto *dest = (struct amrex_iparser*) std::malloc(sizeof(struct amrex_iparser));
    dest->sz_mempool = source->sz_mempool;
    dest->p_root = std::malloc(dest->sz_mempool);
    dest->p_free = dest->p_root;

    dest->ast = iparser_ast_dup(dest, source->ast, 0); /* 0: don't free the source */

    return dest;
}

std::size_t
iparser_ast_size (struct iparser_node* node)
{
    std::size_t result = 0;

    switch (node->type)
    {
    case IPARSER_NUMBER:
        result = iparser_aligned_size(sizeof(struct iparser_number));
        break;
    case IPARSER_SYMBOL:
        result = iparser_aligned_size(    sizeof(struct iparser_symbol))
            + iparser_aligned_size(std::strlen(((struct iparser_symbol*)node)->name)+1);
        break;
    case IPARSER_ADD:
    case IPARSER_SUB:
    case IPARSER_MUL:
    case IPARSER_DIV:
    case IPARSER_ADD_PP:
    case IPARSER_SUB_PP:
    case IPARSER_MUL_PP:
    case IPARSER_DIV_PP:
    case IPARSER_LIST:
        result = iparser_aligned_size(sizeof(struct iparser_node))
            + iparser_ast_size(node->l) + iparser_ast_size(node->r);
        break;
    case IPARSER_NEG:
        result = iparser_aligned_size(sizeof(struct iparser_node))
            + iparser_ast_size(node->l);
        break;
    case IPARSER_F1:
        result = iparser_aligned_size(sizeof(struct iparser_f1))
            +             iparser_ast_size(((struct iparser_f1*)node)->l);
        break;
    case IPARSER_F2:
        result = iparser_aligned_size(sizeof(struct iparser_f2))
            +             iparser_ast_size(((struct iparser_f2*)node)->l)
            +             iparser_ast_size(((struct iparser_f2*)node)->r);
        break;
    case IPARSER_F3:
        result = iparser_aligned_size(sizeof(struct iparser_f3))
            +             iparser_ast_size(((struct iparser_f3*)node)->n1)
            +             iparser_ast_size(((struct iparser_f3*)node)->n2)
            +             iparser_ast_size(((struct iparser_f3*)node)->n3);
        break;
    case IPARSER_ASSIGN:
        result += iparser_aligned_size(sizeof(struct iparser_assign))
            + iparser_ast_size((struct iparser_node*)(((struct iparser_assign*)node)->s))
            + iparser_ast_size(((struct iparser_assign*)node)->v);
        break;
    case IPARSER_ADD_VP:
    case IPARSER_SUB_VP:
    case IPARSER_MUL_VP:
    case IPARSER_DIV_VP:
    case IPARSER_DIV_PV:
        result = iparser_aligned_size(sizeof(struct iparser_node))
            + iparser_ast_size(node->r);
        break;
    case IPARSER_NEG_P:
        result = iparser_aligned_size(sizeof(struct iparser_node))
            + iparser_ast_size(node->l);
        break;
    default:
        amrex::Abort("iparser_ast_size: unknown node type " + std::to_string(node->type));
    }

    return result;
}

struct iparser_node*
iparser_ast_dup (struct amrex_iparser* my_iparser, struct iparser_node* node, int move)
{
    void* result = nullptr;

    switch (node->type)
    {
    case IPARSER_NUMBER:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_number));
        std::memcpy(result, node            , sizeof(struct iparser_number));
        break;
    case IPARSER_SYMBOL:
    {
        result = iparser_allocate(my_iparser, sizeof(struct iparser_symbol));
        std::memcpy(result, node            , sizeof(struct iparser_symbol));
        const auto len = std::strlen(((struct iparser_symbol*)node)->name)+1;
        ((struct iparser_symbol*)result)->name = (char*) iparser_allocate
            (my_iparser, len);
        std::strncpy(((struct iparser_symbol*)result)->name,
                     ((struct iparser_symbol*)node  )->name, len);
        break;
    }
    case IPARSER_ADD:
    case IPARSER_SUB:
    case IPARSER_MUL:
    case IPARSER_DIV:
    case IPARSER_ADD_PP:
    case IPARSER_SUB_PP:
    case IPARSER_MUL_PP:
    case IPARSER_DIV_PP:
    case IPARSER_LIST:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_node));
        std::memcpy(result, node            , sizeof(struct iparser_node));
        ((struct iparser_node*)result)->l = iparser_ast_dup(my_iparser, node->l, move);
        ((struct iparser_node*)result)->r = iparser_ast_dup(my_iparser, node->r, move);
        break;
    case IPARSER_NEG:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_node));
        std::memcpy(result, node            , sizeof(struct iparser_node));
        ((struct iparser_node*)result)->l = iparser_ast_dup(my_iparser, node->l, move);
        ((struct iparser_node*)result)->r = nullptr;
        break;
    case IPARSER_F1:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_f1));
        std::memcpy(result, node            , sizeof(struct iparser_f1));
        ((struct iparser_f1*)result)->l = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_f1*)node)->l, move);
        break;
    case IPARSER_F2:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_f2));
        std::memcpy(result, node            , sizeof(struct iparser_f2));
        ((struct iparser_f2*)result)->l = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_f2*)node)->l, move);
        ((struct iparser_f2*)result)->r = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_f2*)node)->r, move);
        break;
    case IPARSER_F3:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_f3));
        std::memcpy(result, node            , sizeof(struct iparser_f3));
        ((struct iparser_f3*)result)->n1 = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_f3*)node)->n1, move);
        ((struct iparser_f3*)result)->n2 = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_f3*)node)->n2, move);
        ((struct iparser_f3*)result)->n3 = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_f3*)node)->n3, move);
        break;
    case IPARSER_ASSIGN:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_assign));
        std::memcpy(result, node            , sizeof(struct iparser_assign));
        ((struct iparser_assign*)result)->s = (struct iparser_symbol*)
            iparser_ast_dup(my_iparser, (struct iparser_node*)
                                                  (((struct iparser_assign*)node)->s), move);
        ((struct iparser_assign*)result)->v = iparser_ast_dup(my_iparser,
                                                   ((struct iparser_assign*)node)->v, move);
        break;
    case IPARSER_ADD_VP:
    case IPARSER_SUB_VP:
    case IPARSER_MUL_VP:
    case IPARSER_DIV_VP:
    case IPARSER_DIV_PV:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_node));
        std::memcpy(result, node            , sizeof(struct iparser_node));
        ((struct iparser_node*)result)->r = iparser_ast_dup(my_iparser, node->r, move);
        break;
    case IPARSER_NEG_P:
        result = iparser_allocate(my_iparser, sizeof(struct iparser_node));
        std::memcpy(result, node            , sizeof(struct iparser_node));
        ((struct iparser_node*)result)->l = iparser_ast_dup(my_iparser, node->l, move);
        break;
    default:
        amrex::Abort("iparser_ast_dup: unknown node type " + std::to_string(node->type));
    }
    if (move) {
        /* Note that we only do this on the original AST.  We do not
         * need to call free for AST stored in amrex_iparser because the
         * memory is not allocated with std::malloc directly.
         */
        if (node->type == IPARSER_SYMBOL) {
            std::free(((struct iparser_symbol*)node)->name);
        }
        std::free((void*)node);
    }
    return (struct iparser_node*)result;
}

#define IPARSER_MOVEUP_R(node, v) \
    struct iparser_node* n = (node)->r->r; \
    int ip = (node)->r->rip; \
    (node)->r = n; \
    (node)->lvp.v = v; \
    (node)->rip   = ip;
#define IPARSER_MOVEUP_L(node, v) \
    struct iparser_node* n = (node)->l->r; \
    int ip = (node)->l->rip; \
    (node)->r = n; \
    (node)->lvp.v = v; \
    (node)->rip   = ip;
#define IPARSER_EVAL_R(node) (node)->r->lvp.v
#define IPARSER_EVAL_L(node) (node)->l->lvp.v

#define IPARSER_NEG_MOVEUP(node) \
    (node)->r = (node)->l->r; \
    (node)->lvp.v = -(node)->l->lvp.v; \
    (node)->rip = (node)->l->rip;

void
iparser_ast_optimize (struct iparser_node* node)
{
    /* No need to free memory because we only call this on ASTs in
     * amrex_iparser that are allocated from the memory pool.
     */
    switch (node->type)
    {
    case IPARSER_NUMBER:
    case IPARSER_SYMBOL:
        break;
    case IPARSER_ADD:
    case IPARSER_ADD_PP:
        iparser_ast_optimize(node->l);
        iparser_ast_optimize(node->r);
        if (node->l->type == IPARSER_NUMBER &&
            node->r->type == IPARSER_NUMBER)
        {
            int v = ((struct iparser_number*)(node->l))->value
                +   ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.v = ((struct iparser_number*)(node->l))->value;
            node->rip   = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_NUMBER)
        {
            node->lvp.v = ((struct iparser_number*)(node->r))->value;
            node->rip   = ((struct iparser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
            node->rip    = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_ADD_PP; // For *_PP, the names are stored in the l and r nodes.
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_ADD_VP)
        {
            int v = ((struct iparser_number*)(node->l))->value + IPARSER_EVAL_R(node);
            IPARSER_MOVEUP_R(node, v);
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SUB_VP)
        {
            int v = ((struct iparser_number*)(node->l))->value + IPARSER_EVAL_R(node);
            IPARSER_MOVEUP_R(node, v);
            node->type = IPARSER_SUB_VP;
        }
        else if (node->l->type == IPARSER_ADD_VP &&
                 node->r->type == IPARSER_NUMBER)
        {
            int v = IPARSER_EVAL_L(node) + ((struct iparser_number*)(node->r))->value;
            IPARSER_MOVEUP_L(node, v);
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_SUB_VP &&
                 node->r->type == IPARSER_NUMBER)
        {
            int v = IPARSER_EVAL_L(node) + ((struct iparser_number*)(node->r))->value;
            IPARSER_MOVEUP_L(node, v);
            node->type = IPARSER_SUB_VP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_ADD)
        {
            if (node->r->l->type == IPARSER_NUMBER)
            { // #l + (#rl + node_rr) -> (#l + #rl) + node_rr, same type
                int v = ((struct iparser_number*)(node->l))->value
                    +   ((struct iparser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct iparser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == IPARSER_NUMBER)
            { // #l + (node_rl + #rr) -> (#l + #rr) + node_rl, same type
                int v = ((struct iparser_number*)(node->l))->value
                    +   ((struct iparser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct iparser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SUB)
        {
            if (node->r->l->type == IPARSER_NUMBER)
            { // #l + (#rl - node_rr) -> (#l + #rl) - node_rr, type change
                int v = ((struct iparser_number*)(node->l))->value
                    +   ((struct iparser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct iparser_number*)(node->l))->value = v;
                node->type = IPARSER_SUB;
            }
            else if (node->r->r->type == IPARSER_NUMBER)
            { // #l + (node_rl - #rr) -> (#l - #rr) + node_rl, same type
                int v = ((struct iparser_number*)(node->l))->value
                    -   ((struct iparser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct iparser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == IPARSER_ADD &&
                 node->r->type == IPARSER_NUMBER)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // (#ll + node_lr) + #r -> nodel_lr + (#ll + #r), same type
                int v = ((struct iparser_number*)(node->l->l))->value
                    +   ((struct iparser_number*)(node->r))->value;
                node->l = node->l->r;
                ((struct iparser_number*)(node->r))->value = v;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // (node_ll + #lr) + #r -> node_ll + (#lr + #r), same type
                int v = ((struct iparser_number*)(node->l->r))->value
                    +   ((struct iparser_number*)(node->r))->value;
                node->l = node->l->l;
                ((struct iparser_number*)(node->r))->value = v;
            }
        }
        else if (node->l->type == IPARSER_SUB &&
                 node->r->type == IPARSER_NUMBER)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // (#ll - node_lr) + #r -> (#ll + #r) - node_lr, type change
                int v = ((struct iparser_number*)(node->l->l))->value
                    +   ((struct iparser_number*)(node->r))->value;
                node->r = node->l->r;
                ((struct iparser_number*)(node->l))->type = IPARSER_NUMBER;
                ((struct iparser_number*)(node->l))->value = v;
                node->type = IPARSER_SUB;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // (node_ll - #lr) + #r -> node_ll + (#r - #lr), same type
                int v = ((struct iparser_number*)(node->r))->value
                    -   ((struct iparser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct iparser_number*)(node->r))->value = v;
            }
        }
        break;
    case IPARSER_SUB:
    case IPARSER_SUB_PP:
        iparser_ast_optimize(node->l);
        iparser_ast_optimize(node->r);
        if (node->l->type == IPARSER_NUMBER &&
            node->r->type == IPARSER_NUMBER)
        {
            int v = ((struct iparser_number*)(node->l))->value
                -   ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.v = ((struct iparser_number*)(node->l))->value;
            node->rip   = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_SUB_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_NUMBER)
        {
            node->lvp.v = -((struct iparser_number*)(node->r))->value;
            node->rip   =  ((struct iparser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
            node->rip    = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_SUB_PP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_ADD_VP)
        {
            int v = ((struct iparser_number*)(node->l))->value - IPARSER_EVAL_R(node);
            IPARSER_MOVEUP_R(node, v);
            node->type = IPARSER_SUB_VP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SUB_VP)
        {
            int v = ((struct iparser_number*)(node->l))->value - IPARSER_EVAL_R(node);
            IPARSER_MOVEUP_R(node, v);
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_ADD_VP &&
                 node->r->type == IPARSER_NUMBER)
        {
            int v = IPARSER_EVAL_L(node) - ((struct iparser_number*)(node->r))->value;
            IPARSER_MOVEUP_L(node, v);
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_SUB_VP &&
                 node->r->type == IPARSER_NUMBER)
        {
            int v = IPARSER_EVAL_L(node) - ((struct iparser_number*)(node->r))->value;
            IPARSER_MOVEUP_L(node, v);
            node->type = IPARSER_SUB_VP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_ADD)
        {
            if (node->r->l->type == IPARSER_NUMBER)
            { // #l - (#rl + node_rr) -> (#l - #rl) - node_rr, same type
                int v = ((struct iparser_number*)(node->l))->value
                    -   ((struct iparser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct iparser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == IPARSER_NUMBER)
            { // #l - (node_rl + #rr) -> (#l - #rr) - node_rl, same type
                int v = ((struct iparser_number*)(node->l))->value
                    -   ((struct iparser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct iparser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SUB)
        {
            if (node->r->l->type == IPARSER_NUMBER)
            { // #l - (#rl - node_rr) -> (#l - #rl) + node_rr, type change
                int v = ((struct iparser_number*)(node->l))->value
                    -   ((struct iparser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct iparser_number*)(node->l))->value = v;
                node->type = IPARSER_ADD;
            }
            else if (node->r->r->type == IPARSER_NUMBER)
            { // #l - (node_rl - #rr) -> (#l + #rr) - node_rl, same type
                int v = ((struct iparser_number*)(node->l))->value
                    +   ((struct iparser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct iparser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == IPARSER_ADD &&
                 node->r->type == IPARSER_NUMBER)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // (#ll + node_lr) - #r -> node_lr - (#r - #ll), same type
                int v = ((struct iparser_number*)(node->r))->value
                    -   ((struct iparser_number*)(node->l->l))->value;
                node->l = node->l->r;
                ((struct iparser_number*)(node->r))->value = v;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // (node_ll + #lr) - #r -> node_ll - (#r - #lr), same type
                int v = ((struct iparser_number*)(node->r))->value
                    -   ((struct iparser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct iparser_number*)(node->r))->value = v;
            }
        }
        else if (node->l->type == IPARSER_SUB &&
                 node->r->type == IPARSER_NUMBER)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // (#ll - node_lr) - #r -> (#ll - #r) - node_lr, type change
                int v = ((struct iparser_number*)(node->l->l))->value
                    -   ((struct iparser_number*)(node->r))->value;
                node->r = node->l->r;
                node->l->type = IPARSER_NUMBER;
                ((struct iparser_number*)(node->l))->value = v;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // (node_ll - #lr) - #r -> node_ll - (#r + #lr), same type
                int v = ((struct iparser_number*)(node->r))->value
                    +   ((struct iparser_number*)(node->l->r))->value;
                node->l = node->l->l;
                ((struct iparser_number*)(node->r))->value = v;
            }
        }
        break;
    case IPARSER_MUL:
    case IPARSER_MUL_PP:
        iparser_ast_optimize(node->l);
        iparser_ast_optimize(node->r);
        if (node->l->type == IPARSER_NUMBER &&
            node->r->type == IPARSER_NUMBER)
        {
            int v = ((struct iparser_number*)(node->l))->value
                *   ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.v = ((struct iparser_number*)(node->l))->value;
            node->rip   = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_MUL_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_NUMBER)
        {
            node->lvp.v = ((struct iparser_number*)(node->r))->value;
            node->rip   = ((struct iparser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = IPARSER_MUL_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
            node->rip    = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_MUL_PP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_MUL_VP)
        {
            int v = ((struct iparser_number*)(node->l))->value * IPARSER_EVAL_R(node);
            IPARSER_MOVEUP_R(node, v);
            node->type = IPARSER_MUL_VP;
        }
        else if (node->l->type == IPARSER_MUL_VP &&
                 node->r->type == IPARSER_NUMBER)
        {
            int v = IPARSER_EVAL_L(node) * ((struct iparser_number*)(node->r))->value;
            IPARSER_MOVEUP_L(node, v);
            node->type = IPARSER_MUL_VP;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_MUL)
        {
            if (node->r->l->type == IPARSER_NUMBER)
            { // #l * (#rl * node_rr) -> (#l * #rl) * node_rr, same type
                int v = ((struct iparser_number*)(node->l))->value
                    *   ((struct iparser_number*)(node->r->l))->value;
                node->r = node->r->r;
                ((struct iparser_number*)(node->l))->value = v;
            }
            else if (node->r->r->type == IPARSER_NUMBER)
            { // #l * (node_rl * #rr) -> (#l * #rr) * node_rl, same type
                int v = ((struct iparser_number*)(node->l))->value
                    *   ((struct iparser_number*)(node->r->r))->value;
                node->r = node->r->l;
                ((struct iparser_number*)(node->l))->value = v;
            }
        }
        else if (node->l->type == IPARSER_MUL &&
                 node->r->type == IPARSER_NUMBER)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // (#ll * node_lr) * #r -> nodel_lr * (#ll * #r), same type
                int v = ((struct iparser_number*)(node->l->l))->value
                    *   ((struct iparser_number*)(node->r))->value;
                node->l = node->l->r;
                ((struct iparser_number*)(node->r))->value = v;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // (node_ll * #lr) * #r -> node_ll + (#lr * #r), same type
                int v = ((struct iparser_number*)(node->l->r))->value
                    *   ((struct iparser_number*)(node->r))->value;
                node->l = node->l->l;
                ((struct iparser_number*)(node->r))->value = v;
            }
        }
        break;
    case IPARSER_DIV:
    case IPARSER_DIV_PP:
        iparser_ast_optimize(node->l);
        iparser_ast_optimize(node->r);
        if (node->l->type == IPARSER_NUMBER &&
            node->r->type == IPARSER_NUMBER)
        {
            int v = ((struct iparser_number*)(node->l))->value
                /   ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        else if (node->l->type == IPARSER_NUMBER &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.v = ((struct iparser_number*)(node->l))->value;
            node->rip   = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_DIV_VP;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_NUMBER)
        {
            node->lvp.v = ((struct iparser_number*)(node->r))->value;
            node->rip   = ((struct iparser_symbol*)(node->l))->ip;
            node->r = node->l;
            node->type = IPARSER_DIV_PV;
        }
        else if (node->l->type == IPARSER_SYMBOL &&
                 node->r->type == IPARSER_SYMBOL)
        {
            node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
            node->rip    = ((struct iparser_symbol*)(node->r))->ip;
            node->type = IPARSER_DIV_PP;
        }
        break;
    case IPARSER_NEG:
        iparser_ast_optimize(node->l);
        if (node->l->type == IPARSER_NUMBER)
        {
            int v = -((struct iparser_number*)(node->l))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        else if (node->l->type == IPARSER_SYMBOL)
        {
            node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
            node->type = IPARSER_NEG_P;
        }
        else if (node->l->type == IPARSER_ADD_VP)
        {
            IPARSER_NEG_MOVEUP(node);
            node->type = IPARSER_SUB_VP;
        }
        else if (node->l->type == IPARSER_SUB_VP)
        {
            IPARSER_NEG_MOVEUP(node);
            node->type = IPARSER_ADD_VP;
        }
        else if (node->l->type == IPARSER_MUL_VP)
        {
            IPARSER_NEG_MOVEUP(node);
            node->type = IPARSER_MUL_VP;
        }
        else if (node->l->type == IPARSER_ADD)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // -(#ll + node_lr) -> -#ll - node_lr
                node->r = node->l->r;
                ((struct iparser_number*)(node->l))->value =
                    -((struct iparser_number*)(node->l->l))->value;
                node->l->type = IPARSER_NUMBER;
                node->type = IPARSER_SUB;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // -(node_ll + #lr) -> -#lr - node_ll
                node->r = node->l->l;
                ((struct iparser_number*)(node->l))->value =
                    -((struct iparser_number*)(node->l->r))->value;
                node->l->type = IPARSER_NUMBER;
                node->type = IPARSER_SUB;
            }
        }
        else if (node->l->type == IPARSER_SUB)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // -(#ll - node_lr) -> -#ll + node_lr
                node->r = node->l->r;
                ((struct iparser_number*)(node->l))->value =
                    -((struct iparser_number*)(node->l->l))->value;
                node->l->type = IPARSER_NUMBER;
                node->type = IPARSER_ADD;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // -(node_ll - #lr) -> #lr - node_ll
                node->r = node->l->l;
                ((struct iparser_number*)(node->l))->value =
                    ((struct iparser_number*)(node->l->r))->value;
                node->l->type = IPARSER_NUMBER;
                node->type = IPARSER_SUB;
            }
        }
        else if (node->l->type == IPARSER_MUL)
        {
            if (node->l->l->type == IPARSER_NUMBER)
            { // -(#ll * node_lr) -> -#ll * node_lr
                node->r = node->l->r;
                ((struct iparser_number*)(node->l))->value =
                    -((struct iparser_number*)(node->l->l))->value;
                node->l->type = IPARSER_NUMBER;
                node->type = IPARSER_MUL;
            }
            else if (node->l->r->type == IPARSER_NUMBER)
            { // -(node_ll * #lr) -> -#lr * node_ll
                node->r = node->l->l;
                ((struct iparser_number*)(node->l))->value =
                    -((struct iparser_number*)(node->l->r))->value;
                node->l->type = IPARSER_NUMBER;
                node->type = IPARSER_MUL;
            }
        }
        break;
    case IPARSER_F1:
        iparser_ast_optimize(node->l);
        if (node->l->type == IPARSER_NUMBER)
        {
            int v = iparser_call_f1
                (((struct iparser_f1*)node)->ftype,
                 ((struct iparser_number*)(((struct iparser_f1*)node)->l))->value);
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_F2:
        iparser_ast_optimize(node->l);
        iparser_ast_optimize(node->r);
        if (node->l->type == IPARSER_NUMBER &&
            node->r->type == IPARSER_NUMBER)
        {
            int v = iparser_call_f2
                (((struct iparser_f2*)node)->ftype,
                 ((struct iparser_number*)(((struct iparser_f2*)node)->l))->value,
                 ((struct iparser_number*)(((struct iparser_f2*)node)->r))->value);
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_F3:
        iparser_ast_optimize(((struct iparser_f3*)node)->n1);
        iparser_ast_optimize(((struct iparser_f3*)node)->n2);
        iparser_ast_optimize(((struct iparser_f3*)node)->n3);
        if (((struct iparser_f3*)node)->n1->type == IPARSER_NUMBER &&
            ((struct iparser_f3*)node)->n2->type == IPARSER_NUMBER &&
            ((struct iparser_f3*)node)->n3->type == IPARSER_NUMBER)
        {
            int v = iparser_call_f3
                (((struct iparser_f3*)node)->ftype,
                 ((struct iparser_number*)(((struct iparser_f3*)node)->n1))->value,
                 ((struct iparser_number*)(((struct iparser_f3*)node)->n2))->value,
                 ((struct iparser_number*)(((struct iparser_f3*)node)->n3))->value);
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_ADD_VP:
        iparser_ast_optimize(node->r);
        if (node->r->type == IPARSER_NUMBER)
        {
            int v = node->lvp.v + ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_SUB_VP:
        iparser_ast_optimize(node->r);
        if (node->r->type == IPARSER_NUMBER)
        {
            int v = node->lvp.v - ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_MUL_VP:
        iparser_ast_optimize(node->r);
        if (node->r->type == IPARSER_NUMBER)
        {
            int v = node->lvp.v * ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_DIV_VP:
        iparser_ast_optimize(node->r);
        if (node->r->type == IPARSER_NUMBER)
        {
            int v = node->lvp.v / ((struct iparser_number*)(node->r))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_DIV_PV:
        iparser_ast_optimize(node->r);
        if (node->r->type == IPARSER_NUMBER)
        {
            int v = ((struct iparser_number*)(node->r))->value / node->lvp.v;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_NEG_P:
        iparser_ast_optimize(node->l);
        if (node->l->type == IPARSER_NUMBER)
        {
            int v = -((struct iparser_number*)(node->l))->value;
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = v;
        }
        break;
    case IPARSER_ASSIGN:
        iparser_ast_optimize(((struct iparser_assign*)node)->v);
        break;
    case IPARSER_LIST:
        iparser_ast_optimize(node->l);
        iparser_ast_optimize(node->r);
        break;
    default:
        amrex::Abort("iparser_ast_optimize: unknown node type " + std::to_string(node->type));
    }
}

static
void
iparser_ast_print_f1 (struct iparser_f1* f1, std::string const& space, AllPrint& printer)
{
    printer << space;
    switch (f1->ftype) {
    case IPARSER_ABS:         printer << "ABS\n";         break;
    default:
        amrex::AllPrint() << "iparser_ast_print_f1: Unknown function " << f1->ftype << "\n";
    }
    iparser_ast_print(f1->l, space+"  ", printer);
}

static
void
iparser_ast_print_f2 (struct iparser_f2* f2, std::string const& space, AllPrint& printer)
{
    printer << space;
    switch (f2->ftype) {
    case IPARSER_FLRDIV:
        printer << "FLRDIV\n";
        break;
    case IPARSER_POW:
        printer << "POW\n";
        break;
    case IPARSER_GT:
        printer << "GT\n";
        break;
    case IPARSER_LT:
        printer << "LT\n";
        break;
    case IPARSER_GEQ:
        printer << "GEQ\n";
        break;
    case IPARSER_LEQ:
        printer << "LEQ\n";
        break;
    case IPARSER_EQ:
        printer << "EQ\n";
        break;
    case IPARSER_NEQ:
        printer << "NEQ\n";
        break;
    case IPARSER_AND:
        printer << "AND\n";
        break;
    case IPARSER_OR:
        printer << "OR\n";
        break;
    case IPARSER_MIN:
        printer << "MIN\n";
        break;
    case IPARSER_MAX:
        printer << "MAX\n";
        break;
    default:
        amrex::AllPrint() << "iparser_ast_print_f2: Unknown function " << f2->ftype << "\n";
    }
    iparser_ast_print(f2->l, space+"  ", printer);
    iparser_ast_print(f2->r, space+"  ", printer);
}

static
void
iparser_ast_print_f3 (struct iparser_f3* f3, std::string const& space, AllPrint& printer)
{
    std::string const& more_space = space + "  ";
    switch (f3->ftype) {
    case IPARSER_IF:
        printer << space << "IF\n";
        break;
    default:
        amrex::AllPrint() << "iparser_ast_print_f3: Unknown function " << f3->ftype << "\n";
    }
    iparser_ast_print(f3->n1, more_space, printer);
    iparser_ast_print(f3->n2, more_space, printer);
    iparser_ast_print(f3->n3, more_space, printer);
}

void
iparser_ast_print (struct iparser_node* node, std::string const& space, AllPrint& printer)
{
    std::string const& more_space = space + "  ";
    switch (node->type)
    {
    case IPARSER_NUMBER:
        printer << space << "NUMBER: " << ((struct iparser_number*)node)->value << "\n";
        break;
    case IPARSER_SYMBOL:
        printer << space << "VARIABLE: " << ((struct iparser_symbol*)node)->name << "\n";
        break;
    case IPARSER_ADD:
        printer << space << "ADD\n";
        iparser_ast_print(node->l, more_space, printer);
        iparser_ast_print(node->r, more_space, printer);
        break;
    case IPARSER_SUB:
        printer << space << "SUB\n";
        iparser_ast_print(node->l, more_space, printer);
        iparser_ast_print(node->r, more_space, printer);
        break;
    case IPARSER_MUL:
        printer << space << "MUL\n";
        iparser_ast_print(node->l, more_space, printer);
        iparser_ast_print(node->r, more_space, printer);
        break;
    case IPARSER_DIV:
        printer << space << "DIV\n";
        iparser_ast_print(node->l, more_space, printer);
        iparser_ast_print(node->r, more_space, printer);
        break;
    case IPARSER_NEG:
        printer << space << "NEG\n";
        iparser_ast_print(node->l, more_space, printer);
        break;
    case IPARSER_F1:
        iparser_ast_print_f1((struct iparser_f1*)node, space, printer);
        break;
    case IPARSER_F2:
        iparser_ast_print_f2((struct iparser_f2*)node, space, printer);
        break;
    case IPARSER_F3:
        iparser_ast_print_f3((struct iparser_f3*)node, space, printer);
        break;
    case IPARSER_ASSIGN:
        printer << space << "=: " << ((struct iparser_assign*)node)->s->name << " =\n";
        iparser_ast_print(((struct iparser_assign*)node)->v, more_space, printer);
        break;
    case IPARSER_LIST:
        printer << space <<"LIST\n";
        iparser_ast_print(node->l, more_space, printer);
        iparser_ast_print(node->r, more_space, printer);
        break;
    case IPARSER_ADD_VP:
        printer << space << "ADD: " << node->lvp.v << " "
                << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_SUB_VP:
        printer << space << "SUB: " << node->lvp.v << " "
                << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_MUL_VP:
        printer << space << "MUL: " << node->lvp.v << " "
                << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_DIV_VP:
        printer << space << "DIV: " << node->lvp.v << " "
                << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_DIV_PV:
        printer << space << "DIV: " << ((struct iparser_symbol*)(node->r))->name
                << " " << node->lvp.v << "\n";
        break;
    case IPARSER_NEG_P:
        printer << space << "NEG: " << ((struct iparser_symbol*)(node->l))->name << "\n";
        break;
    case IPARSER_ADD_PP:
        printer << space << "ADD: " << ((struct iparser_symbol*)(node->l))->name
                << "  "             << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_SUB_PP:
        printer << space << "SUB: " << ((struct iparser_symbol*)(node->l))->name
                << "  "             << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_MUL_PP:
        printer << space << "MUL: " << ((struct iparser_symbol*)(node->l))->name
                << "  "             << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    case IPARSER_DIV_PP:
        printer << space << "DIV: " << ((struct iparser_symbol*)(node->l))->name
                << "  "             << ((struct iparser_symbol*)(node->r))->name << "\n";
        break;
    default:
        amrex::Abort("iparser_ast_print: unknown node type " + std::to_string(node->type));
    }
}

int
iparser_ast_depth (struct iparser_node* node)
{
    switch (node->type)
    {
    case IPARSER_NUMBER:
    case IPARSER_SYMBOL:
        return 1;
    case IPARSER_ADD:
    case IPARSER_SUB:
    case IPARSER_MUL:
    case IPARSER_DIV:
    case IPARSER_LIST:
    {
        int d1 = iparser_ast_depth(node->l);
        int d2 = iparser_ast_depth(node->r);
        return std::max(d1,d2)+1;
    }
    case IPARSER_NEG:
        return iparser_ast_depth(node->l)+1;
    case IPARSER_F1:
        return iparser_ast_depth(((struct iparser_f1*)node)->l) + 1;
    case IPARSER_F2:
    {
        int d1 = iparser_ast_depth(((struct iparser_f2*)node)->l);
        int d2 = iparser_ast_depth(((struct iparser_f2*)node)->r);
        return std::max(d1,d2)+1;
    }
    case IPARSER_F3:
    {
        int d1 = iparser_ast_depth(((struct iparser_f3*)node)->n1);
        int d2 = iparser_ast_depth(((struct iparser_f3*)node)->n2);
        int d3 = iparser_ast_depth(((struct iparser_f3*)node)->n3);
        return std::max({d1,d2,d3})+1;
    }
    case IPARSER_ASSIGN:
    {
        int d = iparser_ast_depth(((struct iparser_assign*)node)->v);
        return d+1;
    }
    case IPARSER_ADD_VP:
    case IPARSER_SUB_VP:
    case IPARSER_MUL_VP:
    case IPARSER_DIV_VP:
    case IPARSER_DIV_PV:
    case IPARSER_NEG_P:
    case IPARSER_ADD_PP:
    case IPARSER_SUB_PP:
    case IPARSER_MUL_PP:
    case IPARSER_DIV_PP:
        return 1;
    default:
        amrex::Abort("iparser_ast_print: unknown node type " + std::to_string(node->type));
        return 0;
    }
}

void
iparser_ast_regvar (struct iparser_node* node, char const* name, int i)
{
    switch (node->type)
    {
    case IPARSER_NUMBER:
        break;
    case IPARSER_SYMBOL:
        if (std::strcmp(name, ((struct iparser_symbol*)node)->name) == 0) {
            ((struct iparser_symbol*)node)->ip = i;
        }
        break;
    case IPARSER_ADD:
    case IPARSER_SUB:
    case IPARSER_MUL:
    case IPARSER_DIV:
    case IPARSER_LIST:
        iparser_ast_regvar(node->l, name, i);
        iparser_ast_regvar(node->r, name, i);
        break;
    case IPARSER_NEG:
        iparser_ast_regvar(node->l, name, i);
        break;
    case IPARSER_F1:
        iparser_ast_regvar(((struct iparser_f1*)node)->l, name, i);
        break;
    case IPARSER_F2:
        iparser_ast_regvar(((struct iparser_f2*)node)->l, name, i);
        iparser_ast_regvar(((struct iparser_f2*)node)->r, name, i);
        break;
    case IPARSER_F3:
        iparser_ast_regvar(((struct iparser_f3*)node)->n1, name, i);
        iparser_ast_regvar(((struct iparser_f3*)node)->n2, name, i);
        iparser_ast_regvar(((struct iparser_f3*)node)->n3, name, i);
        break;
    case IPARSER_ASSIGN:
        iparser_ast_regvar(((struct iparser_assign*)node)->v, name, i);
        break;
    case IPARSER_ADD_VP:
    case IPARSER_SUB_VP:
    case IPARSER_MUL_VP:
    case IPARSER_DIV_VP:
    case IPARSER_DIV_PV:
        iparser_ast_regvar(node->r, name, i);
        node->rip = ((struct iparser_symbol*)(node->r))->ip;
        break;
    case IPARSER_NEG_P:
        iparser_ast_regvar(node->l, name, i);
        node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
        break;
    case IPARSER_ADD_PP:
    case IPARSER_SUB_PP:
    case IPARSER_MUL_PP:
    case IPARSER_DIV_PP:
        iparser_ast_regvar(node->l, name, i);
        iparser_ast_regvar(node->r, name, i);
        node->lvp.ip = ((struct iparser_symbol*)(node->l))->ip;
        node->rip = ((struct iparser_symbol*)(node->r))->ip;
        break;
    default:
        amrex::AllPrint() << "iparser_ast_regvar: unknown node type " << node->type << "\n";
        amrex::Abort();
    }
}

void iparser_ast_setconst (struct iparser_node* node, char const* name, int c)
{
    switch (node->type)
    {
    case IPARSER_NUMBER:
        break;
    case IPARSER_SYMBOL:
        if (std::strcmp(name, ((struct iparser_symbol*)node)->name) == 0) {
            ((struct iparser_number*)node)->type = IPARSER_NUMBER;
            ((struct iparser_number*)node)->value = c;
        }
        break;
    case IPARSER_ADD:
    case IPARSER_SUB:
    case IPARSER_MUL:
    case IPARSER_DIV:
    case IPARSER_ADD_PP:
    case IPARSER_SUB_PP:
    case IPARSER_MUL_PP:
    case IPARSER_DIV_PP:
    case IPARSER_LIST:
        iparser_ast_setconst(node->l, name, c);
        iparser_ast_setconst(node->r, name, c);
        break;
    case IPARSER_NEG:
    case IPARSER_NEG_P:
        iparser_ast_setconst(node->l, name, c);
        break;
    case IPARSER_F1:
        iparser_ast_setconst(((struct iparser_f1*)node)->l, name, c);
        break;
    case IPARSER_F2:
        iparser_ast_setconst(((struct iparser_f2*)node)->l, name, c);
        iparser_ast_setconst(((struct iparser_f2*)node)->r, name, c);
        break;
    case IPARSER_F3:
        iparser_ast_setconst(((struct iparser_f3*)node)->n1, name, c);
        iparser_ast_setconst(((struct iparser_f3*)node)->n2, name, c);
        iparser_ast_setconst(((struct iparser_f3*)node)->n3, name, c);
        break;
    case IPARSER_ASSIGN:
        iparser_ast_setconst(((struct iparser_assign*)node)->v, name, c);
        break;
    case IPARSER_ADD_VP:
    case IPARSER_SUB_VP:
    case IPARSER_MUL_VP:
    case IPARSER_DIV_VP:
    case IPARSER_DIV_PV:
        iparser_ast_setconst(node->r, name, c);
        break;
    default:
        amrex::Abort("iparser_ast_setconst: unknown node type " + std::to_string(node->type));
    }
}

void iparser_ast_get_symbols (struct iparser_node* node, std::set<std::string>& symbols,
                              std::set<std::string>& local_symbols)
{
    switch (node->type)
    {
    case IPARSER_NUMBER:
        break;
    case IPARSER_SYMBOL:
        symbols.emplace(((struct iparser_symbol*)node)->name);
        break;
    case IPARSER_ADD:
    case IPARSER_SUB:
    case IPARSER_MUL:
    case IPARSER_DIV:
    case IPARSER_ADD_PP:
    case IPARSER_SUB_PP:
    case IPARSER_MUL_PP:
    case IPARSER_DIV_PP:
    case IPARSER_LIST:
        iparser_ast_get_symbols(node->l, symbols, local_symbols);
        iparser_ast_get_symbols(node->r, symbols, local_symbols);
        break;
    case IPARSER_NEG:
    case IPARSER_NEG_P:
        iparser_ast_get_symbols(node->l, symbols, local_symbols);
        break;
    case IPARSER_F1:
        iparser_ast_get_symbols(((struct iparser_f1*)node)->l, symbols, local_symbols);
        break;
    case IPARSER_F2:
        iparser_ast_get_symbols(((struct iparser_f2*)node)->l, symbols, local_symbols);
        iparser_ast_get_symbols(((struct iparser_f2*)node)->r, symbols, local_symbols);
        break;
    case IPARSER_F3:
        iparser_ast_get_symbols(((struct iparser_f3*)node)->n1, symbols, local_symbols);
        iparser_ast_get_symbols(((struct iparser_f3*)node)->n2, symbols, local_symbols);
        iparser_ast_get_symbols(((struct iparser_f3*)node)->n3, symbols, local_symbols);
        break;
    case IPARSER_ASSIGN:
        local_symbols.emplace(((struct iparser_assign*)node)->s->name);
        iparser_ast_get_symbols(((struct iparser_assign*)node)->v, symbols, local_symbols);
        break;
    case IPARSER_ADD_VP:
    case IPARSER_SUB_VP:
    case IPARSER_MUL_VP:
    case IPARSER_DIV_VP:
    case IPARSER_DIV_PV:
        iparser_ast_get_symbols(node->r, symbols, local_symbols);
        break;
    default:
        amrex::Abort("iparser_ast_get_symbols: unknown node type " + std::to_string(node->type));
    }
}

void
iparser_regvar (struct amrex_iparser* iparser, char const* name, int i)
{
    iparser_ast_regvar(iparser->ast, name, i);
}

void
iparser_setconst (struct amrex_iparser* iparser, char const* name, int c)
{
    iparser_ast_setconst(iparser->ast, name, c);
    iparser_ast_optimize(iparser->ast);
}

void
iparser_print (struct amrex_iparser* iparser)
{
    amrex::AllPrint printer{};
    iparser_ast_print(iparser->ast, std::string("  "), printer);
}

std::set<std::string>
iparser_get_symbols (struct amrex_iparser* iparser)
{
    std::set<std::string> symbols;
    std::set<std::string> local_symbols;
    iparser_ast_get_symbols(iparser->ast, symbols, local_symbols);
    for (auto const& ls : local_symbols) {
        symbols.erase(ls);
    }
    return symbols;
}

int
iparser_depth (struct amrex_iparser* iparser)
{
    return iparser_ast_depth(iparser->ast);
}

}
