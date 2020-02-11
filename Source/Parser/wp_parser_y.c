#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "wp_parser_y.h"
#include "wp_parser.tab.h"

static struct wp_node* wp_root = NULL;

/* This is called by a bison rule to store the original AST in a
 * static variable.  Accessing this directly is not thread safe.  So
 * this will be duplicated later for each thread.
 */
void
wp_defexpr (struct wp_node* body)
{
    wp_root = body;
}

struct wp_node*
wp_newnumber (amrex_real d)
{
    struct wp_number* r = (struct wp_number*) malloc(sizeof(struct wp_number));
    r->type = WP_NUMBER;
    r->value = d;
    return (struct wp_node*) r;
}

struct wp_symbol*
wp_makesymbol (char* name)
{
    struct wp_symbol* symbol = (struct wp_symbol*) malloc(sizeof(struct wp_symbol));
    symbol->type = WP_SYMBOL;
    symbol->name = strdup(name);
    symbol->ip.p = NULL;
    return symbol;
}

struct wp_node*
wp_newsymbol (struct wp_symbol* symbol)
{
    return (struct wp_node*) symbol;
}

struct wp_node*
wp_newnode (enum wp_node_t type, struct wp_node* l, struct wp_node* r)
{
    struct wp_node* tmp = (struct wp_node*) malloc(sizeof(struct wp_node));
    tmp->type = type;
    tmp->l = l;
    tmp->r = r;
    return tmp;
}

struct wp_node*
wp_newf1 (enum wp_f1_t ftype, struct wp_node* l)
{
    struct wp_f1* tmp = (struct wp_f1*) malloc(sizeof(struct wp_f1));
    tmp->type = WP_F1;
    tmp->l = l;
    tmp->ftype = ftype;
    return (struct wp_node*) tmp;
}

struct wp_node*
wp_newf2 (enum wp_f2_t ftype, struct wp_node* l, struct wp_node* r)
{
    struct wp_f2* tmp = (struct wp_f2*) malloc(sizeof(struct wp_f2));
    tmp->type = WP_F2;
    tmp->l = l;
    tmp->r = r;
    tmp->ftype = ftype;
    return (struct wp_node*) tmp;
}

AMREX_GPU_HOST_DEVICE
void
yyerror (char const *s, ...)
{
    va_list vl;
    va_start(vl, s);
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    printf(s,"\n");
    assert(0);
#else
    vfprintf(stderr, s, vl);
    fprintf(stderr, "\n");
#endif
    va_end(vl);
}

/*******************************************************************/

struct wp_parser*
wp_parser_new (void)
{
    struct wp_parser* my_parser = (struct wp_parser*) malloc(sizeof(struct wp_parser));

    my_parser->sz_mempool = wp_ast_size(wp_root);
    my_parser->p_root = malloc(my_parser->sz_mempool);
    my_parser->p_free = my_parser->p_root;

    my_parser->ast = wp_parser_ast_dup(my_parser, wp_root,1); /* 1: free the source wp_root */

    if ((char*)my_parser->p_root + my_parser->sz_mempool != (char*)my_parser->p_free) {
        yyerror("wp_parser_new: error in memory size");
        exit(1);
    }

    wp_ast_optimize(my_parser->ast);

    return my_parser;
}

void
wp_parser_delete (struct wp_parser* parser)
{
    free(parser->p_root);
    free(parser);
}

static size_t
wp_aligned_size (size_t N)
{
    const unsigned int align_size = 16;
    size_t x = N + (align_size-1);
    x -= x & (align_size-1);
    return x;
}

static
void*
wp_parser_allocate (struct wp_parser* my_parser, size_t N)
{
    void* r = my_parser->p_free;
    my_parser->p_free = (char*)r + wp_aligned_size(N);
    return r;
}

struct wp_parser*
wp_parser_dup (struct wp_parser* source)
{
    struct wp_parser* dest = (struct wp_parser*) malloc(sizeof(struct wp_parser));
    dest->sz_mempool = source->sz_mempool;
    dest->p_root = malloc(dest->sz_mempool);
    dest->p_free = dest->p_root;

    dest->ast = wp_parser_ast_dup(dest, source->ast, 0); /* 0: don't free the source */

    return dest;
}

AMREX_GPU_HOST_DEVICE
amrex_real
wp_call_f1 (enum wp_f1_t type, amrex_real a)
{
    switch (type) {
    case WP_SQRT:        return sqrt(a);
    case WP_EXP:         return exp(a);
    case WP_LOG:         return log(a);
    case WP_LOG10:       return log10(a);
    case WP_SIN:         return sin(a);
    case WP_COS:         return cos(a);
    case WP_TAN:         return tan(a);
    case WP_ASIN:        return asin(a);
    case WP_ACOS:        return acos(a);
    case WP_ATAN:        return atan(a);
    case WP_SINH:        return sinh(a);
    case WP_COSH:        return cosh(a);
    case WP_TANH:        return tanh(a);
    case WP_ABS:         return fabs(a);
    case WP_POW_M3:      return 1.0/(a*a*a);
    case WP_POW_M2:      return 1.0/(a*a);
    case WP_POW_M1:      return 1.0/a;
    case WP_POW_P1:      return a;
    case WP_POW_P2:      return a*a;
    case WP_POW_P3:      return a*a*a;
    default:
        yyerror("wp_call_f1: Unknow function %d", type);
        return 0.0;
    }
}

AMREX_GPU_HOST_DEVICE
amrex_real
wp_call_f2 (enum wp_f2_t type, amrex_real a, amrex_real b)
{
    switch (type) {
    case WP_POW:
        return pow(a,b);
    case WP_GT:
        return (a > b) ? 1.0 : 0.0;
    case WP_LT:
        return (a < b) ? 1.0 : 0.0;
    case WP_GEQ:
        return (a >= b) ? 1.0 : 0.0;
    case WP_LEQ:
        return (a <= b) ? 1.0 : 0.0;
    case WP_EQ:
        return (a == b) ? 1.0 : 0.0;
    case WP_NEQ:
        return (a != b) ? 1.0 : 0.0;
    case WP_AND:
        return (a && b) ? 1.0 : 0.0;
    case WP_OR:
        return (a || b) ? 1.0 : 0.0;
    case WP_HEAVISIDE:
        return (a < 0.0) ? 0.0 : ((a > 0.0) ? 1.0 : b);
    case WP_MIN:
        return (a < b) ? a : b;
    case WP_MAX:
        return (a > b) ? a : b;
    default:
        yyerror("wp_call_f2: Unknow function %d", type);
        return 0.0;
    }
}

size_t
wp_ast_size (struct wp_node* node)
{
    size_t result;

    switch (node->type)
    {
    case WP_NUMBER:
        result = wp_aligned_size(sizeof(struct wp_number));
        break;
    case WP_SYMBOL:
        result = wp_aligned_size(sizeof(struct wp_symbol))
            + wp_aligned_size(strlen(((struct wp_symbol*)node)->name)+1);
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->l) + wp_ast_size(node->r);
        break;
    case WP_NEG:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->l);
        break;
    case WP_F1:
        result = wp_aligned_size(sizeof(struct wp_f1))
            + wp_ast_size(node->l);
        break;
    case WP_F2:
        result = wp_aligned_size(sizeof(struct wp_f2))
            + wp_ast_size(node->l) + wp_ast_size(node->r);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->r);
        break;
    case WP_NEG_P:
        result = wp_aligned_size(sizeof(struct wp_node))
            + wp_ast_size(node->l);
        break;
    default:
        yyerror("wp_ast_size: unknown node type %d\n", node->type);
        exit(1);
    }

    return result;
}

struct wp_node*
wp_parser_ast_dup (struct wp_parser* my_parser, struct wp_node* node, int move)
{
    void* result;

    switch (node->type)
    {
    case WP_NUMBER:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_number));
        memcpy(result, node                  , sizeof(struct wp_number));
        break;
    case WP_SYMBOL:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_symbol));
        memcpy(result, node                  , sizeof(struct wp_symbol));
        ((struct wp_symbol*)result)->name = (char*) wp_parser_allocate
            (my_parser, strlen(((struct wp_symbol*)node)->name)+1);
        strcpy(((struct wp_symbol*)result)->name,
               ((struct wp_symbol*)node  )->name);
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->l = wp_parser_ast_dup(my_parser, node->l, move);
        ((struct wp_node*)result)->r = wp_parser_ast_dup(my_parser, node->r, move);
        break;
    case WP_NEG:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->l = wp_parser_ast_dup(my_parser, node->l, move);
        ((struct wp_node*)result)->r = NULL;
        break;
    case WP_F1:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_f1));
        memcpy(result, node                  , sizeof(struct wp_f1));
        ((struct wp_f1*)result)->l = wp_parser_ast_dup(my_parser, ((struct wp_f1*)node)->l, move);
        break;
    case WP_F2:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_f2));
        memcpy(result, node                  , sizeof(struct wp_f2));
        ((struct wp_f2*)result)->l = wp_parser_ast_dup(my_parser, ((struct wp_f2*)node)->l, move);
        ((struct wp_f2*)result)->r = wp_parser_ast_dup(my_parser, ((struct wp_f2*)node)->r, move);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->r = wp_parser_ast_dup(my_parser, node->r, move);
        break;
    case WP_NEG_P:
        result = wp_parser_allocate(my_parser, sizeof(struct wp_node));
        memcpy(result, node                  , sizeof(struct wp_node));
        ((struct wp_node*)result)->l = wp_parser_ast_dup(my_parser, node->l, move);
        break;
    default:
        yyerror("wp_ast_dup: unknown node type %d\n", node->type);
        exit(1);
    }
    if (move) {
        /* Note that we only do this on the original AST.  We do not
         * need to call free for AST stored in wp_parser because the
         * memory is not allocated with malloc directly.
         */
        if (node->type == WP_SYMBOL) {
            free(((struct wp_symbol*)node)->name);
        }
        free((void*)node);
    }
    return (struct wp_node*)result;
}

#define WP_MOVEUP_R(node, v) \
    struct wp_node* n = node->r->r; \
    amrex_real* p = node->r->rip.p; \
    node->r = n; \
    node->lvp.v = v; \
    node->rip.p = p;
#define WP_MOVEUP_L(node, v) \
    struct wp_node* n = node->l->r; \
    amrex_real* p = node->l->rip.p; \
    node->r = n; \
    node->lvp.v = v; \
    node->rip.p = p;
#define WP_EVAL_R(node) node->r->lvp.v
#define WP_EVAL_L(node) node->l->lvp.v

#define WP_NEG_MOVEUP(node) \
    node->r = node->l->r; \
    node->lvp.v = -node->l->lvp.v; \
    node->rip.p = node->l->rip.p;

void
wp_ast_optimize (struct wp_node* node)
{
    /* No need to free memory because we only call this on ASTs in
     * wp_parser that are allocated from the memory pool.
     */
    switch (node->type)
    {
    case WP_NUMBER:
    case WP_SYMBOL:
        break;
    case WP_ADD:
    case WP_ADD_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                +      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = ((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_ADD_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_ADD_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value + WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SUB_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value + WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_ADD_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) + ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SUB_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) + ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_SUB_VP;
        }
        break;
    case WP_SUB:
    case WP_SUB_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                -      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = -((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_SUB_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_ADD_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value - WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SUB_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value - WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_ADD_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) - ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_SUB_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) - ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_SUB_VP;
        }
        break;
    case WP_MUL:
    case WP_MUL_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                *      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = ((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_MUL_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_MUL_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value * WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_DIV_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value * WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_DIV_VP;
        }
        else if (node->l->type == WP_MUL_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) * ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_DIV_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) * ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_DIV_VP;
        }
        break;
    case WP_DIV:
    case WP_DIV_PP:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value
                /      ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.v = ((struct wp_number*)(node->l))->value;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_DIV_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_NUMBER)
        {
            node->lvp.v = 1./((struct wp_number*)(node->r))->value;
            node->rip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->r = node->l;
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_SYMBOL &&
                 node->r->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
            node->type = WP_DIV_PP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_MUL_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value / WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_DIV_VP;
        }
        else if (node->l->type == WP_NUMBER &&
                 node->r->type == WP_DIV_VP)
        {
            amrex_real v = ((struct wp_number*)(node->l))->value / WP_EVAL_R(node);
            WP_MOVEUP_R(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_MUL_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) / ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_DIV_VP &&
                 node->r->type == WP_NUMBER)
        {
            amrex_real v = WP_EVAL_L(node) / ((struct wp_number*)(node->r))->value;
            WP_MOVEUP_L(node, v);
            node->type = WP_DIV_VP;
        }
        break;
    case WP_NEG:
        wp_ast_optimize(node->l);
        if (node->l->type == WP_NUMBER)
        {
            amrex_real v = -((struct wp_number*)(node->l))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->l->type == WP_SYMBOL)
        {
            node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
            node->type = WP_NEG_P;
        }
        else if (node->l->type == WP_ADD_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_SUB_VP;
        }
        else if (node->l->type == WP_SUB_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_ADD_VP;
        }
        else if (node->l->type == WP_MUL_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_MUL_VP;
        }
        else if (node->l->type == WP_DIV_VP)
        {
            WP_NEG_MOVEUP(node);
            node->type = WP_DIV_VP;
        }
        break;
    case WP_F1:
        wp_ast_optimize(node->l);
        if (node->l->type == WP_NUMBER)
        {
            amrex_real v = wp_call_f1
                (((struct wp_f1*)node)->ftype,
                 ((struct wp_number*)(((struct wp_f1*)node)->l))->value);
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_F2:
        wp_ast_optimize(node->l);
        wp_ast_optimize(node->r);
        if (node->l->type == WP_NUMBER &&
            node->r->type == WP_NUMBER)
        {
            amrex_real v = wp_call_f2
                (((struct wp_f2*)node)->ftype,
                 ((struct wp_number*)(((struct wp_f2*)node)->l))->value,
                 ((struct wp_number*)(((struct wp_f2*)node)->r))->value);
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        else if (node->r->type == WP_NUMBER && ((struct wp_f2*)node)->ftype == WP_POW)
        {
            struct wp_node* n = node->l;
            amrex_real v = ((struct wp_number*)(node->r))->value;
            if (-3.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_M3;
            } else if (-2.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_M2;
            } else if (-1.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_M1;
            } else if (0.0 == v) {
                ((struct wp_number*)node)->type = WP_NUMBER;
                ((struct wp_number*)node)->value = 1.0;
            } else if (1.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_P1;
            } else if (2.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_P2;
            } else if (3.0 == v) {
                ((struct wp_f1*)node)->type = WP_F1;
                ((struct wp_f1*)node)->l = n;
                ((struct wp_f1*)node)->ftype = WP_POW_P3;
            }
        }
        break;
    case WP_ADD_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v + ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_SUB_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v - ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_MUL_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v * ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_DIV_VP:
        wp_ast_optimize(node->r);
        if (node->r->type == WP_NUMBER)
        {
            amrex_real v = node->lvp.v / ((struct wp_number*)(node->r))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    case WP_NEG_P:
        wp_ast_optimize(node->l);
        if (node->l->type == WP_NUMBER)
        {
            amrex_real v = -((struct wp_number*)(node->l))->value;
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = v;
        }
        break;
    default:
        yyerror("wp_ast_optimize: unknown node type %d\n", node->type);
        exit(1);
    }
}

static
void
wp_ast_print_f1 (struct wp_f1* f1)
{
    wp_ast_print(f1->l);
    switch (f1->ftype) {
    case WP_SQRT:        printf("SQRT\n");        break;
    case WP_EXP:         printf("EXP\n");         break;
    case WP_LOG:         printf("LOG\n");         break;
    case WP_LOG10:       printf("LOG10\n");       break;
    case WP_SIN:         printf("SIN\n");         break;
    case WP_COS:         printf("COS\n");         break;
    case WP_TAN:         printf("TAN\n");         break;
    case WP_ASIN:        printf("ASIN\n");        break;
    case WP_ACOS:        printf("ACOS\n");        break;
    case WP_ATAN:        printf("ATAN\n");        break;
    case WP_SINH:        printf("SINH\n");        break;
    case WP_COSH:        printf("COSH\n");        break;
    case WP_TANH:        printf("TANH\n");        break;
    case WP_ABS:         printf("ABS\n");         break;
    case WP_POW_M3:      printf("POW(,-3)\n");    break;
    case WP_POW_M2:      printf("POW(,-2)\n");    break;
    case WP_POW_M1:      printf("POW(,-1)\n");    break;
    case WP_POW_P1:      printf("POW(,1)\n");     break;
    case WP_POW_P2:      printf("POW(,2)\n");     break;
    case WP_POW_P3:      printf("POW(,3)\n");     break;
    default:
        yyerror("wp_ast+print_f1: Unknow function %d", f1->ftype);
    }
}

static
void
wp_ast_print_f2 (struct wp_f2* f2)
{
    wp_ast_print(f2->l);
    wp_ast_print(f2->r);
    switch (f2->ftype) {
    case WP_POW:
        printf("POW\n");
        break;
    case WP_GT:
        printf("GT\n");
        break;
    case WP_LT:
        printf("LT\n");
        break;
    case WP_GEQ:
        printf("GEQ\n");
        break;
    case WP_LEQ:
        printf("LEQ\n");
        break;
    case WP_EQ:
        printf("EQ\n");
        break;
    case WP_NEQ:
        printf("NEQ\n");
        break;
    case WP_AND:
        printf("AND\n");
        break;
    case WP_OR:
        printf("OR\n");
        break;
    case WP_HEAVISIDE:
        printf("HEAVISIDE\n");
        break;
    case WP_MIN:
        printf("MIN\n");
        break;
    case WP_MAX:
        printf("MAX\n");
        break;
    default:
        yyerror("wp_ast_print_f2: Unknow function %d", f2->ftype);
    }
}

void
wp_ast_print (struct wp_node* node)
{
    switch (node->type)
    {
    case WP_NUMBER:
        printf("NUMBER:  %.17g\n", ((struct wp_number*)node)->value);
        break;
    case WP_SYMBOL:
        printf("VARIABLE:  %s\n", ((struct wp_symbol*)node)->name);
        break;
    case WP_ADD:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        printf("ADD\n");
        break;
    case WP_SUB:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        printf("SUB\n");
        break;
    case WP_MUL:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        printf("MUL\n");
        break;
    case WP_DIV:
        wp_ast_print(node->l);
        wp_ast_print(node->r);
        printf("DIV\n");
        break;
    case WP_NEG:
        wp_ast_print(node->l);
        printf("NEG\n");
        break;
    case WP_F1:
        wp_ast_print_f1((struct wp_f1*)node);
        break;
    case WP_F2:
        wp_ast_print_f2((struct wp_f2*)node);
        break;
    case WP_ADD_VP:
        printf("ADD:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_SUB_VP:
        printf("SUM:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_MUL_VP:
        printf("MUL:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_DIV_VP:
        printf("DIV:  %.17g  %s\n", node->lvp.v, ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_NEG_P:
        printf("NEG:  %s\n", ((struct wp_symbol*)(node->l))->name);
        break;
    case WP_ADD_PP:
        printf("ADD:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                 ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_SUB_PP:
        printf("SUB:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                 ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_MUL_PP:
        printf("MUL:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                 ((struct wp_symbol*)(node->r))->name);
        break;
    case WP_DIV_PP:
        printf("DIV:  %s  %s\n", ((struct wp_symbol*)(node->l))->name,
                                 ((struct wp_symbol*)(node->r))->name);
        break;
    default:
        yyerror("wp_ast_print: unknown node type %d\n", node->type);
        exit(1);
    }
}

void
wp_ast_regvar (struct wp_node* node, char const* name, amrex_real* p)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        if (strcmp(name, ((struct wp_symbol*)node)->name) == 0) {
            ((struct wp_symbol*)node)->ip.p = p;
        }
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
        wp_ast_regvar(node->l, name, p);
        wp_ast_regvar(node->r, name, p);
        break;
    case WP_NEG:
        wp_ast_regvar(node->l, name, p);
        break;
    case WP_F1:
        wp_ast_regvar(node->l, name, p);
        break;
    case WP_F2:
        wp_ast_regvar(node->l, name, p);
        wp_ast_regvar(node->r, name, p);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_regvar(node->r, name, p);
        node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
        break;
    case WP_NEG_P:
        wp_ast_regvar(node->l, name, p);
        node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
        break;
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_regvar(node->l, name, p);
        wp_ast_regvar(node->r, name, p);
        node->lvp.ip.p = ((struct wp_symbol*)(node->l))->ip.p;
        node->rip.p = ((struct wp_symbol*)(node->r))->ip.p;
        break;
    default:
        yyerror("wp_ast_regvar: unknown node type %d\n", node->type);
        exit(1);
    }
}

void
wp_ast_regvar_gpu (struct wp_node* node, char const* name, int i)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        if (strcmp(name, ((struct wp_symbol*)node)->name) == 0) {
            ((struct wp_symbol*)node)->ip.i = i;
        }
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
        wp_ast_regvar_gpu(node->l, name, i);
        wp_ast_regvar_gpu(node->r, name, i);
        break;
    case WP_NEG:
        wp_ast_regvar_gpu(node->l, name, i);
        break;
    case WP_F1:
        wp_ast_regvar_gpu(node->l, name, i);
        break;
    case WP_F2:
        wp_ast_regvar_gpu(node->l, name, i);
        wp_ast_regvar_gpu(node->r, name, i);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_regvar_gpu(node->r, name, i);
        node->rip.i = ((struct wp_symbol*)(node->r))->ip.i;
        break;
    case WP_NEG_P:
        wp_ast_regvar_gpu(node->l, name, i);
        node->lvp.ip.i = ((struct wp_symbol*)(node->l))->ip.i;
        break;
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_regvar_gpu(node->l, name, i);
        wp_ast_regvar_gpu(node->r, name, i);
        node->lvp.ip.i = ((struct wp_symbol*)(node->l))->ip.i;
        node->rip.i = ((struct wp_symbol*)(node->r))->ip.i;
        break;
    default:
        yyerror("wp_ast_regvar_gpu: unknown node type %d\n", node->type);
        exit(1);
    }
}

void wp_ast_setconst (struct wp_node* node, char const* name, amrex_real c)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        if (strcmp(name, ((struct wp_symbol*)node)->name) == 0) {
            ((struct wp_number*)node)->type = WP_NUMBER;
            ((struct wp_number*)node)->value = c;
        }
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
        wp_ast_setconst(node->l, name, c);
        wp_ast_setconst(node->r, name, c);
        break;
    case WP_NEG:
        wp_ast_setconst(node->l, name, c);
        break;
    case WP_F1:
        wp_ast_setconst(node->l, name, c);
        break;
    case WP_F2:
        wp_ast_setconst(node->l, name, c);
        wp_ast_setconst(node->r, name, c);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_setconst(node->r, name, c);
        break;
    case WP_NEG_P:
        wp_ast_setconst(node->l, name, c);
        break;
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_setconst(node->l, name, c);
        wp_ast_setconst(node->r, name, c);
        break;
    default:
        yyerror("wp_ast_setconst: unknown node type %d\n", node->type);
        exit(1);
    }
}

void
wp_parser_regvar (struct wp_parser* parser, char const* name, amrex_real* p)
{
    wp_ast_regvar(parser->ast, name, p);
}

void
wp_parser_regvar_gpu (struct wp_parser* parser, char const* name, int i)
{
    wp_ast_regvar_gpu(parser->ast, name, i);
}

void
wp_parser_setconst (struct wp_parser* parser, char const* name, amrex_real c)
{
    wp_ast_setconst(parser->ast, name, c);
    wp_ast_optimize(parser->ast);
}
