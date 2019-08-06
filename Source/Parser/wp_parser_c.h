#ifndef WP_PARSER_C_H_
#define WP_PARSER_C_H_

#include "wp_parser_y.h"
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>

#ifdef __cplusplus
extern "C" {
#endif

    struct wp_parser* wp_c_parser_new (char const* function_body);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

#include <set>
#include <string>

AMREX_GPU_HOST_DEVICE
inline double
wp_ast_eval (struct wp_node* node)
{
    double result;

#ifdef AMREX_DEVICE_COMPILE
    extern __shared__ double extern_xyz[];
    int tid = threadIdx.x + threadIdx.y*blockDim.x + threadIdx.z*(blockDim.x*blockDim.y);
    double* x = extern_xyz + tid*3;
#endif

    switch (node->type)
    {
    case WP_NUMBER:
    {
        result = ((struct wp_number*)node)->value;
        break;
    }
    case WP_SYMBOL:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i =((struct wp_symbol*)node)->ip.i;
        result = x[i];
#else
        result = *(((struct wp_symbol*)node)->ip.p);
#endif
        break;
    }
    case WP_ADD:
    {
        result = wp_ast_eval(node->l) + wp_ast_eval(node->r);
        break;
    }
    case WP_SUB:
    {
        result = wp_ast_eval(node->l) - wp_ast_eval(node->r);
        break;
    }
    case WP_MUL:
    {
        result = wp_ast_eval(node->l) * wp_ast_eval(node->r);
        break;
    }
    case WP_DIV:
    {
        result = wp_ast_eval(node->l) / wp_ast_eval(node->r);
        break;
    }
    case WP_NEG:
    {
        result = -wp_ast_eval(node->l);
        break;
    }
    case WP_F1:
    {
        result = wp_call_f1(((struct wp_f1*)node)->ftype,
                wp_ast_eval(((struct wp_f1*)node)->l));
        break;
    }
    case WP_F2:
    {
        result = wp_call_f2(((struct wp_f2*)node)->ftype,
                wp_ast_eval(((struct wp_f2*)node)->l),
                wp_ast_eval(((struct wp_f2*)node)->r));
        break;
    }
    case WP_ADD_VP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v + x[i];
#else
        result = node->lvp.v + *(node->rip.p);
#endif
        break;
    }
    case WP_ADD_PP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] + x[j];
#else
        result = *(node->lvp.ip.p) + *(node->rip.p);
#endif
        break;
    }
    case WP_SUB_VP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v - x[i];
#else
        result = node->lvp.v - *(node->rip.p);
#endif
        break;
    }
    case WP_SUB_PP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] - x[j];
#else
        result = *(node->lvp.ip.p) - *(node->rip.p);
#endif
        break;
    }
    case WP_MUL_VP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v * x[i];
#else
        result = node->lvp.v * *(node->rip.p);
#endif
        break;
    }
    case WP_MUL_PP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] * x[j];
#else
        result = *(node->lvp.ip.p) * *(node->rip.p);
#endif
        break;
    }
    case WP_DIV_VP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = node->lvp.v / x[i];
#else
        result = node->lvp.v / *(node->rip.p);
#endif
        break;
    }
    case WP_DIV_PP:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->lvp.ip.i;
        int j = node->rip.i;
        result = x[i] / x[j];
#else
        result = *(node->lvp.ip.p) / *(node->rip.p);
#endif
        break;
    }
    case WP_NEG_P:
    {
#ifdef AMREX_DEVICE_COMPILE
        int i = node->rip.i;
        result = -x[i];
#else
        result = -*(node->lvp.ip.p);
#endif
        break;
    }
    default:
        yyerror("wp_ast_eval: unknown node type %d\n", node->type);
    }

    return result;
}

inline
void
wp_ast_get_symbols (struct wp_node* node, std::set<std::string>& symbols)
{
    switch (node->type)
    {
    case WP_NUMBER:
        break;
    case WP_SYMBOL:
        symbols.emplace(((struct wp_symbol*)node)->name);
        break;
    case WP_ADD:
    case WP_SUB:
    case WP_MUL:
    case WP_DIV:
    case WP_ADD_PP:
    case WP_SUB_PP:
    case WP_MUL_PP:
    case WP_DIV_PP:
        wp_ast_get_symbols(node->l, symbols);
        wp_ast_get_symbols(node->r, symbols);
        break;
    case WP_NEG:
    case WP_NEG_P:
        wp_ast_get_symbols(node->l, symbols);
        break;
    case WP_F1:
        wp_ast_get_symbols(((struct wp_f1*)node)->l, symbols);
        break;
    case WP_F2:
        wp_ast_get_symbols(((struct wp_f2*)node)->l, symbols);
        wp_ast_get_symbols(((struct wp_f2*)node)->r, symbols);
        break;
    case WP_ADD_VP:
    case WP_SUB_VP:
    case WP_MUL_VP:
    case WP_DIV_VP:
        wp_ast_get_symbols(node->r, symbols);
        break;
    default:
        yyerror("wp_ast_get_symbols: unknown node type %d\n", node->type);
    }
}

#endif

#endif
