#ifndef WP_PARSER_C_H_
#define WP_PARSER_C_H_

#include "wp_parser_y.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct wp_parser* wp_c_parser_new (char const* function_body);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

inline
double
wp_ast_eval (struct wp_node* node)
{
    double result;

    switch (node->type)
    {
    case WP_NUMBER:
        result = ((struct wp_number*)node)->value;
        break;
    case WP_SYMBOL:
        result = *(((struct wp_symbol*)node)->pointer);
        break;
    case WP_ADD:
        result = wp_ast_eval(node->l) + wp_ast_eval(node->r);
        break;
    case WP_SUB:
        result = wp_ast_eval(node->l) - wp_ast_eval(node->r);
        break;
    case WP_MUL:
        result = wp_ast_eval(node->l) * wp_ast_eval(node->r);
        break;
    case WP_DIV:
        result = wp_ast_eval(node->l) / wp_ast_eval(node->r);
        break;
    case WP_NEG:
        result = -wp_ast_eval(node->l);
        break;
    case WP_F1:
        result = wp_call_f1(((struct wp_f1*)node)->ftype,
                wp_ast_eval(((struct wp_f1*)node)->l));
        break;
    case WP_F2:
        result = wp_call_f2(((struct wp_f2*)node)->ftype,
                wp_ast_eval(((struct wp_f2*)node)->l),
                wp_ast_eval(((struct wp_f2*)node)->r));
        break;
    case WP_ADD_VP:
        result = node->lvp.v + *(node->rp);
        break;
    case WP_ADD_PP:
        result = *(node->lvp.p) + *(node->rp);
        break;
    case WP_SUB_VP:
        result = node->lvp.v - *(node->rp);
        break;
    case WP_SUB_PP:
        result = *(node->lvp.p) - *(node->rp);
        break;
    case WP_MUL_VP:
        result = node->lvp.v * *(node->rp);
        break;
    case WP_MUL_PP:
        result = *(node->lvp.p) * *(node->rp);
        break;
    case WP_DIV_VP:
        result = node->lvp.v / *(node->rp);
        break;
    case WP_DIV_PP:
        result = *(node->lvp.p) / *(node->rp);
        break;
    case WP_NEG_P:
        result = -*(node->lvp.p);
        break;
    default:
        yyerror("wp_ast_eval: unknown node type %d\n", node->type);
    }

    return result;
}

#endif

#endif
