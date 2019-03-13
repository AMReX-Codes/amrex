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

#endif
