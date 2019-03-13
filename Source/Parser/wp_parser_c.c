#include "wp_parser_c.h"
#include "wp_parser.lex.h"
#include "wp_parser.tab.h"

struct wp_parser*
wp_c_parser_new (char const* body)
{
    YY_BUFFER_STATE buffer = yy_scan_string(body);
    yyparse();
    struct wp_parser* parser = wp_parser_new();
    yy_delete_buffer(buffer);
    return parser;
}
