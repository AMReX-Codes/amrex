/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_AMREX_PARSER_AMREX_PARSER_TAB_H_INCLUDED
# define YY_AMREX_PARSER_AMREX_PARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef AMREX_PARSERDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define AMREX_PARSERDEBUG 1
#  else
#   define AMREX_PARSERDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define AMREX_PARSERDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined AMREX_PARSERDEBUG */
#if AMREX_PARSERDEBUG
extern int amrex_parserdebug;
#endif

/* Token kinds.  */
#ifndef AMREX_PARSERTOKENTYPE
# define AMREX_PARSERTOKENTYPE
  enum amrex_parsertokentype
  {
    AMREX_PARSEREMPTY = -2,
    AMREX_PARSEREOF = 0,           /* "end of file"  */
    AMREX_PARSERerror = 256,       /* error  */
    AMREX_PARSERUNDEF = 257,       /* "invalid token"  */
    NODE = 258,                    /* NODE  */
    NUMBER = 259,                  /* NUMBER  */
    SYMBOL = 260,                  /* SYMBOL  */
    F1 = 261,                      /* F1  */
    F2 = 262,                      /* F2  */
    F3 = 263,                      /* F3  */
    EOL = 264,                     /* EOL  */
    POW = 265,                     /* "**"  */
    GEQ = 266,                     /* ">="  */
    LEQ = 267,                     /* "<="  */
    EQ = 268,                      /* "=="  */
    NEQ = 269,                     /* "!="  */
    AND = 270,                     /* "and"  */
    OR = 271,                      /* "or"  */
    NEG = 272,                     /* NEG  */
    UPLUS = 273                    /* UPLUS  */
  };
  typedef enum amrex_parsertokentype amrex_parsertoken_kind_t;
#endif

/* Value type.  */
#if ! defined AMREX_PARSERSTYPE && ! defined AMREX_PARSERSTYPE_IS_DECLARED
union AMREX_PARSERSTYPE
{

    struct amrex::parser_node* n;
    double d;
    struct amrex::parser_symbol* s;
    enum amrex::parser_f1_t f1;
    enum amrex::parser_f2_t f2;
    enum amrex::parser_f3_t f3;


};
typedef union AMREX_PARSERSTYPE AMREX_PARSERSTYPE;
# define AMREX_PARSERSTYPE_IS_TRIVIAL 1
# define AMREX_PARSERSTYPE_IS_DECLARED 1
#endif


extern AMREX_PARSERSTYPE amrex_parserlval;


int amrex_parserparse (void);


#endif /* !YY_AMREX_PARSER_AMREX_PARSER_TAB_H_INCLUDED  */
