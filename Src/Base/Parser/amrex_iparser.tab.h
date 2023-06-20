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

#ifndef YY_AMREX_IPARSER_AMREX_IPARSER_TAB_H_INCLUDED
# define YY_AMREX_IPARSER_AMREX_IPARSER_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef AMREX_IPARSERDEBUG
# if defined YYDEBUG
#if YYDEBUG
#   define AMREX_IPARSERDEBUG 1
#  else
#   define AMREX_IPARSERDEBUG 0
#  endif
# else /* ! defined YYDEBUG */
#  define AMREX_IPARSERDEBUG 0
# endif /* ! defined YYDEBUG */
#endif  /* ! defined AMREX_IPARSERDEBUG */
#if AMREX_IPARSERDEBUG
extern int amrex_iparserdebug;
#endif

/* Token kinds.  */
#ifndef AMREX_IPARSERTOKENTYPE
# define AMREX_IPARSERTOKENTYPE
  enum amrex_iparsertokentype
  {
    AMREX_IPARSEREMPTY = -2,
    AMREX_IPARSEREOF = 0,          /* "end of file"  */
    AMREX_IPARSERerror = 256,      /* error  */
    AMREX_IPARSERUNDEF = 257,      /* "invalid token"  */
    NODE = 258,                    /* NODE  */
    NUMBER = 259,                  /* NUMBER  */
    SYMBOL = 260,                  /* SYMBOL  */
    F1 = 261,                      /* F1  */
    F2 = 262,                      /* F2  */
    F3 = 263,                      /* F3  */
    EOL = 264,                     /* EOL  */
    FLRDIV = 265,                  /* "//"  */
    POW = 266,                     /* "**"  */
    GEQ = 267,                     /* ">="  */
    LEQ = 268,                     /* "<="  */
    EQ = 269,                      /* "=="  */
    NEQ = 270,                     /* "!="  */
    AND = 271,                     /* "and"  */
    OR = 272,                      /* "or"  */
    NEG = 273,                     /* NEG  */
    UPLUS = 274                    /* UPLUS  */
  };
  typedef enum amrex_iparsertokentype amrex_iparsertoken_kind_t;
#endif

/* Value type.  */
#if ! defined AMREX_IPARSERSTYPE && ! defined AMREX_IPARSERSTYPE_IS_DECLARED
union AMREX_IPARSERSTYPE
{

    struct amrex::iparser_node* n;
    int d;
    struct amrex::iparser_symbol* s;
    enum amrex::iparser_f1_t f1;
    enum amrex::iparser_f2_t f2;
    enum amrex::iparser_f3_t f3;


};
typedef union AMREX_IPARSERSTYPE AMREX_IPARSERSTYPE;
# define AMREX_IPARSERSTYPE_IS_TRIVIAL 1
# define AMREX_IPARSERSTYPE_IS_DECLARED 1
#endif


extern AMREX_IPARSERSTYPE amrex_iparserlval;


int amrex_iparserparse (void);


#endif /* !YY_AMREX_IPARSER_AMREX_IPARSER_TAB_H_INCLUDED  */
