/****************************************************************************
  Lexer and recursive-descent parser for kc's input language.
****************************************************************************/

#ifndef _PARSER_
#define _PARSER_

/* Reads a model from standard input into the tables of the table manager.
   Returns the number of errors found; the tables are complete only if it
   returns 0. */
extern int parse(void);

#endif /* _PARSER_ */
