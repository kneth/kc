/****************************************************************************
  Lexer and recursive-descent parser for kc's input language.

  A model has three sections, in this order:

    model      = { declaration } reaction { reaction } { initial }
    declaration= name "=" ( expr | '"' name '"' ) ";"
               | "#parameter" ( name | conc ) "=" expr { "," expr } ";"
               | "#print" ( name | conc ) { "," ( name | conc ) } ";"
    reaction   = number ":" side ( "->" | "<->" | "<=>" | "=" ) side ";" rates
               | name "'" "=" expr ";"
    side       = [ number ] species { "+" [ number ] species }
    rates      = { "c" "(" species ")" "=" expr ";" }
                 ( ( "k>" | "v>" ) "=" expr ";" [ ( "k<" | "v<" ) "=" expr ";" ]
                 | "K" "=" expr )
    initial    = conc [ "(0)" ] "=" expr ";"
               | name "(" species ")" "=" expr ";"
               | name "(0)" "=" expr ";"
    species    = name [ "(" ( "." | ( "+" | "-" ) [ number ] ) ")" ]
    conc       = "[" species "]"

  In expressions, unary minus binds tighter than "^", and all binary
  operators are left-associative.

  Errors are reported with their line number and counted. After a syntax
  error the parser skips to the end of the statement and carries on, so one
  run reports all errors.
****************************************************************************/

#include <ctype.h>
#include <float.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "misc.h"
#include "parser.h"
#include "symbmath.h"
#include "tableman.h"

static int errors;

/* Reports an error in the model at LINE. */
static void report_error(int line, const char *format, ...) {
  va_list args;

  fprintf(stderr, "Line %d: ", line);
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
  fputc('\n', stderr);
  errors++;
}

/*** Lexer ***/

typedef enum {
  T_EOF, T_NAME, T_NUMBER, T_FUNC, T_PARAM, T_PRINT, T_PRIME, T_LPAR, T_RPAR,
  T_TIME0, T_ONEWAY, T_TWOWAYS, T_EQUAL, T_RADICAL, T_LEFTARR, T_RIGHTARR,
  T_PLUS, T_MINUS, T_POW, T_MUL, T_DIV, T_LCONC, T_RCONC, T_COMMA, T_SEMI,
  T_COLON, T_QUOTE
} kind_t;

typedef struct {
  kind_t      kind;
  int         line;
  char        text[128];    /* the token as written */
  double      value;        /* T_NUMBER */
  expr_func_t func;         /* T_FUNC */
} token_t;

static const struct { const char *name; expr_func_t func; } functions[] = {
  {"exp", Exp}, {"log", Log}, {"ln", Ln}, {"sin", Sin}, {"cos", Cos},
  {"tan", Tan}, {"sinh", Sinh}, {"cosh", Cosh}, {"tanh", Tanh},
  {"asin", Asin}, {"acos", Acos}, {"atan", Atan}, {"asinh", Asinh},
  {"acosh", Acosh}, {"atanh", Atanh}
};

static const struct { const char *text; kind_t kind; } symbols[] = {
  {"#parameter", T_PARAM}, {"#print", T_PRINT}, {"<->", T_TWOWAYS},
  {"<=>", T_TWOWAYS}, {"(0)", T_TIME0}, {"->", T_ONEWAY}, {"**", T_POW},
  {"'", T_PRIME}, {"(", T_LPAR}, {")", T_RPAR}, {"=", T_EQUAL},
  {".", T_RADICAL}, {"+", T_PLUS}, {"-", T_MINUS}, {"^", T_POW},
  {"*", T_MUL}, {"/", T_DIV}, {"[", T_LCONC}, {"]", T_RCONC},
  {",", T_COMMA}, {";", T_SEMI}, {":", T_COLON}, {"\"", T_QUOTE}
};

static int lineno;
static int pushed[16], npushed;

static int getch(void) { return npushed > 0 ? pushed[--npushed] : getchar(); }
static void ungetch(int c) { pushed[npushed++] = c; }
static int  peekch(void) { int c = getch(); ungetch(c); return c; }

/* Consumes the characters of S if the input continues with them. */
static bool follows(const char *s) {
  size_t n = 0;
  int c;

  while (s[n] != '\0' && (c = getch()) == s[n])
    n++;
  if (s[n] == '\0')
    return true;
  ungetch(c);
  while (n > 0)
    ungetch(s[--n]);
  return false;
}

/* Skips a comment up to and including the closing '*' followed by CLOSE. */
static void skip_comment(int close) {
  int c, prev = 0, start = lineno;

  while ((c = getch()) != EOF) {
    if (c == '\n')
      lineno++;
    if (prev == '*' && c == close)
      return;
    prev = c;
  }
  report_error(start, "comment is not closed");
}

static void add_char(token_t *t, size_t *n, int c) {
  if (*n < sizeof t->text - 1)
    t->text[(*n)++] = (char)c;
  t->text[*n] = '\0';
}

/* number = digit { digit } | { digit } "." digit { digit } [ exponent ],
   where the exponent marker is e, E, l or L. */
static void scan_number(token_t *t, int c) {
  size_t n = 0;

  while (isdigit(c)) {
    add_char(t, &n, c);
    c = getch();
  }
  if (c == '.') {
    int d = getch();
    if (isdigit(d)) {
      add_char(t, &n, '.');
      for (c = d; isdigit(c); c = getch())
        add_char(t, &n, c);
      if (c == 'e' || c == 'E' || c == 'l' || c == 'L') {
        int marker = c, sign = getch(), digit = sign;
        if (sign == '-' || sign == '+')
          digit = getch();
        if (isdigit(digit)) {
          add_char(t, &n, marker);
          if (digit != sign)
            add_char(t, &n, sign);
          for (c = digit; isdigit(c); c = getch())
            add_char(t, &n, c);
        } else {
          ungetch(digit);
          if (digit != sign)
            ungetch(sign);
        }
      }
    } else {
      ungetch(d);
    }
  }
  ungetch(c);
  t->kind = T_NUMBER;
  char number[sizeof t->text];
  strcpy(number, t->text);
  for (char *p = number; *p != '\0'; p++)
    if (*p == 'l' || *p == 'L')
      *p = 'e';
  t->value = (double)strtof(number, NULL);
}

static void scan_name(token_t *t, int c) {
  size_t n = 0, len = 0;

  for (; isalnum(c); c = getch(), len++)
    add_char(t, &n, c);
  if (len >= STRING_LENGTH) {
    report_error(lineno, "name %s is longer than %d characters", t->text, STRING_LENGTH-1);
    t->text[STRING_LENGTH-1] = '\0';
  }
  t->kind = T_NAME;
  if ((strcmp(t->text, "k") == 0 || strcmp(t->text, "v") == 0) && (c == '>' || c == '<')) {
    add_char(t, &n, c);
    t->kind = c == '>' ? T_LEFTARR : T_RIGHTARR;
    return;
  }
  ungetch(c);
  for (size_t i = 0; i < sizeof functions / sizeof functions[0]; i++)
    if (strcmp(t->text, functions[i].name) == 0) {
      t->kind = T_FUNC;
      t->func = functions[i].func;
    }
}

/* Reads the next token into T. */
static void lex(token_t *t) {
  for (;;) {
    int c = getch();

    t->line = lineno;
    t->text[0] = '\0';
    if (c == EOF) {
      t->kind = T_EOF;
      return;
    }
    if (c == '\n') {
      lineno++;
      continue;
    }
    if (c == ' ' || c == '\t' || c == '\r')
      continue;
    if ((c == '/' || c == '(') && follows("*")) {
      skip_comment(c == '/' ? '/' : ')');
      continue;
    }
    if (isdigit(c) || (c == '.' && isdigit(peekch()))) {
      scan_number(t, c);
      return;
    }
    if (isalpha(c)) {
      scan_name(t, c);
      return;
    }
    ungetch(c);
    for (size_t i = 0; i < sizeof symbols / sizeof symbols[0]; i++)
      if (follows(symbols[i].text)) {
        strcpy(t->text, symbols[i].text);
        t->kind = symbols[i].kind;
        return;
      }
    report_error(lineno, "unexpected character '%c'", getch());
  }
}

/*** Parser ***/

typedef struct {
  char   name[STRING_LENGTH];
  double charge;            /* FLT_MAX for a radical */
} species_t;

static token_t tok, ahead;  /* the current token and, if has_ahead, the next */
static bool    has_ahead;
static int     last_line;   /* line of the token before tok */
static jmp_buf on_error;    /* where a syntax error continues */

static void next(void) {
  last_line = tok.line;
  if (has_ahead) {
    tok = ahead;
    has_ahead = false;
  } else
    lex(&tok);
}

static kind_t peek(void) {
  if (!has_ahead) {
    lex(&ahead);
    has_ahead = true;
  }
  return ahead.kind;
}

/* Reports at LINE that EXPECTED was expected at the current token, and
   abandons the statement. */
static _Noreturn void syntax_error_at(int line, const char *expected) {
  if (tok.kind == T_EOF)
    report_error(line, "expected %s, found end of input", expected);
  else
    report_error(line, "expected %s, found '%s'", expected, tok.text);
  longjmp(on_error, 1);
}

static _Noreturn void syntax_error(const char *expected) {
  syntax_error_at(tok.line, expected);
}

static bool accept(kind_t kind) {
  if (tok.kind != kind)
    return false;
  next();
  return true;
}

/* KIND must be a symbol. A missing symbol, typically a ';' at the end of a
   line, is reported on the line of the token before it. */
static void expect(kind_t kind) {
  char expected[16] = "";

  if (accept(kind))
    return;
  for (size_t i = 0; expected[0] == '\0'; i++)
    if (symbols[i].kind == kind)
      snprintf(expected, sizeof expected, "'%s'", symbols[i].text);
  syntax_error_at(last_line < tok.line ? last_line : tok.line, expected);
}

static void expect_name(char name[STRING_LENGTH]) {
  if (tok.kind != T_NAME)
    syntax_error("a name");
  strcpy(name, tok.text);
  next();
}

/* charge = [ "(" ( "." | ( "+" | "-" ) [ number ] ) ")" ] */
static double charge(void) {
  double q;

  if (!accept(T_LPAR))
    return 0.0;
  if (accept(T_RADICAL))
    q = FLT_MAX;
  else {
    if (accept(T_PLUS))
      q = 1.0;
    else if (accept(T_MINUS))
      q = -1.0;
    else
      syntax_error("'+', '-' or '.'");
    if (tok.kind == T_NUMBER) {
      q *= tok.value;
      next();
    }
  }
  expect(T_RPAR);
  return q;
}

static species_t species(void) {
  species_t s;

  expect_name(s.name);
  s.charge = charge();
  return s;
}

static species_t conc(void) {
  species_t s;

  expect(T_LCONC);
  s = species();
  expect(T_RCONC);
  return s;
}

static expr_t expr(void);

static expr_t primary(void) {
  expr_t e;

  switch (tok.kind) {
  case T_LPAR:
    next();
    e = expr();
    expect(T_RPAR);
    return e;
  case T_NUMBER:
    e = expr_const(tok.value);
    next();
    return e;
  case T_FUNC: {
    expr_func_t f = tok.func;
    next();
    expect(T_LPAR);
    e = expr();
    expect(T_RPAR);
    return expr_apply_func(f, e);
  }
  case T_NAME: {
    char name[STRING_LENGTH];
    double value;

    expect_name(name);
    value = GetConstant(name);
    if (GetError() != NotFound)
      return expr_const(value);
    if (strcmp(name, "time") == 0)
      NonAutoSystem();
    return expr_var(name);
  }
  case T_LCONC: {
    species_t s = conc();
    char name[STRING_LENGTH];

    if (accept(T_TIME0)) {
      double value = GetBeginConc(s.name, s.charge);
      if (GetError() == NoError)
        return expr_const(value);
      fprintf(stderr, "[%s(%e)] not found.\n", s.name, s.charge);
      return expr_const(0.0);
    }
    RenameSpec(name, s.name, s.charge);
    return expr_var(name);
  }
  default:
    syntax_error("an expression");
  }
}

static expr_t unary(void) {
  return accept(T_MINUS) ? expr_neg(unary()) : primary();
}

static expr_t power(void) {
  expr_t e = unary();

  while (accept(T_POW))
    e = expr_pow(e, unary());
  return e;
}

static expr_t term(void) {
  expr_t e = power();

  for (;;)
    if (accept(T_MUL))
      e = expr_mul(e, power());
    else if (accept(T_DIV))
      e = expr_div(e, power());
    else
      return e;
}

static expr_t expr(void) {
  expr_t e = term();

  for (;;)
    if (accept(T_PLUS))
      e = expr_add(e, term());
    else if (accept(T_MINUS))
      e = expr_sub(e, term());
    else
      return e;
}

/* Evaluates E, which must be a number; OK tells whether it was. */
static double eval(expr_t e, bool *ok) {
  double value = expr_eval(e);

  *ok = expr_get_error() != EXPR_NO_EVAL;
  return value;
}

/*** Declarations ***/

static void constant(void) {
  char name[STRING_LENGTH];
  int line;
  bool ok;
  double value;
  expr_t e;

  expect_name(name);
  expect(T_EQUAL);
  if (accept(T_QUOTE)) {
    char string[STRING_LENGTH];

    expect_name(string);
    expect(T_QUOTE);
    expect(T_SEMI);
    NewStrConst(name, string);
    return;
  }
  e = expr();
  line = tok.line;
  expect(T_SEMI);
  value = eval(e, &ok);
  if (ok)
    NewConstant(name, value);
  else
    fprintf(stderr, "Unable to evaluate expression in line %d.\n", line);
}

/* "#parameter" ( name | conc ) "=" init, low, high, delta, pref, direction ";" */
static void parameter(void) {
  enum { INIT, LOW, HIGH, DELTA, PREF, DIRECTION, FIELDS };
  bool is_conc = tok.kind == T_LCONC, ok[FIELDS];
  species_t s;
  double v[FIELDS];
  int line;

  if (is_conc)
    s = conc();
  else
    expect_name(s.name);
  expect(T_EQUAL);
  for (int i = 0; i < FIELDS; i++) {
    if (i > 0)
      expect(T_COMMA);
    line = tok.line;
    v[i] = eval(expr(), &ok[i]);
    if (!ok[i])
      fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", line);
  }
  expect(T_SEMI);
  if (is_conc) {
    if (ok[INIT])
      NewParamConc(s.name, s.charge, v[INIT]);
    NewLowHighPrefConc(s.name, s.charge, v[LOW], v[HIGH], v[PREF]);
    NewDirectForConc(s.name, s.charge, (int)v[DIRECTION]);
    NewDeltaConc(s.name, s.charge, v[DELTA]);
  } else {
    if (ok[INIT])
      NewParameter(s.name, v[INIT]);
    NewLowHighPrefParam(s.name, v[LOW], v[HIGH], v[PREF]);
    NewDirectForParam(s.name, (int)v[DIRECTION]);
    NewDeltaParam(s.name, v[DELTA]);
  }
}

static void print_list(void) {
  do {
    if (tok.kind == T_LCONC) {
      species_t s = conc();
      NewPrintConc(s.name, s.charge);
    } else {
      char name[STRING_LENGTH];
      expect_name(name);
      NewPrintVar(name);
    }
  } while (accept(T_COMMA));
  expect(T_SEMI);
}

/*** Reactions ***/

/* SIDE is 1 for the reactants and -1 for the products. */
static void reaction_side(int side) {
  do {
    double coeff = 1.0;
    species_t s;

    if (tok.kind == T_NUMBER) {
      coeff = tok.value;
      next();
    }
    s = species();
    SpecieInReaction(GetCurrentReaction(), s.name, s.charge);
    NewCoeff(GetCurrentReaction(), s.name, s.charge, coeff, side);
    NewSpecie(s.name, s.charge);
  } while (accept(T_PLUS));
}

/* ( "k" | "v" ) arrow "=" expr: DIRECT is -1 forwards and 1 backwards. */
static void rate(int direct) {
  char kind = tok.text[0];
  expr_t e;

  next();
  expect(T_EQUAL);
  e = expr();
  if (kind == 'k')
    NewRateConst(GetCurrentReaction(), direct, e);
  else
    NewRateExpr(GetCurrentReaction(), direct, e);
}

static void rates(void) {
  int r = GetCurrentReaction();
  char name[STRING_LENGTH];
  bool ok;

  while (tok.kind == T_NAME && peek() == T_LPAR) {
    species_t s;
    double value;
    expr_t e;
    int line;

    expect_name(name);
    expect(T_LPAR);
    s = species();
    expect(T_RPAR);
    expect(T_EQUAL);
    e = expr();
    line = tok.line;
    expect(T_SEMI);
    value = eval(e, &ok);
    if (strcmp(name, "c") != 0)
      fprintf(stderr, "c expected!\n");
    else if (!ok)
      fprintf(stderr, "Unable to evaluate expression in line %d.\n", line);
    else
      NewPowerConst(r, s.name, s.charge, value, -1);
  }
  if (tok.kind == T_LEFTARR) {
    rate(-1);
    expect(T_SEMI);
    if (tok.kind != T_RIGHTARR) {
      if (GetReactKind(r) == bi)
        fprintf(stderr, "The reaction in line %d is two-ways.\n", tok.line);
      return;
    }
    if (GetReactKind(r) != bi) {
      int line = tok.line;
      next();
      expect(T_EQUAL);
      (void) expr();
      expect(T_SEMI);
      fprintf(stderr, "The reaction in line %d is a >one-way< reaction or >equilibrium.<\n", line);
      return;
    }
    rate(1);
    expect(T_SEMI);
    return;
  }
  if (tok.kind != T_NAME)
    syntax_error("'k>', 'v>' or 'K'");
  expect_name(name);
  expect(T_EQUAL);
  if (strcmp(name, "K") != 0) {
    (void) expr();
    fprintf(stderr, "Syntax error: K expected in line %d.\n", tok.line);
  } else
    NewRateConst(r, 0, expr());
}

static void reaction(void) {
  Direc kind;

  NewReaction((int)tok.value);
  next();
  expect(T_COLON);
  reaction_side(1);
  if (accept(T_ONEWAY))
    kind = uni;
  else if (accept(T_TWOWAYS))
    kind = bi;
  else if (accept(T_EQUAL))
    kind = equi;
  else
    syntax_error("'->', '<->' or '='");
  AddReactionKind(GetCurrentReaction(), kind);
  reaction_side(-1);
  expect(T_SEMI);
  rates();
}

/* name "'" "=" expr ";" */
static void ode(void) {
  char name[STRING_LENGTH];
  expr_t e;

  expect_name(name);
  expect(T_PRIME);
  expect(T_EQUAL);
  e = expr();
  expect(T_SEMI);
  NewExpr(name, e);
  NewDynVar(name);
}

/*** Initial values ***/

static void initial(void) {
  char name[STRING_LENGTH];
  bool time0, ok;
  double value;
  expr_t e;

  if (tok.kind == T_LCONC) {
    species_t s = conc();

    time0 = accept(T_TIME0);
    expect(T_EQUAL);
    e = expr();
    expect(T_SEMI);
    value = eval(e, &ok);
    if (!time0)
      NewConstraint(s.name, s.charge, e);
    else if (!ok)
      fprintf(stderr, "Unable to evaluate expression.\n");
    else
      NewBeginConc(s.name, s.charge, value);
    return;
  }
  expect_name(name);
  if (accept(T_LPAR)) {
    species_t s = species();

    expect(T_RPAR);
    expect(T_EQUAL);
    e = expr();
    expect(T_SEMI);
    value = eval(e, &ok);
    NewSpecConst(s.name, s.charge, name, value);
    if (GetError() == NotFound && s.charge == 0.0)
      NewDynVarConst(s.name, name, value);
    return;
  }
  time0 = accept(T_TIME0);
  expect(T_EQUAL);
  e = expr();
  expect(T_SEMI);
  value = eval(e, &ok);
  if (!ok)
    fprintf(stderr, "Line %d: Unable to evaluate expression\n", tok.line);
  else if (!time0)
    fprintf(stderr, "Line %d: (0) expected\n", tok.line);
  else
    NewInitValue(name, value);
}

/*** Model ***/

/* Prints CHARGE as kc writes it after a species name. */
static void print_charge(double charge) {
  if (charge == FLT_MAX)
    printf(".");
  else if (charge > 0.0)
    printf("%d+", (int)charge);
  else if (charge < 0.0)
    printf("%d-", -(int)charge);
}

static void report(void) {
  char name[STRING_LENGTH], string[STRING_LENGTH];
  double charge;
  int kind;

  printf("Species used (*=dyn. var.):\n");
  for (int i = 1; i <= NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    printf(IsSpecInConstraint(name, charge) == 0 ? "  *%s(" : "   %s(", name);
    print_charge(charge);
    printf(")\n");
  }
  for (int i = 1; i <= NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    printf("  *%s\n", name);
  }
  if (NumOfParameter() != 0) {
    printf("Parameters declared:\n");
    for (int i = 1; i <= NumOfParameter(); i++) {
      GetParamNo(i, name, &charge, &kind);
      if (charge == 0.0 || kind == 1)
        printf("  %s\n", name);
      else {
        printf("  %s(", name);
        print_charge(charge);
        printf(")\n");
      }
    }
  }
  printf("Constants declared:\n");
  for (int i = 1; i <= NumOfConstants(); i++) {
    GetConstantNo(i, name);
    printf("  %s = %e\n", name, GetConstant(name));
  }
  for (int i = 1; i <= NumOfStrConst(); i++) {
    GetStrConstNo(i, name, string);
    printf("  %s = \"%s\"\n", name, string);
  }
}

typedef enum { DECLARATIONS, REACTIONS, INITIALS } section_t;

static section_t section;       /* the furthest section seen so far */
static bool      has_reactions;

/* Notes that a statement of section S, described by WHAT, starts here. */
static void enter(section_t s, const char *what) {
  if (s < section)
    report_error(tok.line, "%s must come before the %s", what,
                 section == REACTIONS ? "reactions" : "initial values");
  else
    section = s;
  has_reactions |= s == REACTIONS;
}

static void statement(void) {
  if (tok.kind == T_PARAM || tok.kind == T_PRINT || (tok.kind == T_NAME && peek() == T_EQUAL)) {
    enter(DECLARATIONS, "a declaration");
    if (accept(T_PARAM))
      parameter();
    else if (accept(T_PRINT))
      print_list();
    else
      constant();
  } else if (tok.kind == T_NUMBER) {
    enter(REACTIONS, "a reaction");
    reaction();
  } else if (tok.kind == T_NAME && peek() == T_PRIME) {
    enter(REACTIONS, "a differential equation");
    ode();
  } else if (tok.kind == T_NAME || tok.kind == T_LCONC) {
    enter(INITIALS, "an initial value");
    initial();
  } else
    syntax_error("a declaration, reaction or initial value");
}

/* Skips the rest of a statement after a syntax error: up to and including
   the next ';', or up to a reaction or '#' declaration on a later line.
   The rate constants of a broken reaction are skipped too. */
static void recover(void) {
  int line = last_line;

  while (tok.kind != T_EOF) {
    if (tok.line > line && (tok.kind == T_PARAM || tok.kind == T_PRINT ||
                            (tok.kind == T_NUMBER && peek() == T_COLON)))
      return;
    if (accept(T_SEMI)) {
      if (tok.kind != T_LEFTARR && tok.kind != T_RIGHTARR)
        return;
    } else
      next();
  }
}

int parse(void) {
  lineno = 1;
  errors = 0;
  section = DECLARATIONS;
  has_reactions = false;
  SetupTableMan();
  next();
  while (tok.kind != T_EOF)
    if (setjmp(on_error) == 0)
      statement();
    else
      recover();
  if (!has_reactions)
    report_error(tok.line, "the model has no reactions or differential equations");
  if (errors > 0)
    return errors;
  if (verbose == 1)
    report();
  if (IsNonAutoSystem() == 1) {
    NewDynVar("time");
    NewExpr("time", expr_const(1.0));
  }
  return 0;
}
