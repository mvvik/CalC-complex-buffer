/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2019 Victor Matveev
 *
 *                            syntax.h
 *
 *   Common CalC definitions and main script parser/interpreter:
 *
 *              Script parser (class TokenString)
 *              Expression interpreter (class ExpressionObj)
 *              Utility error procedures
 *
 ************************************************************************
 
    This file is part of Calcium Calculator (CalC).

    CalC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CalC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CalC.  If not, see <https://www.gnu.org/licenses/>

 ************************************************************************/

#ifndef CALC_SYNTAX_H_included
#define CALC_SYNTAX_H_included

extern int  VERBOSE;
extern int  GEOMETRY;
extern int  GEOMETRY1;
extern int  GEOMETRY2;
extern int  GEOMETRY3;
extern int  DIMENSIONALITY;
extern char LABEL_DIM1[2];
extern char LABEL_DIM2[6];
extern char LABEL_DIM3[4];

//********************************************************************************************

#define ERR_NAN_DT_TINY_ERR 2

//********************************************************************************************

// Geometry flag: bits 1 & 2 - dimension; 
//                bits 3 & 4 - dimension 1 (x=0, r(1)=1, r(2)=2) 
//                bits 5 & 6 - dimension 2 (y=0, z=1, phi=2, theta=3) 
//                bits 7 & 8 - dimension 3 (z=0, phi=1) 

#define CARTESIAN1D   1   // 1 
#define CARTESIAN2D   2   // 2
#define CARTESIAN3D   3   // 3

#define SPHERICAL     9   // 1 + 2*4  (r2)
#define DISC          5   // 1 + 1*4  (r1)

#define CONICAL       58  // 2 + 2*4 + 3*16 (r2 theta)
#define POLAR         38  // 2 + 1*4 + 2*16 (r1 phi)
#define CYLINDRICAL   22  // 2 + 1*4 + 1*16 (r1 z)

#define SPHERICAL3D   123 // 3 + 2*4 + 3*16 + 1*64 (r2 theta phi) 
#define CYLINDRICAL3D 39  // 3 + 1*4 + 2*16 + 0*64 (r1 phi   z)

#define DIMENSION_MASK  3

//********************************************************************************************

#define MAX_STORAGE_LENGTH  1600000
#define MAX_TOKEN_NUM       500000

#define MAX_TOKEN_LENGTH  512
#define MAX_LINE_LENGTH   2048
#define MAX_STRING_LENGTH 2048

#define MAX_FORMULA_LENGTH 2048    // Maximal length of the formula (ExpressionObj) string

#define IF_NEST_LEVELS 10

#define MAX_SCRIPT_FILES 31 // 2 ^ SCRIPT_ID_BITS - 1
#define SCRIPT_ID_BITS   5    

#define LINE_CONTINUE_TOKEN "..."
#define COMMENT_TOKEN  "%"
#define CARR_RET_TOKEN ";"
#define ODE_TOKEN      "`"
#define VAR_TOKEN      ":="
#define ASSIGN_TOKEN   "="
#define REACTION_TOKEN "->"

#define LOOP_TOKEN  "for"
#define IF_TOKEN    "if"
#define THEN_TOKEN  "then"
#define ELSE_TOKEN  "else"
#define ENDIF_TOKEN "endif"

#define ALIAS_NUM 15     /* aliases */

const char aliasArray[2*ALIAS_NUM][18] =  // in each alias pair, the first element is the obsolete token 
  { "&&", "and",   "||", "or",   "Adaptive", "adaptive",   "input", "include",   "bc_type", "bc.define",
    "point.mute", "mute",  "peak", "max", "run", "Run", "ln(", "log(",  "Neumann", "Noflux", "4D", "binary",
    "plot.point.steps", "plot.steps.point",  "plot.1D.steps", "plot.steps.1D", "sphere", "volume", "sobstacle", "obstacle" };

#define SPECIAL_NUM 7   /* non-delimited tokens: (these have higher parsing priority than single-char tokens) */
const char specialArray[SPECIAL_NUM][6] = { ODE_TOKEN, VAR_TOKEN, ASSIGN_TOKEN, REACTION_TOKEN,
                                            LINE_CONTINUE_TOKEN, COMMENT_TOKEN, CARR_RET_TOKEN };  

//********************************************************************************************

#define MAX_CALL_LEVEL 500

#define NUMBER_TYPE     50
#define POINTER_TYPE    51

#define T_POWER 0
#define T_MULT  1
#define T_DIV   2
#define T_MOD   3
#define T_PLUS  4
#define T_MINUS 5
#define T_LE 6
#define T_LT 7
#define T_GE 8
#define T_GT 9
#define T_EQ 10
#define T_NEQ 11
#define T_AND 12
#define T_OR 13

#define TOP_BINARY 13

#define BR_CLOSE 14
#define BRACKET  15

#define T_SINH  16
#define T_COSH  17
#define T_TANH  18
#define T_SIN   19
#define T_COS   20
#define T_TAN   21
#define T_EXP   22
#define T_SQR   23
#define T_SQRT  24
#define T_THETA 25
#define T_SIGMA 26
#define T_LOG   27
#define T_LOG10 28
#define T_ABS   29
#define T_ATAN  30
#define T_INT   31
#define T_NOT   32
#define T_RAND  33

#define T_TOP_FN  33

#define T_UNARY_NOT   34
#define T_UNARY_PLUS  35
#define T_UNARY_MINUS 36

#define TYPE_NUM 37

#define priorityLevels 7  // highest priority below plus 1
const char priority[TOP_BINARY+1] = { 0, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 6 };

//*** make sure a token is not a prefix of any token further in the type_id array
//*** (otherwise make it "special" token)

const char type_id[TYPE_NUM][7] = { "^", "*", "/", "mod", "+", "-", "<=",
				    "<", ">=", ">", "==", "!=", "and", "or", ")", "(",
                    "sinh(", "cosh(", "tanh(", "sin(", "cos(", "tan(", 
                    "exp(", "sqr(", "sqrt(", "theta(", "sigma(", "log(", "log10(", "abs(", 
                    "atan(", "int(", "not(", "rand(", "!", "+", "-" };

//***************************************************************************

inline bool is_letter(const char c) {
  if ( (c >= 'a' && c <='z') || (c >= 'A' && c <='Z') || c == '_' ) return true;
  else return false; }

inline bool is_numeral(const char c) {
  return (c >= '0' && c <= '9') ?  true : false; }

inline bool is_sign(const char c) {
  if (c == '+' || c == '-') return true;  else return false; }

inline bool is_short_token(const char c) {
  if ( strchr("~!@#%^&*()[]=-+\\|;:/<>?,`",c) ) return true;
  else return false;
} 
 
inline bool is_delim(const char c) {
  if ( !is_letter(c) && !is_numeral(c) && !is_short_token(c) && !strchr(".{}\"\'$",c) ) return true;
  else return false;
} 

inline bool is_allowed_first(const char c) {
  if ( is_letter(c) || is_numeral(c) ) return true;
  else return false;
} 

char equal(const char *s1, const char *s2);
char equal_(const char *s1, const char *s2);  // returns "1" if s2 is a prefix of s1

inline bool isBinary(char tp) { if ( tp <= TOP_BINARY ) return true; else return false; }
inline bool isUnary(char tp)  { if ( tp >= T_UNARY_NOT && tp < TYPE_NUM ) return true; else return false; }
inline bool isOperand(char tp)    { return ((tp == NUMBER_TYPE) || (tp == POINTER_TYPE) || (tp == BR_CLOSE )) ? true : false; }
inline bool isFunction(char tp)   { return ( (tp >= BRACKET) && (tp <= T_TOP_FN) ) ?  true : false; }

inline bool isExpressionStart(char *s) {  // to signal the parser that the minus sign is a separate token
  if ( equal(s,ASSIGN_TOKEN) || equal(s,VAR_TOKEN) || equal(s,ODE_TOKEN) || equal(s,IF_TOKEN) || 
       equal(s, "[" ) || equal(s, "print" ) || equal(s,"append") || equal(s,"current") || equal(s,"tortuosity") ) 
  return true; else return false; }

inline bool isExpressionStop(char *s) {
  if (equal(s,CARR_RET_TOKEN) || equal(s,"to") || equal(s,"step") || equal(s, "then" )  || equal(s,"else") || 
      equal(s,"endif") || equal(s, "," ) ) 
  return true; else return false;}

inline bool isLineStart(char *s) {
  if (equal(s,CARR_RET_TOKEN) || equal(s,"then") || equal(s,"else") ) return true; else return false;}

inline bool isLineEnd(char *s) {
  if (equal(s,CARR_RET_TOKEN) || equal(s,"else") || equal(s,"endif") ) return true; else return false;}


//***************************************************************************

bool isFloat(const char *s);
bool isInt(const char *s);

inline char *StrCpy(const char *s) 
  { char *d = new char[strlen(s) + 1]; return strcpy(d, s); }


inline char *StrCpy(const char *s1, const char *s2) 
  { char *d = new char[strlen(s1) + strlen(s2) + 2]; strcpy(d,s1); strcat(d,"."); return strcat(d, s2); }


inline char *StrCpy(const char *s, int n) 
  { char *d = new char[n+1]; strncpy(d, s, n); d[n] = '\0'; return d; }

//*******************************************************************************


FILE *fopenAssure(const char *fname, const char *mode, const char *id="", const char *str = "");

char *makeMessage(const char *fmt, ...);

//**************************************************************************************************
                             
void globalError(char *message); //, const char *constMessage = "");

//***********************************************************************************************
//                         V A R I A B L E   L I S T   C L A S S
//***********************************************************************************************


class VarList          // An abstract class which contains a collection of varables with IDs
{                      // that can be retreaved using member ResolvePtr(), and whose addresses 
 public:               // can be retreaved using member ResolveID()
  VarList()  { };
  virtual ~VarList() { };
  virtual double *ResolveID(const char *, double **t=0) = 0;
  virtual const char *ResolvePtr(double *ptr) = 0;
};

//***********************************************************************************************
//                                 F I L E   P A R S E R
//***********************************************************************************************


class TokenString
{
  friend class LoopObj;
  friend class InterpolObj;

private:

 long length;
 char storage[MAX_STORAGE_LENGTH];
 char *token_ptr[MAX_TOKEN_NUM];
 long lineNum[MAX_TOKEN_NUM];
 char *scriptFileNames[MAX_SCRIPT_FILES];
 int  scriptNum;

public:

 long  token_num;

 // constructor and parsing routines:

 TokenString(const char *fname, const char *, const char *, int argc, char **argv);
 TokenString(const TokenString &);
 ~TokenString(); 

 void parse(const char *fname, const char *extra1, const char *extra2, int argc, char **argv);
 void getToken(char *&ptr0, char *&ptr1, char *, bool);
 void addToken(const char *token, int line=0, int script=0);
 void deleteToken(long, long);

 // routines returning the number of occurences of an argument token in the TokenString:

 int token_count(const char *token, long *p=0);
 int token_count(const char *t1, const char *t2, long *p=0)  // for tokens of form "t1.t2"
    { char *t = ::StrCpy(t1,t2); int n = token_count(t, p); delete [] t; return n; }
 int token2_count(const char *token1, const char *token2, long *p=0);
 int token3_count(const char *token1, const char *token2, const char *token3, long *p=0);

 // routines returning the location index of an argument token within the TokenString:

 long token_index(const char *token, int i = 1);
 long token_index(const char *t1, const char *t2, int i = 1)  // for tokens of form "t1.t2"
    { char *t = ::StrCpy(t1,t2); long n = token_index(t, i); delete [] t; return n; }
 long token2_index(const char *token1, const char *token2, int i = 1); 
 long token3_index(const char *token1, const char *token2, const char *token3, int i = 1); 

 // routines interpreting the numerical/string value of a token argument:

 bool isParam      (const char *s, long *);                // is followed by an assignment token
 bool isNumberParam(const char *s, long &, double *r = 0); // is left-hand side of an assignment evaluating to double
 bool isConst      (const char *s, double *r = 0);         // is constant numerical value or is GoodParam

 // routines interpreting the numerical/string value of a token specified by its location;
 // error message thrown on failure

 double get_double(long p);      
 int    get_int   (long p);         
 long   get_long  (long p);        
 char  *get_string(const char *);
 char  *get_string(const char *, char *);
 char  *get_string(long p, char *s);
 char  *get_string(long p);

 // Assert

 bool Assert(const char *var, const char *val);

 // routines interpreting an assigned value; preserve value of the argument if no assignment is found

 int    get_int_param   (const char *,   int  *   = 0,  const char *message = 0);
 long   get_long_param  (const char *,   long *   = 0,  const char *message = 0);
 char  *get_string_param(const char *,   char *   = 0);
 double get_param       (const char *t,  double * = 0,  const char *message = 0);
 double get2_param      (const char *t1, const char *t2, double *d=0, const char *message=0)  // for tokens of form "t1.t2"
    { char *t = ::StrCpy(t1,t2); double x = get_param(t, d, message); delete [] t; return x; }

 // utilities

 char  *line_string(long, char *, class VarList *VL = 0);
 void  printResults(class VarList *VL = 0);

 // routines parsing the arguments/variables following a given token on the same line:

 int   trail_pars(long ind, va_list p);
 int   trail_pars(long ind, ... );
 int   trail_pars(const char *token, int m, ... );
 int   trail_pars(const char *t1, const char *t2, int m, ... ); // for tokens of form "t1.t2"

 // auxiliary routines:

 int   tokens_to_eol(long p);  // tokens to end of statement, with #p included in the count
 long  lastInLine(long p)  { return  p - 1 + tokens_to_eol(p); } // position of last token in statement

 char equal(long p, const char *s) { return ::equal(token_ptr[p], s); }
 char equal(const char *s, long p) { return ::equal(s, token_ptr[p]); }
 // char equal_(long p, const char *s) { return ::equal_(token_ptr[p], s); }

 char *StrCpy(long p) { return ::StrCpy(token_ptr[p]); } 

 char *operator[](long index) const  { if (index < 0) return token_ptr[0];
                                      else if (index > token_num ) return token_ptr[token_num]; 
                                      else return token_ptr[index]; }
   
 void print(const char *fname);
 void print(FILE * = (FILE *)stderr);
 char *linePrint(long p, FILE *f = (FILE *)stderr);

 // exception handling routines

 void errorMessage(long p, char *message1, const char *message2="");
 char *getVarName( long p, const char *what, int offset = 0); // check if token number "p" can serve as a valid name tag
 void checkName( long p, const char *what); 

};

//***********************************************************************************************
//                       E X P R E S S I O N    I N T E R P R E T E R
//***********************************************************************************************

struct TermStruct
{
  char   type;
  double val;
  double *ptr;
};


double Evaluate(struct TermStruct TS);
int get_max_term_num(TokenString &, long);


//***********************************************************************************************

class ExpressionObj
{

  static int callLevel;

 public:

  struct TermStruct *term_array;
  int    term_num;
  char   formula[MAX_FORMULA_LENGTH+2];

  ExpressionObj(int n = 0)  { Allocate(n); }

  void Allocate(int n) {  term_array = new TermStruct[term_num = n]; strcpy(formula,""); }

  ExpressionObj(const ExpressionObj &EO) {  Allocate( EO.term_num );
                                            *this = EO; }
  

  ExpressionObj(TokenString &, long p0, const char * message = ">>> Syntax error in expression", 
                class VarList * = 0, long *pLast = 0,
                double * = 0, const char * = "", double * = 0, const char * = "", double * = 0, const char * = "");

  void operator=(const ExpressionObj &EO) {
    for (int i = 0; i < term_num; i++) term_array[i] = EO.term_array[i];
    strcpy(formula, EO.formula); 
  }

  ~ExpressionObj()    { delete [] term_array; }

  int    eliminate(int, int);
  void   pushOp(int op, int *opStack, int *indStack0, int *indStack1, int &stackHead, int &j);
  void   Optimize(int * = 0);
  double Evaluate(int * = 0);
  void   print(FILE *f = (FILE *)stderr);

  void  buildFormulaString(class VarList *, double *, const char *, double *, const char *, 
                                            double *, const char *, const char * = "");
};


//***********************************************************************************************



#endif
