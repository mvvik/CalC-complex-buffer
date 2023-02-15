/************************************************************************
*
*                   Calcium Calculator (CalC)
*            Copyright (C) 2001-2019 Victor Matveev
*
*                            syntax.cpp
*
*   Script parser (class TokenString)
*   Expression interpreter (class ExpressionObj)
*   Utility error procedures
*
****************************************************************************
 
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

#define _CRT_SECURE_NO_WARNINGS

#include "PlatformSpecific.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>  // for compatibility with Visual C++
#include <math.h>
#include "syntax.h"

extern int VERBOSE;

int ExpressionObj::callLevel = 0;

//**************************************************************************

char equal(const char *s1, const char *s2)  {
	if ( !s1 || !s2 ) return 0; 
	return ( (strcmp(s1, s2) == 0) && (strlen(s1) == strlen(s2)) ) ? 1 : 0; 
}

char equal_(const char *s1, const char *s2) { // returns "1" if s2 is a prefix of s1
	if ( !s1 || !s2 ) return 0; 
	return ( strncmp(s2, s1, strlen(s2)) == 0 ); 
}

//**************************************************************************

double Evaluate(struct TermStruct TS) { if (TS.type == POINTER_TYPE) return *TS.ptr;
else if (TS.type == NUMBER_TYPE) return TS.val;
else return 0.0; } 

//**************************************************************************

bool isFloat(const char *s)
{
	bool pointFlag = false, expFlag = false;
	size_t last = strlen(s) - 1;
	unsigned int  i;

	for (i=0; i < strlen(s); i++) {
		if ( is_numeral(s[i]) ) continue;
		else if ( s[i] == '.' ) { 
			if ( pointFlag || expFlag || (last == 0) ) return false; // a single dot character is not a number
			else pointFlag = true; 
		} 
		else if ( is_sign(s[i]) ) {
			if ( i == last ) return false;
			if ( (i == 0) && !( is_numeral(s[1]) || (s[1] == '.') ) ) return false;  // mantissa sign
			if ( (i > 0) && ( !( s[i-1] == 'e' || s[i-1] == 'E' ) || !is_numeral(s[i+1]) ) ) return false; // exponent sign
		}
		else if ( s[i] == 'e' || s[i] == 'E' ) {
			if ( (i == last) || (i == 0) ) return false;
			if ( !is_numeral(s[i-1]) || !(is_numeral(s[i+1]) || is_sign(s[i+1])) ) return false;
			if (expFlag) return false;
			expFlag = true;
		}
		else return false;
	}
	return true;
}

//***************************************************************************

bool isInt(const char *s)
{
	unsigned int  i;

	for (i=0; i < strlen(s); i++) {
		if ( is_numeral(s[i]) ) continue;
		else if ( is_sign(s[i]) ) { if (i > 0 ) return false; }
		else return false;
	}
	return true;
}

//*************************************************************************************

void TokenString::errorMessage(long p, char *message1, const char *message2) { 
	char message[1024];
	char *lineInfo = linePrint(p);

	snprintf(message, 1023, "\n\n*** Error on line %ld of script \"%s\":\n", 
		lineNum[p] >> SCRIPT_ID_BITS, scriptFileNames[ lineNum[p] & MAX_SCRIPT_FILES ]);
	if ( message1 ) { strcat(message, "    "); strcat(message, message1); strcat(message,"\n");
	delete [] message1; }
	if ( strlen(message2) )  { strcat(message, "    "); strcat(message, message2); strcat(message,"\n"); }
	strcat(message, lineInfo);
	strcat(message, "\n");

	delete [] lineInfo;

	fprintf(stderr, "%s\n", message);
	throw 1;
}

//**************************************************************************************************

void globalError(char *message) { 
	fprintf(stderr, "\n\n*** Error: %s\n", message);
	delete [] message;
	throw 1;
}

//*************************************************************************************

char *TokenString::getVarName(long p, const char *what, int offset) {

	char  *str = StrCpy(p); //get_string( p );
	if (offset) strcpy(str, token_ptr[p] + offset);

	try {  
		long pos;
		if ( isParam( str, &pos ) ) 
			errorMessage( pos, makeMessage("Bad %s name \"%s\": redefined here", what, str ) );

		checkName(p, what);

		if (p > 0) if ( this->equal(p - 1, "buffer") ) p--;
		if (p > 0) if ( this->equal(p - 1, "cooperative") ) p--;

		if (p > 0) 
			if ( !isLineStart( token_ptr[p - 1] ) )  {
				if ( !equal(p + 1, "min") && !equal(p + 1, "max") )
					errorMessage( p - 1, 0, "Expected a new line at this position" );
				else if (p > 1)
					if ( !isLineStart( token_ptr[p - 2]) )
						errorMessage( p - 2, 0, "Expected a new line at this position");
			}					
	}   
	catch (int ERR) {   perror("Error in getVarName: ");
	                    errorMessage(p, makeMessage("ERROR=%d: Bad %s definition", ERR, what) ); }

	return str;
}
//*************************************************************************************

void TokenString::checkName( long pos, const char *what) {

	if ( isFloat( token_ptr[pos] ) )
		errorMessage( pos, makeMessage("Bad %s name: cannot be a number", what ) );
	else if ( !is_letter(token_ptr[pos][0]) && !is_numeral(token_ptr[pos][0]) )
		errorMessage( pos, makeMessage("Bad %s name: must start with an alphanumeral", what ) );
}

//**************************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//***************************************************************************

TokenString::TokenString(const char *fname, const char *extra1, const char *extra2, int argc, char **argv)
{
	scriptNum = 0;
	length = 0;
	token_num = 0;
	token_ptr[0] = storage;
	lineNum[0] = 1;

	parse(fname, extra1, extra2, argc, argv);
	strcpy(token_ptr[token_num], CARR_RET_TOKEN);
	if (VERBOSE > 6) print();
}

//***************************************************************************

TokenString::~TokenString()  { 
	for (int i=0; i < scriptNum; i++) delete [] scriptFileNames[i];
}


//***************************************************************************

void TokenString::parse(const char *fname, const char *extra1, const char *extra2, int argc, char **argv) {

	FILE   *f;
	bool exprFlag = false;   // in expression, minus/plus signs are treated as separate tokens

	long   ifStart[IF_NEST_LEVELS], ifLevel = -1;
	bool  ifClause[IF_NEST_LEVELS], waitForEndif[IF_NEST_LEVELS], waitForElse[IF_NEST_LEVELS], reachedThen[IF_NEST_LEVELS];

	int currentLine = 0, currentScript = scriptNum;

	f = fopenAssure(fname, "r", "simulation script", "");
	scriptFileNames[scriptNum] = ::StrCpy(fname);
	if (++scriptNum > MAX_SCRIPT_FILES) 
		globalError( makeMessage("Cannot include script file \"%s\": the maximal number of allowed scripts is %d", 
		fname, MAX_SCRIPT_FILES) );

	//if (VERBOSE) fprintf(stderr, "========> Parsing script file \"%s\":\n", fname);

	char lineString[MAX_LINE_LENGTH], token[MAX_TOKEN_LENGTH];

	while (1)   {
		if ( token_num == 0 && extra1) { strcpy(lineString, extra1); strcat(lineString, extra2); }
		else { fgets(lineString,MAX_LINE_LENGTH,f); currentLine ++; if ( feof(f) ) break; }

		char cminus = -106;  // European keyboard long minus code
		char *ptr0 = lineString, *ptr1, *last = lineString + strlen(lineString), *pminus;
		long tnum_old = token_num;

		pminus = strchr(lineString, cminus);
		while (pminus) {
			*pminus = '-';
			pminus = strchr(lineString, cminus);
		}

		try {

			if ( strlen(lineString) >= MAX_LINE_LENGTH-1 ) 
				throw makeMessage("Exceeded the maximal line length of %d", MAX_LINE_LENGTH);

			do {  
				getToken(ptr0,ptr1,last, exprFlag);
				if ( (ptr1 - ptr0) >= MAX_TOKEN_LENGTH-1 ) 
					throw makeMessage("Exceeded the maximal token length of %d", MAX_TOKEN_LENGTH);

				strncpy(token,ptr0,ptr1-ptr0); token[ptr1-ptr0] = '\0';
				if ( ::equal(token, "debug") ) VERBOSE = 7;

				for (int i = 0; i < ALIAS_NUM; i++) 
					if (::equal(aliasArray[2*i], token) ) strcpy( token, aliasArray[2*i+1] );

				if ( ::equal(token, COMMENT_TOKEN) || ::equal(token, LINE_CONTINUE_TOKEN) )  break;

				if ( isExpressionStart(token) ) 
					exprFlag = true;  //***  An expression

				if ( ::equal(token, IF_TOKEN) )  {                //*** IF conditional statement
					if ( ++ifLevel == IF_NEST_LEVELS ) throw ::StrCpy("Exceeded \"if\" nesting limit");
					ifStart[ifLevel] = token_num + 1;
					waitForEndif[ifLevel] = waitForElse[ifLevel] = reachedThen[ifLevel] = false;
					if (ifLevel > 0 ) {
						if ( !reachedThen[ifLevel-1] ) throw ::StrCpy("Missing \"then\" after previous \"if\""); 
						if ( waitForEndif[ifLevel-1] || waitForElse[ifLevel-1] )  
							waitForEndif[ifLevel] = reachedThen[ifLevel] = true;  // skip this if
					}
				}

				if ( ifLevel >= 0 )                       //*** Skip token 
					if ( ( waitForEndif[ifLevel]  && !::equal(token, ENDIF_TOKEN) ) ||
						( waitForElse[ifLevel]   && !(::equal(token, ELSE_TOKEN) || ::equal(token, ENDIF_TOKEN) ) ) )  
					{ ptr0 = ptr1; continue; }

				if ( ::equal(token, THEN_TOKEN) ) { 
					if ( ifLevel < 0 ) throw ::StrCpy("\"then\" without an \"if\"");
					else if ( reachedThen[ifLevel] )  throw ::StrCpy("Missplaced \"then\" statement");
					addToken(token, currentLine, currentScript);
					exprFlag = false;
					ptr0 = ptr1;
					double condition = ExpressionObj(*this, ifStart[ifLevel], "Cannot evaluate the \"if\" condition").Evaluate();
					reachedThen[ifLevel] = true;
					ifClause[ifLevel] = ( condition > 0 ) ? true : false;
					if ( !ifClause[ifLevel] ) waitForElse[ifLevel] = true;
					continue;
				}

				else if ( ::equal(token, ELSE_TOKEN) ) {
					if ( ifLevel < 0 ) throw ::StrCpy("\"else\" without an \"if\" ");
					else if ( !reachedThen[ifLevel] ) throw ::StrCpy("\"else\" without a \"then\" "); 
					else if (ifClause[ifLevel]) waitForEndif[ifLevel] = true; else waitForElse[ifLevel] = false; 
				}

				else if ( ::equal(token, ENDIF_TOKEN) ) {
					if ( ifLevel < 0 ) throw ::StrCpy("\"endif\" without an \"if\" ");
					else if ( !reachedThen[ifLevel] ) throw ::StrCpy("\"endif\" without a \"then\" ");
					else ifLevel--;
					if ( (ptr0 = ptr1) >= last ) break;  //*** remove "endif" from token array: it might cause a syntax error 
					else continue;
				}

				else if ( ::equal(token,"include") ) {  //*** Import data file
					if ( (ptr0 = ptr1) >= last ) break;
					getToken(ptr0,ptr1,last, exprFlag);
					char *input_file = ::StrCpy(ptr0, ptr1-ptr0);
					char *filename = get_string(input_file);
					parse(filename, 0, 0, argc, argv);
					delete [] input_file; delete [] filename;
					if ( (ptr0 = ptr1) >= last ) break; 
					else continue;
				}
				else if ( ::equal(token,"echo") ) {    //*** Echo a message
					if ( (ptr0 = ptr1) >= last ) break;
					getToken(ptr0,ptr1,last, exprFlag);
					strncpy(token,ptr0,ptr1-ptr0); token[ptr1-ptr0] = '\0';
					char *s = get_string(token);
					fprintf(stderr,"\n >>> USER MESSAGE: %s\n",s);
					delete [] s;
					if ( (ptr0 = ptr1) >= last ) break; 
					else continue;
				}

				if (token[0] == '$') {
					if (token[1] == '$') snprintf(token, 511, "%d",argc);  //*** Get number of command-line arguments
					else    
					{                //*** Get a command-line argument
						double farg;
						if ( !isConst(token+1, &farg) ) throw makeMessage("Bad command-line argument \"%s\" ", token);
						int arg = int(farg + 0.5);
						if (arg >= argc) throw makeMessage("Command-line argument $%d not specified\n", arg);
						strcpy( token, argv[ arg ] );
					}
				}
				
				addToken(token, currentLine, currentScript);     //!!!  All contingencies checked -> add the token to the list  !!!
				ptr0 = ptr1;

			} while (ptr1 < last ); // loop until end of line

			if (ifLevel >= 0 && !reachedThen[ifLevel]) { // default "then" clause at line break
				addToken(THEN_TOKEN, currentLine, currentScript);
				ptr0 = ptr1;
				double condition = ExpressionObj(*this, ifStart[ifLevel], "Cannot evaluate the \"if\" condition").Evaluate();
				reachedThen[ifLevel] = true;
				ifClause[ifLevel] = ( condition > 0 ) ? true : false;
				if ( !ifClause[ifLevel] ) waitForElse[ifLevel] = true;
				exprFlag = false;
			}

			if ( !::equal(token, LINE_CONTINUE_TOKEN) ) {
				exprFlag = false;   // by default each new line is not an expression
				if  (token_num > tnum_old ) addToken(CARR_RET_TOKEN, currentLine, currentScript);
			}
		} catch(char *str) {
			fclose(f);
			globalError(makeMessage("%s on Line %d of file \"%s\":\n %s\n", str, currentLine, fname, lineString) ); }
	}  // while (1) - loop over lines of script

	fclose(f);
	//this->print();
}

//***************************************************************************

void TokenString::getToken(char *&ptr0, char *&ptr1, char *last, bool expression_flag)
{
	int i;
	ptr1 = ptr0;

	while ( is_delim( *ptr0 ) && ( ptr0 < last ) )  ptr1 = ++ptr0;
	if (ptr1 == last) return;

	if (*ptr0 == '\"' || *ptr0 == '\'' ) {
		char quote = *ptr0;
		ptr1 = ptr0; ++ ptr0;
		while ( *(++ptr1) != quote ) 
			if (ptr1 >= last) throw makeMessage("Unterminated string: %s", ptr0);
		*ptr1 = ' ';  // erase closing quote
		return;
	}

	for (i = 0; i < ALIAS_NUM; i++)      // *** 1. Look for alias tokens first
		if (::equal_(ptr0,aliasArray[2*i]) ) { ptr1 += strlen(aliasArray[2*i]); return; }

	if (expression_flag)                 // *** 2. Look for numerical expression tokens next
		for (i=0; i< TYPE_NUM; i++)         // except for "mod", "and" and "or" - no priority to these!
			if ( ::equal_(ptr0, type_id[i]) && i != 3 && i != 12 && i != 13) { ptr1 += strlen(type_id[i]); return; }

	for (i = 0; i < SPECIAL_NUM; i++)    // *** 3. Then look for other special (non-delimited) tokens
		if (::equal_(ptr0,specialArray[i]) ) { ptr1 += strlen(specialArray[i]); return; }

	if ( is_short_token(*ptr1) && !is_sign(*ptr1) ) { ptr1++; return; }   // 4. *** Single-character token

	if ( is_sign(*ptr1) ) {
		ptr1++;         // *** 5. A sign is a separate token only in an expression
	}

	if ( is_numeral(*ptr1) || (*ptr1 == '.' && is_numeral(ptr1[1]) ) ) {    // getting a number
		while ( is_numeral(*ptr1) ) ptr1++;
		if (*ptr1 == '.') ptr1++;
		while ( is_numeral(*ptr1) ) ptr1++;
		if ( (*ptr1 == 'e' || *ptr1 == 'E') && 
			( is_numeral(ptr1[1]) || ( is_sign(ptr1[1]) && is_numeral(ptr1[2]) ) ) ) {  // getting an exponent
				ptr1++;
				if ( is_sign(*ptr1) ) ptr1++;
				while ( is_numeral(*ptr1) ) ptr1++;
		}
		else                                                          // otherwise it's a weird identifier  
			while (!is_short_token(*ptr1) && !is_delim(*ptr1)) ptr1 ++; // (as in "plot 2D", for instance)
		return;
	}

	while (!is_short_token(*ptr1) && !is_delim(*ptr1)) ptr1 ++;  // some sort of an identifie
}

//***************************************************************************

void TokenString::addToken(const char *token, int line, int script)
{
	size_t tlength = strlen(token);
	if (tlength == 0) return;

	if (length + tlength + 2 >= MAX_STORAGE_LENGTH)
		throw("******** Scripts is too long: increase MAX_STORAGE_LENGTH in file \"source/syntax.h\"");

	if (token_num + 2 >= MAX_TOKEN_NUM)
		throw("******** Script is too long: increase MAX_TOKEN_NUM in file \"source/syntax.h\"");

	lineNum[token_num] = (line << SCRIPT_ID_BITS) + script;
	strncpy( token_ptr[token_num], token, tlength + 1 );
	*(token_ptr[token_num] + tlength) = 0;

	if (VERBOSE > 6) fprintf(stderr, "%s ", token_ptr[token_num]);

	token_ptr[token_num+1] = token_ptr[token_num] + tlength + 1;
	token_num++;
	length += tlength + 1;
}

//***************************************************************************

void TokenString::deleteToken(long p0, long p1)
{

	if ( p1 >= token_num || p0 < 0 )
		throw makeMessage("Cannot delete tokens %ld through %ld (token number = %ld)", p0, p1, token_num);

	long dp = p1 - p0 + 1 ;
	token_num -= dp;

	for ( long p = p0 ; p < token_num; p++ ) { 
		token_ptr[p] = token_ptr[p + dp]; 
		lineNum[p] = lineNum[p + dp]; 
	}

}
//***************************************************************************

TokenString::TokenString(const TokenString &tstring)
{
	token_num = tstring.token_num;
	length    = tstring.length;
}

//***************************************************************************

int TokenString::token_count(const char *token, long *p)
{
	int cnt = 0;
	for (long i=0; i<token_num; i++)  
		if ( equal(i, token) ) { if (p && !cnt) *p = i; cnt++; }

	return cnt;
}

//***************************************************************************

int TokenString::token2_count(const char *token1, const char *token2, long *p)
{
	int cnt = 0;
	for (long i = 1; i < token_num; i++)
		if ( equal(i-1, token1) && equal(i, token2) ) { if (p && !cnt) *p = i; cnt++; }
	return cnt;
}

//***************************************************************************

int TokenString::token3_count(const char *token1, const char *token2, const char *token3, long *p)
{
	int cnt = 0;
	for (long i = 2; i < token_num; i++)
		if ( equal(i-2, token1) && equal(i-1, token2) && equal(i, token3) ) { if (p&& !cnt) *p = i; cnt++; }
	return cnt;
}

//***************************************************************************

void TokenString::print(const char *fname)
{
	FILE *f = fopenAssure(fname, "w", "TokenString printing","");
	print(f);
}

//***************************************************************************

void TokenString::print(FILE *f)
{

	fprintf(f, "\n### Parameter file parsed: { ");
	for (long i = 0; i < token_num; i++)
		fprintf(f, "%s ", token_ptr[i]);
	fprintf(f, " }\n");
}

//***************************************************************************

char *TokenString::linePrint(long p, FILE *f)
{
	char message[1024];  strcpy(message, "");
	long p0 = p - 1, p1 = p + 1;
	long i, j;

	while (p0 >= 0) 
		if ( lineNum[p0] == lineNum[p] ) p0--; else break;
	p0 ++;

	while (p1 < token_num) 
		if ( !isLineEnd( token_ptr[ p1 ] ) )   p1++; else break;
	// if (lineNum[p1] == lineNum[p] ) p1++; else break;
	//p1 -= 2; // ignore the carriage return 
	p1 --; 

	for (i = p0; i <= p1; i++) {
		strcat(message," ");
		strcat(message,token_ptr[i]);
	}

	strcat(message,"\n ");

	for (i = p0; i <= lastInLine(p); i++) 
		for (j = 0; j < long(strlen(token_ptr[i]) + 1); j++) 
			if (i < p) strcat(message," ");
			else strcat(message,"^");

	return ::StrCpy(message);
}
//***************************************************************************

int TokenString::trail_pars(long index, ... )
{
	va_list p;
	va_start(p, index);

	return trail_pars(index, p);
}

//***************************************************************************

int TokenString::trail_pars(const char *token, int m, ... )
{
	long ind = token_index(token, m);

	va_list p;
	va_start(p, m);

	return trail_pars(ind, p);

}

//***************************************************************************

int TokenString::trail_pars(const char *token1, const char *token2, int m, ... )
{

	long ind = token_index(token1, token2, m);

	va_list p;
	va_start(p, m);

	return trail_pars(ind, p);

}
//***************************************************************************

int TokenString::trail_pars(long ind0, va_list p)
{
	char    type;
	int    *ipar;
	double *dpar;
	long   *lpar;
	char   *cpar;
	char   *spar;
	long    ind = ind0;

	while (1)
	{
		ind ++; 
		if ( isLineEnd( token_ptr[ ind ] ) ) break;
		type = char( va_arg(p, int) );
		if (type == 'E' ) break;

		switch (type)
		{
		case 'd': dpar =  va_arg(p, double *);
			if (dpar) *dpar = get_double(ind); 
			break;
		case 'i': ipar =  va_arg(p, int *);   
			if (ipar) *ipar = get_int(ind); 
			break;
		case 'l': lpar =  va_arg(p, long *);
			if (lpar) *lpar = get_long(ind); 
			break;
		case 'c': cpar =  va_arg(p, char *); 
		    if (cpar) *cpar = *token_ptr[ind];
			break;
		case 's': spar =  va_arg(p, char *); 
		    if (spar) get_string(ind,spar);
			break;
		case 'S': spar =  va_arg(p, char *); 
		    if (spar) strcpy(spar,token_ptr[ind]);
			break;
		default:  va_arg(p, int); 
		          break;  // ignore argument 
		}
	}
	va_end(p);

	return  ind - ind0 - 1; 
}

//***************************************************************************

long TokenString::token_index(const char *token, int n) 
{ 
	int cnt = 0;
	for (long i = 0; i < token_num; i++)
	{
		if ( equal(i, token) ) cnt ++;
		if (cnt == n) return i;
	}
	throw makeMessage("Less than %d instances of \"%s\"",n,token);
}

//***************************************************************************

long TokenString::token2_index(const char *token1, const char *token2, int n)
{ 
	int cnt = 0;
	for (long i = 1; i < token_num; i++)
	{
		if ( equal(i-1, token1) && equal(i, token2) ) cnt ++;
		if (cnt == n) return i;
	}
	throw makeMessage("Less than %d instances of \"%s %s\"",n,token1,token2);
}


//***************************************************************************

long TokenString::token3_index(const char *token1, const char *token2, const char *token3, int n)
{ 
	int cnt = 0;
	for (long i = 2; i < token_num; i++)
	{
		if ( equal(i-2, token1) && equal(i-1, token2) && equal(i, token3) ) cnt ++;
		if (cnt == n) return i;
	}
	throw makeMessage("Less than %d instances of \"%s %s %s\"",n,token1,token2,token3);
}

//***************************************************************************

bool TokenString::isConst(const char *s, double *result)
{ 
	long pos;

	if (::isFloat(s)) { 
		if (result) sscanf(s, "%lf", result);
		return true; }
	else return isNumberParam(s, pos, result);
}

//***************************************************************************

bool TokenString::isParam(const char *s, long *pFirst) { 

	if ( !( is_letter(s[0]) || is_numeral(s[0]) ) ) return false;

	int n = token2_count(s, ASSIGN_TOKEN, pFirst);

	if (n == 0) return false;

	if (n > 1 && VERBOSE > 2 )
		if ( n - token3_count("for",s,ASSIGN_TOKEN) > 1 )
			fprintf(stderr,"\n*** WARNING: more than one definition for \"%s\" \n",s);

	checkName(*pFirst - 1, "parameter");

	// have to check for loop token so that loop parameter would be evaluatable, in case
	// it is encountered in preprocessor statements like the "if" statement

	if (*pFirst > 1) 
		if ( !isLineStart( token_ptr[ *pFirst - 2 ] )  && !equal(*pFirst - 2, LOOP_TOKEN) )  
			errorMessage( *pFirst-2, 0, "Expected a new line at this position" );

	(*pFirst)++;
	return true;
}

//***************************************************************************

bool TokenString::isNumberParam(const char *s, long &pFirst, double *result)
{ 
	long  pLast;

	if ( strchr( s, '{') ) {   // get an array element
		char arg[1024];
		strcpy(arg, s);
		char *name = strtok( arg, "{");
		if ( !isParam(name, &pFirst) ) {
			if ( token_count(name, &pFirst) ) pFirst++; // access definition keywords (like "grid n m")
			else return false;
		}

		char  *sind;
		if ( !( sind = strtok( NULL, "}") ) ) errorMessage( token_index(s), 0, "No closing curly bracket" );  

		double dind;
		if ( !isConst(sind, &dind) )  errorMessage( token_index(s), 0, "Bad array index");  

		long signed ind   = pFirst + long( dind + 0.5 ) - 1;
		pLast = lastInLine(pFirst);
		if ( ind > pLast) 
			errorMessage( token_index(s), makeMessage("Array index is too large: %d > %d", 
			int( ind-pFirst ), int(pLast-pFirst)) );
		else if ( ind < pFirst - 1 ) 
			errorMessage( token_index(s), 0, "Array index cannot be negative");

		if (result) *result = ( ind == pFirst - 1 ? double(pLast-pFirst+1) : get_double(ind) );
		return true;
	}

	if ( !isParam(s, &pFirst) )  return false;

	double x = ExpressionObj( *this, pFirst, "", 0, &pLast).Evaluate();
	if (result) *result = x;

	if ( pLast >= pFirst && isExpressionStop( token_ptr[pLast+1] ) ) return true;
	else return false;
}

//***************************************************************************

double TokenString::get_param(const char *s, double *res, const char *message)
{
	double r;
	long pos;

	if ( !isNumberParam(s, pos, &r) )  {
		if (message) errorMessage( token_index(s), 0, message );
		else return -1.0;
	}

	if (res) *res = r;
	return r;
}

//***************************************************************************

int TokenString::get_int_param(const char *s, int *ires, const char *message) {

	double r, res;
	int ir;

	if (ires) res = double( *ires );
	r = get_param(s, &res, message);
	if (ires) r = res;
	ir = int(r);

	if ( (r - ir) < 1e-12 ) { if (ires) *ires = ir; return ir; }
	else  errorMessage( token_index(s), 0, "Expected an integer value");
	return 1;
}

//***************************************************************************

long TokenString::get_long_param(const char *s, long *ires, const char *message) {

	double r, res;
	long ir;

	if (ires) res = double( *ires );
	r = get_param(s, &res, message);
	if (ires) r = res;
	ir = long(r);

	if ( (r - ir) < 1e-12 ) { if (ires) *ires = ir; return ir; }
	else errorMessage( token_index(s), 0, "Expected an integer value");
	return 1;
}

//***************************************************************************

char *TokenString::get_string_param(const char *str, char *result) {

	long  pos;

	if ( isParam( str, &pos ) ) { 
		if (result) return get_string( pos, result);
		else return get_string(pos);
	}
	else if (result) return result;
	else return 0;
}

//***************************************************************************

bool TokenString::Assert(const char *var, const char *val) {

	if (token2_count(var,val) || token3_count(var,ASSIGN_TOKEN,val)) return true;  
	else return false;
}

//***************************************************************************

double TokenString::get_double(long p)
{
	double r;

	if (::isFloat(token_ptr[p]) ) { sscanf(token_ptr[p],"%lf",&r); return r; }
	checkName(p, "float parameter");
	return get_param(token_ptr[p], 0, "Expected a float value");
}

//***************************************************************************

int TokenString::get_int(long p) {
	int ir;
	if (::isInt(token_ptr[p]) ) { sscanf(token_ptr[p],"%d",&ir); return ir; }
	checkName(p, "integer parameter");
	return get_int_param(token_ptr[p], 0, "Expected an integer");
}

//***************************************************************************

long TokenString::get_long(long p) {
	long   lr;
	if (::isInt(token_ptr[p]) ) { sscanf(token_ptr[p],"%ld",&lr); return lr; }
	checkName(p, "long integer parameter");
	return get_long_param(token_ptr[p], 0, "Expected a long integer");
}

//***************************************************************************

char  *TokenString::get_string(long p, char *s) { 
	if ( equal(p, CARR_RET_TOKEN) ) errorMessage(p, 0, "String value expected");
	//   if ( is_delim(token_ptr[p][0] ) ) errorMessage(p, 0, "String value expected");
	else return get_string(token_ptr[p], s); 
	return 0;
}

//***************************************************************************

char  *TokenString::get_string(long p)          {  
	if ( equal(p, CARR_RET_TOKEN) ) errorMessage(p, 0, "String value expected");
	//   if ( is_delim(token_ptr[p][0] ) ) errorMessage(p, 0, "String value expected");
	else  return get_string(token_ptr[p]); 
	return 0;
}

//***************************************************************************

char *TokenString::get_string(const char *sin, char *sout)
{
	long p;
	if ( ::equal(sin, CARR_RET_TOKEN) ) throw ::StrCpy("String value expected");
	//   if ( is_delim( sin[0] ) ) throw ::StrCpy("String value expected");

	if ( !isParam(sin, &p) )  strcpy(sout,sin);  // simple string
	else  line_string(p, sout);

	return sout;
}

//***************************************************************************

char *TokenString::get_string(const char *sin)
{
	char sout[MAX_STRING_LENGTH];
	get_string(sin, sout);
	return ::StrCpy(sout);
}

//***************************************************************************

char *TokenString::line_string(long p, char *sout, VarList *VL )
{
	long pp, p1 = lastInLine(p);

	strcpy(sout,"");

	while (p <= p1) {

		double x = ExpressionObj(*this, p, "", VL, &pp).Evaluate(); 

		if (pp >= p) {
			char temp[40];
			p = pp + 1;
			snprintf(temp, 39, "%.12g", x);
			strcat(sout,temp);
		}
		else {
			char *temp = get_string( p );
			strcat(sout, temp);
			p++;
			delete [] temp;
		}
	}

	return sout;
}
//**********************************************************************************************


void TokenString::printResults(class VarList *VL) {

	int    i;
	long   p;
	char   *fileStr = NULL;
	FILE   *f;
	char   temp[4000];
	bool   fileSpec = false;

	if (token_count("print.file", &p)) {
		p++;
		fileStr = ::StrCpy( line_string( p + ( equal(p, ASSIGN_TOKEN) ? 1 : 0), temp, VL ) );
		fileSpec = true;
	}

	for (i = 0; i < token_count("print"); i++) {
		p = token_index("print", i+1) + 1;
		if ( !fileSpec ) fileStr = get_string(p++); 
		f = fopenAssure(fileStr, "w", "the output of \"print\" statements", "");
		fprintf( f, "%s\n",  line_string( p, temp, VL ) );
		if ( (f != (FILE *)stdout) && (f != (FILE *)stderr) ) fclose(f); 
	}

	for (i = 0; i < token_count("append"); i++) {
		p = token_index("append", i+1) + 1;
		if  ( !fileSpec ) fileStr = get_string(p++);
		f = fopenAssure(fileStr, "a", "the output of \"append\" statements", "");
		fprintf( f, "%s\n", line_string( p, temp, VL ) );
		if ( (f != (FILE *)stdout) && (f != (FILE *)stderr) ) fclose(f); 
	}
}

//***************************************************************************

int TokenString::tokens_to_eol(long p)
{
	int    num = 0;
	while ( !isLineEnd(token_ptr[ p ]) && p < token_num ) { p++; num++; }
	return num;
}

//**************************************************************************************************
//                        E X P R E S S I O N    I N T E R P R E T E R
//**************************************************************************************************

inline double binary(char op, double &a, double b)
{
	switch(op) {
	case T_PLUS:  return a += b;
	case T_MINUS: return a -= b;
	case T_OR:    return a = ( (a > 0) || (b > 0) ? 1: 0 );
	case T_AND:   return a = ( (a > 0) && (b > 0) ? 1: 0 );
	case T_GT:    return a = (a > b  ? 1 : 0);
	case T_GE:    return a = (a >= b ? 1 : 0);
	case T_LT:    return a = (a < b  ? 1 : 0);
	case T_LE:    return a = (a <= b ? 1 : 0);
	case T_EQ:    return a = (a == b ? 1 : 0);
	case T_NEQ:   return a = (a != b ? 1 : 0);
	case T_MULT:  return a == 0.0 || b == 0.0 ? (a = 0.0) : a *= b;
	case T_DIV:   return a == 0.0 && b == 0.0 ? 1.0 : ( b != 0 ? a /= b : 0.0 );
	case T_MOD:   return a = double ( int(a) % int(b) );
	case T_POWER: return a = pow(a,b);
	default:      throw makeMessage("Unknown binary operator %d", op);
		return 1;
	}
}

//**************************************************************************

inline double function(char op, double x)
{
	switch (op)
	{
	case BRACKET:   return x;
	case T_RAND:    return rand();
	case T_NOT:     return (x > 0) ? 0 : 1;
	case T_INT:     return double(int(x));
	case T_COSH:    return cosh(x);
	case T_SINH:    return sinh(x);
	case T_COS:     return cos(x);
	case T_SIN:     return sin(x);
	case T_TANH:    return tanh(x);
	case T_TAN:     return tan(x);
	case T_ATAN:    return atan(x);
	case T_EXP:     return exp(x);
	case T_LOG:     return log(x);
	case T_LOG10:   return log10(x);
	case T_SQR:     return x * x;
	case T_ABS:     return fabs( x );
	case T_SQRT:    return sqrt(x);
	case T_THETA:   if (x == 0.0 ) return 0.5; else return (x > 0) ? 1 : 0;
	case T_SIGMA:   return 0.5 * (1 + tanh(x) );
	default:        throw makeMessage("Unknown function %d", op);
		return 1;
	}
}

//**************************************************************************

double ExpressionObj::Evaluate(int *firstTerm)
{
	struct TermStruct term;
	int    pr, unaryOp = 0;
	int    opStack[priorityLevels];
	double val, valStack[priorityLevels + 1];
	int    j, j0, stackHead = 0;

	if (term_num <= 0 ) return 0.0;
	if (firstTerm) j0 = *firstTerm;
	else j0 = 0;

	for (j = j0; j < term_num; j ++) { // *** LOOP OVER TERMS BEGINS HERE ***

		int Type = (term = term_array[j]).type;

		if ( isUnary(Type) ) {   // Catch unary sign operators here
			unaryOp = term.type; continue;
		}

		if ( isBinary(Type) ) {  // ...So this binary operator is not a unary operator
			pr = priority[ Type ];
			while ( stackHead )
				if ( priority[ opStack[stackHead - 1] ] <= pr ) {
					stackHead --; 
					binary( opStack[stackHead], valStack[stackHead], valStack[stackHead+1] );
				} else break;
			opStack[ stackHead++ ] = Type;
			//pushOp( Type, valStack, valHead, opStack, stackHead ); 
			continue;
		}

		if ( Type == BR_CLOSE )  break;   // catch a closing bracket
		// *** The term is a number:
		if ( isFunction(Type) ) {                          // 1. a function call (or nested expression)
			j ++;
			val = function( Type, this->Evaluate(&j) );
		}
		else if ( Type == NUMBER_TYPE )  val = term.val;   // 2. a float      
		else if ( Type == POINTER_TYPE ) val = *term.ptr;  // 3. a variable
		else throw makeMessage("Cannot evaluate expression { %s }: unknown term %d\n", formula, Type); 

		if (unaryOp) {
			if  (unaryOp == T_UNARY_NOT)  val = (val > 0) ? 0 : 1;
			unaryOp = 0;
		}

		valStack[ stackHead ] = val;
	}

	if (firstTerm) *firstTerm = j;
	while ( stackHead-- ) 
		binary( opStack[stackHead], valStack[stackHead], valStack[stackHead+1] ); // *** Collect the result

	if (!_finite(valStack[0]))
		throw makeMessage("\n Not-a-number returned by the following expression: \n\n { %s }\n", formula);
	else return valStack[0] ;
}

//**************************************************************************

int ExpressionObj::eliminate(int kill0, int kill1) {

	if (kill0 < 0 || kill1 >= term_num)   
		throw makeMessage("Cannot kill terms #%d through #%d (max term number = %d)\n", kill0, kill1, term_num);

	int d = kill1 - kill0 + 1;
	for (int i = kill1 + 1; i < term_num; i++) term_array[i-d] = term_array[i];
	term_num -= d;
	return d;
}

// *************************************************************************

void ExpressionObj::pushOp(int op, int *opStack, int *indStack0, int *indStack1, int &stackHead, int &j) {

	int pr = ( op == BR_CLOSE ) ? priorityLevels : priority[ op ];
	int reduce, old, pr_reduce, pr_old;

	while ( stackHead )
		if ( (pr_reduce = priority[ reduce = opStack[stackHead - 1] ]) <= pr ) {

			old = TYPE_NUM; pr_old = priorityLevels;
			if ( stackHead > 1 ) pr_old = priority[old = opStack[stackHead - 2]]; 

			if ( (old == T_DIV)   && (reduce == T_DIV)   )  reduce = T_MULT;  // ** A / B / C = A / (B * C)
			else if ( (old == T_DIV)   && (reduce == T_MULT)  )  reduce = T_DIV;   // ** A / B * C = A / (B / C)
			if ( (old == T_MINUS) && (reduce == T_PLUS)  )  reduce = T_MINUS; // ** A - B + C = A - (B - C)
			else if ( (old == T_MINUS) && (reduce == T_MINUS) )  reduce = T_PLUS;  // ** A - B - C = A - (B + C)

			int j0 = j;
			int shm1 = stackHead-1;
			int j1  = indStack1[shm1];       // Points to the _end_ of the 1st operand
			int j10 = indStack0[shm1];
			int j2  = indStack1[stackHead];  // Points to the _end_ of the 2nd operand

			struct TermStruct t1 = term_array[j1];
			struct TermStruct t2 = term_array[j2];

			// ********************************************************

			if ( t1.type == NUMBER_TYPE && t2.type == NUMBER_TYPE && pr_old >= pr_reduce ) { // *** Binary with two constants
				binary( reduce, t1.val, t2.val );
				term_array[j1] = t1;
				j -= eliminate(j1+1, j2);
			}
			else if ( t2.type == NUMBER_TYPE ) {                   // *** Second operand is a number
				if (t2.val == 0.0 ) {                                // *** Check 0
					if  ( reduce == T_MULT || reduce == T_AND ) {      // {expr} * 0 = 0
						while (stackHead > 1) {
							if (priority[ opStack[stackHead - 2] ] > pr_reduce)  break; 
							stackHead--; shm1--; j10 = indStack0[shm1]; }
						j -= eliminate(j10, j2-1);
						indStack1[shm1] = j10; }
					else if ( reduce == T_PLUS || reduce == T_MINUS || reduce == T_OR ) 
						j -= eliminate(j1+1, j2);                        // ptr (+,-,or) 0 = ptr
					else if ( reduce == T_POWER && old != T_POWER ) {  // ptr^0 = 1
						term_array[j2].val = 1.0;
						j -= eliminate(j10, j2-1);
						indStack1[shm1] = j10;
					}
					else if ( reduce == T_DIV || reduce == T_MOD  )    // ptr (/,%) 0 = error
						throw ::StrCpy("Could not simplify the expression: division by zero");
				}
				else if (t2.val == 1.0) {                            // *** Check 1.0
					if (reduce == T_MULT || reduce == T_DIV || reduce == T_MOD || reduce == T_AND || reduce == T_POWER ) 
						j -= eliminate(j1+1, j2);                        // ptr (*,/,mod,and,^) 1 = ptr
					else if (reduce == T_OR) {                         // {expr} or 1 = 1
						while (stackHead > 1)  {
							if  (priority[ opStack[stackHead - 2] ] > pr_reduce) break; 
							stackHead--; shm1--; j10 = indStack0[shm1]; }
						j -= eliminate(j10, j2-1);
						indStack1[shm1] = j10;
					}
				}
			} //endif second operand is a number
			else if ( t1.type == NUMBER_TYPE ) {                   // *** First operand is a number
				if (t1.val == 0.0 ) {                                // *** Check 0.0
					if (reduce == T_MULT || reduce == T_DIV  || reduce == T_AND || reduce == T_POWER )
						j -= eliminate(j1+1, j2);                        // 0 (*,/,and,^) ptr = 0
					else if (reduce == T_PLUS || reduce == T_OR) {     // 0 (+,or) ptr = ptr
						j -= eliminate(j1, j1+1); 
						indStack1[shm1] = j2 - 2; }
					//else if (reduce == T_MINUS)                      // NO MORE UNARY MINUS!!  // 0 - ptr = - ptr
					//  j -= eliminate(j1, j1);
				}
				else if (t1.val == 1.0) {                            // *** Check 1.0
					if (reduce == T_MULT || reduce == T_AND) {         // 1 (*,and) ptr = ptr
						j -= eliminate(j1, j1+1);
						indStack1[shm1] = j2 - 2; }
					else if (reduce == T_OR && pr_reduce <= pr_old)    // 1 or ptr = 1 (check priority of preceeding op!)
						j -= eliminate(j1+1, j2);
				}
			}    // endif first operand is a number

			//if ( op == BR_CLOSE ) return;
			if ( j < j0 ) stackHead --;
			else break;
		}
		else break;  // end if priority[ reduce = opStack[stackHead - 1] ] <= pr

	opStack[ stackHead++ ] = op;
}

//**************************************************************************

void ExpressionObj::Optimize(int *firstTerm)
{
	struct TermStruct term;
	int    Type, unaryOp = 0;
	int    *opStack   = new int[term_num];
	int    *indStack0 = new int[term_num+1];
	int    *indStack1 = new int[term_num+1];
	int    j, j0, jold = 0, stackHead = 0;

	if (term_num <= 0 ) { delete [] opStack; delete [] indStack0; delete [] indStack1; return; }
	if (firstTerm) j0 = *firstTerm;
	else j0 = 0;

	for (j = j0; j < term_num; j ++) { // *** LOOP OVER TERMS BEGINS HERE ***

		Type = (term = term_array[j]).type;

		if ( isUnary(Type) ) {   // Catch unary sign operators here
			unaryOp = term.type; continue;
		}

		if ( isBinary(Type) ) {  // ...So this binary operator is not a unary operator
			pushOp( Type, opStack, indStack0, indStack1, stackHead, j); 
			continue;
		}

		if ( Type == BR_CLOSE )  break; // catch a closing bracket
		// *** The term is a number:
		jold = j;
		if ( isFunction(Type) ) {         // 1. a function call (or nested expression)
			j ++;
			Optimize(&j); 
			if (jold + 2 == j  && ( term_array[jold+1].type == NUMBER_TYPE || Type == BRACKET ) ) { // kill a function call
				j -= eliminate(jold, jold);
				if (Type != BRACKET) term_array[jold].val = function( Type, term_array[jold].val );
				j -= eliminate(jold+1, jold+1);
				Type = term_array[jold].type;
			}
		}
		else if ( Type != NUMBER_TYPE && Type != POINTER_TYPE)  
			throw makeMessage("Could not simplify expression: unknown term type %d", Type);

		if ( unaryOp != 0 && Type == NUMBER_TYPE) {     // kill a unary operator
			if  (unaryOp == T_UNARY_NOT)  term_array[j].val = (term_array[j].val > 0) ? 0 : 1;
			// else if (unaryOp == T_UNARY_MINUS)  term_array[j].val = -term_array[j].val;
			j -= eliminate(j-1, j-1);
		}
		unaryOp = 0;

		indStack0[ stackHead ] = jold;
		indStack1[ stackHead ] = j;
	}

	pushOp( BR_CLOSE, opStack, indStack0, indStack1, stackHead, j); // *** Collect the result

	if (term_array[j0].type == BRACKET && j0 == jold && term_array[j-1].type == BR_CLOSE) // Simplify "( expression ) = expression"
	{ eliminate(j-1,j-1); eliminate(j0,j0); j-=2; }

	if (firstTerm) *firstTerm = j;

	delete [] opStack; delete [] indStack0; delete [] indStack1;
}



//****************************************************************************************************


void ExpressionObj::print(FILE *f)
{
	fprintf(f,"%s ", formula);
	fflush(f);
}


//**************************************************************************

int get_max_term_num(TokenString &Param, long p)
{
	int terms = 0;

	while ( !isExpressionStop( Param[p] ) && ( p < Param.token_num ) ) 
	{ p++; terms++; }

	return terms;
}

//**************************************************************************************************

ExpressionObj::ExpressionObj(TokenString &Param, long p0, const char *message,  
class VarList *VL, long *pLast,
	double *arg1, const char *arg1id, 
	double *arg2, const char *arg2id, 
	double *arg3, const char *arg3id) {

		Allocate( 2 * get_max_term_num(Param, p0) );
		if ( ++callLevel > MAX_CALL_LEVEL )
			Param.errorMessage(p0, 0, "Check for possible recursive definition (e.g. \"a = cos( b ); b = 4 a\")");

		while (Param.equal(p0, ASSIGN_TOKEN)) p0 ++; 

		int tp;
		long p = p0;
		int term = 0;
		int bracket_num = 0;

		if ( term_num == 0 && pLast == 0 ) Param.errorMessage(p0, ::StrCpy("Empty expression:"), message);

		while ( !isExpressionStop( Param[p] ) )  { // *************  Main expression loop  ***************

			if (VERBOSE > 6)  fprintf(stderr,"expression: term(%d) = %s \n", term , Param[p] );

			for (tp = 0; tp < TYPE_NUM; tp++)
				if (Param.equal(p, type_id[tp])) {  // it's an operator/function token
					if (VERBOSE > 6)  fprintf(stderr,"it's an operator token = %s\n",type_id[tp]);
					term_array[term].type = tp; 

					if ( Param.equal(p, "-") )  { // convert a unary "-" to binary, if appropriate; ELIMINATED UNARY "+"
						if ( term == 0 ) {
							term_array[term].type = NUMBER_TYPE; term_array[term].val = 0;  // insert a zero for unary minus
							term_array[++term].type = tp;
						}
						else if ( !isOperand(term_array[term-1].type) ) {
							if ( isFunction(term_array[term-1].type) || 
								(term_array[term-1].type >= 6 && term_array[term-1].type <= TOP_BINARY) ) { // here's why things like (var < -5) do not work
									term_array[term].type = NUMBER_TYPE; term_array[term].val = 0;  // insert a zero for unary minus
									term_array[++term].type = tp;
							}
							else term_num = term;
						}
					}

					if ( isBinary(tp) )  {  // catch a binary operator as first token, or two binaries in a row
						if (term == 0)  term_num = term; 
						else if ( !isOperand(term_array[term-1].type) ) term_num = term;
					}
					if ( isUnary(tp) && term > 0 )  // catch a unary operator following an operand
						if ( isOperand(term_array[term-1].type) ) term_num = term;

					if ( tp == BR_CLOSE ) 
						if ( --bracket_num < 0 ) term_num = term;  // bracket mismatch

					if ( tp == BRACKET ||  isFunction( tp ) ) bracket_num ++;

					if ( isFunction( tp ) && term > 0 ) 
						if ( isOperand( term_array[term-1].type ) ) { // add a multiplication
							term_array[term + 1].type = tp;
							term_array[term].type = T_MULT;
							term ++;
						}
					break;
				}

			if ( tp >= TYPE_NUM) { // did not match any token

				if (Param.isConst( Param[p], &(term_array[term].val)) ) {  // it's a number
					if (VERBOSE > 6)  fprintf(stderr,"it's a number = %g\n", term_array[term].val );
					term_array[term].type = NUMBER_TYPE;
				}
				else if ( equal(Param[p], arg1id) ) { term_array[term].ptr = arg1; term_array[term].type = POINTER_TYPE; }
				else if ( equal(Param[p], arg2id) ) { term_array[term].ptr = arg2; term_array[term].type = POINTER_TYPE; }
				else if ( equal(Param[p], arg3id) ) { term_array[term].ptr = arg3; term_array[term].type = POINTER_TYPE; }
				else if ( VL ) 
					if (( term_array[term].ptr = VL->ResolveID( Param[p] ) )) {  // it's a pointer // memory leak!!
						if (VERBOSE > 6)  fprintf(stderr,"it's a pointer = %p\n", &term_array[term].ptr);
						term_array[term].type = POINTER_TYPE;
					}
					else  term_num = term;  // does not match anything
				else term_num = term;  // there's no VarList argument

				if ( term > 0 ) if ( isOperand( term_array[term-1].type ) ) { // add a multiplication
					term_array[term + 1] = term_array[term];
					term_array[term].type = T_MULT;
					term ++;
				}
			}

			if (term_num <= term) break;
			p++; term++;
			if ( p >= Param.token_num ) { term_num = term - 1; break; }

		}  // ******************************   Main expression loop 

		if ( isExpressionStop( Param[p] ) )  term_num = term; 

		if ( bracket_num > 0 ) Param.errorMessage(p-1, ::StrCpy("Error before closing bracket"), message);

		if (term_num) {                           // Run some more checks on the last expression:
			// expression can only end in a number or a bracket    
			while ( !isOperand( term_array[term_num-1].type ) && (term_num > 1) )  term_num--;

			char last = term_array[term_num-1].type;
			if (term_num == 1)                     // one-term expression is either a pointer or a number
				if ( (last != NUMBER_TYPE ) && (last != POINTER_TYPE ) ) term_num = 0;
		} 

		if (pLast) *pLast = p - 1;
		else if ( !isExpressionStop( Param[p] ) ) Param.errorMessage(p, ::StrCpy("Cannot parse the expression:"), message);

		// ********************************************************************************************** 
		//  if (VERBOSE > 6) fprintf(stderr,"FORMULA=%d", term_array[0].type); 

		int j=0, jold;

		do {
			jold = j;
			j = 0;                   // Second pass through the optimizer:
			try { Optimize(&j); } 
			catch (char *str) { Param.errorMessage(p0+j-1, str, message); } 
		} while (jold != j);


		if (j < term_num) Param.errorMessage(p0+j-1, ::StrCpy("Bad expression (possibly a bracket nesting error)"), message);
		buildFormulaString(VL, arg1, arg1id, arg2, arg2id, arg3, arg3id, message);
		if (VERBOSE > 6) fprintf(stderr,"\nExpression = { %s}\n", formula);

		Evaluate();
		callLevel--;
}


//****************************************************************************************************


void ExpressionObj::buildFormulaString(class VarList *VL, double *p1, const char *s1, double *p2, const char *s2, double *p3, const char *s3, const char *message) 
{
	int  term;
	int  tp;

	strcpy(formula,"");

	for (term = 0; term < term_num; term++) {
	   if ( strlen(formula) > MAX_FORMULA_LENGTH ) return;

	   tp = term_array[term].type;
	   if (tp < TYPE_NUM) { 
		   if ( tp != T_MULT ) strcat(formula, type_id[tp]); else continue; 
	   }
	   else if (tp == NUMBER_TYPE) {
			if ( term_array[term].val == 0 && (term + 2) < term_num ) // process a "0 - number" combination,
			   if ( term_array[term+1].type == T_MINUS ){                   // which is a unary minus
				   if ( term == 0 ) continue;
				   else if ( isFunction(term_array[term-1].type) ) continue;
			   } 
			char value[64];
			snprintf(value, 63, "%g", term_array[term].val );
			strcat(formula, value);
	   }
	   else if (tp == POINTER_TYPE) {
		   char *str = NULL;
		   if      ( term_array[term].ptr == p1 ) str = ::StrCpy(s1);
		   else if ( term_array[term].ptr == p2 ) str = ::StrCpy(s2); 
		   else if ( term_array[term].ptr == p3 ) str = ::StrCpy(s3); 
		   else if (VL) str = ::StrCpy( VL->ResolvePtr( term_array[term].ptr ) );
		   if (str) {
				strcat(formula, str);
				delete str;
		   }
	   }
	   else throw makeMessage("An error in expression \"%s\": bad type=%d",formula,tp);
	   strcat(formula," ");
	}
}


//****************************************************************************************************


FILE *fopenAssure(const char *fname, const char *mode, const char *action, const char *id) {
	FILE *f;
	     if ( equal(fname, "stderr") )   return (FILE *)stderr;
	else if ( equal(fname, "stdout") )   return (FILE *)stdout;
	else if ( (f = fopen(fname, mode)) )   return f;
	else { perror(">>> fopen Error = ");
	throw makeMessage("Could not open file \"%s\" for %s %s", fname, action, id); }
}

//*******************************************************************************************************

char *makeMessage(const char *fmt, ...) {

	int n, size = 300;
	va_list ap;
	while (1) {
		char *p = new char[size];
		/* Try to print in the allocated space. */
		va_start(ap, fmt);
		n = vsnprintf (p, size, fmt, ap);
		va_end(ap);
		/* If that worked, return the string. */
		if (n > -1 && n < size)  return p;
		/* Else try again with more space. */
		if (n > -1)    /* glibc 2.1 */
			size = n+1; /* precisely what is needed */
		else           /* glibc 2.0 */
			size *= 2;  /* twice the old size */
		delete [] p;
	}
}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************************************
