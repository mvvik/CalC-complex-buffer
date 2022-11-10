/*****************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *                Copyright (C) 2001-2019 Victor Matveev
 *
 *                                table.cpp
 *
 * "table" construct TableObj: imported two-column time dependent data file
 *  Array of tables class "TableArray"
 *  TableObj variables are members of the "KineticObj" defined in gate.h/gate.cpp
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
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "syntax.h"
#include "table.h"

//**************************************************************************

TableObj::TableObj(const char *fname, char *id) : file_name(StrCpy(fname))
{
  ID = id;
  file = fopenAssure(fname, "r", id, "table");
  if (VERBOSE) { fprintf(stderr,"  # table %s (file \"%s\"): ", ID, fname); fflush(stderr); }

  num = index = 0;
  double x, y;
  char s[256];

  while (1)
    {
    fgets(s, 255, file);
    if (strchr(s, '#') || strchr(s, ';') || strchr(s,'%') ) continue;
    sscanf(s,"%lf %lf",&x, &y);    
    if (feof(file)) break;
    num ++;
    }

  xarray = new double[num];
  yarray = new double[num];

  rewind(file);
  num = 0;
  while (1)
    {
    fgets(s, 255, file);
    if (strchr(s, '#') || strchr(s, ';') || strchr(s,'%') ) continue;
    sscanf(s,"%lf %lf",&x, &y);    
    if (feof(file)) break;
    if (num) if (x <= xarray[num-1]) continue;
    xarray[num] = x; yarray[num] = y;
    num ++;
    }

  if (VERBOSE) fprintf(stderr,"%ld entries: first = (%g, %g), last = (%g, %g)\n",
	        num, xarray[0], yarray[0], xarray[num-1], yarray[num-1]);
  fclose(file);
  Evaluate(0.0);
}

//**************************************************************************
/*
void TableObj::copy(const TableObj& TO)
{
  num = TO.num; index = TO.index;

  ID = StrCpy(TO.ID); file_name = StrCpy(TO.file_name);

  xarray = new double[num];
  yarray = new double[num];

  for (int i = 0; i < num; i++) {
    xarray[i] = TO.xarray[i];
    yarray[i] = TO.yarray[i];
  }
}
*/
//**************************************************************************


double TableObj::Evaluate(double x)
{

  while ( (xarray[index] > x) && (index > 0) ) index--;
  if (index < num - 2) while ( (x > xarray[index+1]) && (index < num - 2) ) index++;

  value = ( yarray[index] * (xarray[index+1] - x) + yarray[index+1] * (x - xarray[index]) ) / 
      ( xarray[index+1] - xarray[index] );
  /*
  fprintf(stderr,"TableObj::Evaluate: %g %g index=%ld num=%ld last:(%g,%g) next:(%g,%g) \n",
          x, value,index,num,xarray[num-1],yarray[num-1],xarray[num-2],yarray[num-2]);
  */
  return value;
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//**************************************************************************

  TableArray::TableArray(TokenString &Param) : table_num(Param.token_count(TABLE_TOKEN))
  {
    array = new TableObj *[table_num];
    char fname[1024];

    if (table_num && VERBOSE) fprintf(stderr,"\n### Setting up %d tables:\n", table_num);

    for (int i = 0; i < table_num; i++) {
      long p0 = Param.token_index(TABLE_TOKEN, i+1) - 1;
      try
        {
        Param.line_string( p0 + 2, fname );
        char *id = Param.getVarName( p0, "table" );
        array[i] = new TableObj(fname, id);
        }
        catch (char *str) 
         { Param.errorMessage(p0, str, "Bad \"table\" definition"); }
        catch (int ERR) 
         { Param.errorMessage(p0, 0, makeMessage("ERROR=%d: Bad \"table\" definition", ERR)); }
    }
     
  }

//**************************************************************************
/*
  void TableArray::copy(const TableArray& TA) 
  {
    table_num = TA.table_num;
    array = new TableObj *[table_num];

    for (int i = 0; i < table_num; i++)   *array[i] = *TA.array[i];
  }
 */
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//**************************************************************************
