/*****************************************************************************
 *
 *                          Calcium Calculator (CalC)
 *                  Copyright (C) 2001-2019 Victor Matveev
 *
 *                                  table.h
 *
 * "table" construct TableObj: imported two-column time dependent data file
 *  Array of tables class "TableArray"
 *
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

#ifndef CALC_TABLE_H_included
#define CALC_TABLE_H_included

#define TABLE_TOKEN  "table"

extern int VERBOSE;

//**************************************************************************

class TableObj
{
  public:

  char *file_name, *ID;
  FILE *file;

  double *xarray, *yarray;
  double value, valueOld[2];

  long    num, index, indexOld[2];

  /*
  TableObj() { file_name = ID = 0; xarray = yarray = 0; }

  TableObj(const TableObj& TO)  { copy(TO); }

  const TableObj& operator=(const TableObj& TO)  {
    if (&TO != this) { destroy(); copy(TO); }
    return *this;
  }

  void copy(const TableObj& TO);
  */

  void saveState(int l = 0) { indexOld[l] = index; valueOld[l] = value; }
  void recoverState(int l = 0) { index = indexOld[l]; value = valueOld[l]; }

  TableObj(const char *fname, char *id);

  ~TableObj()  { destroy(); }

  void destroy() { delete [] ID; delete [] xarray; delete [] yarray; }

  double Evaluate(double);
};

//**************************************************************************

class TableArray
{
 public:

  TableObj **array;
  int      table_num;

  /*
  TableArray() { array = 0; table_num = 0; };

  TableArray(int n) : 
    table_num(n), array(new TableObj *[n]) { }

  TableArray(const TableArray& TA)  { copy(TA); }

  const TableArray& operator=(const TableArray& TA)  {
    if (&TA != this) { destroy(); copy(TA); }
    return *this;
  }

  void copy(const TableArray& TA);
  */

  TableArray(TokenString &Param);

  ~TableArray()  {  destroy(); } 


  void Evaluate(double t)  {
    for (int i = 0; i < table_num; i++) array[i]->Evaluate(t);
  }

  void destroy() { for (int i = 0; i < table_num; i++) delete array[i];
                   delete [] array; }

  void saveState(int l = 0) { for (int i = 0; i < table_num; i++) array[i]->saveState(l); }
  void recoverState(int l = 0) { for (int i = 0; i < table_num; i++) array[i]->recoverState(l); }

};

//**************************************************************************

#endif
