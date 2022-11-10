/*****************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *                Copyright (C) 2001-2019 Victor Matveev
 *
 *                               interpol.h
 *
 *  Interpolation/concentration averaging object "InterpolObj" 
 *  (appear in script as "Ca[x,y,z]" or "Ca[]")
 *
 *  Array of interpolation objects class "InterpolArray"
 *  InterpolArray is a member of the "KineticObj" defined in gate.h/gate.cpp
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

#ifndef CALC_INTERPOL_H_included
#define CALC_INTERPOL_H_included

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************


class InterpolObj
{
protected:

  FieldObj *Field;
  bool     average;  // true if InterpolObj simply tracks a field integral over volume

  int      nPointers;
  double **pointers;
  double  *factors;

  double oldres,  newres,  oldtime,  newtime;
  double oldres0, newres0, oldtime0, newtime0;

  char   **token_ptr_source, *token_ptr_target;  // restore 

public:

 double *tptr;
 double result;
 char   *ID;

 InterpolObj() { ID = 0; pointers = 0; }
 
 InterpolObj(TokenString &Param, long p, FieldObj *);

 ~InterpolObj() { if (!average) *token_ptr_source = token_ptr_target;
				  if (pointers) delete [] pointers;
                  if (ID)       delete [] ID; }

 double Evaluate(double t);
 void   reset();

 void saveState()    { oldres0 = oldres; newres0 = newres;  oldtime0 = oldtime; newtime0 = newtime; }
 void recoverState() { oldres  = oldres0; newres = newres0; oldtime = oldtime0; newtime  = newtime0; }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************

class InterpolArray
{
 public:

  InterpolObj  **array;
  int          interpol_num;

  InterpolArray(TokenString &, FieldObj *Ca, BufferArray *Bufs);

  ~InterpolArray()       { for (int i = 0; i < interpol_num; i++) delete array[i]; delete [] array; }

  void Evaluate(double t){ for (int i = 0; i < interpol_num; i++) array[i]->Evaluate(t);   }
  void reset()           { for (int i = 0; i < interpol_num; i++) array[i]->reset();       }
  void saveState()       { for (int i = 0; i < interpol_num; i++) array[i]->saveState();   }
  void recoverState()    { for (int i = 0; i < interpol_num; i++) array[i]->recoverState();}
};

//**************************************************************************

#endif
