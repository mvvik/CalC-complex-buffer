/*****************************************************************************
 *
 *                     Calcium Calculator (CalC)
 *               Copyright (C) 2001-2019 Victor Matveev
 *
 *                               peak.cpp
 *
 *  "min" / "max" extremum-tracking variable object class PeakTrackObj
 *  Array of min/max variables PeakTrackArray
 *  PeakTrackArray is a member of the "KineticObj" defined in gate.h/gate.cpp
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
#include "peak.h"


extern double get_sim_time(TokenString &);

//**************************************************************************

void PeakTrackObj::Evaluate()
   {
     if ( *tptr >= T1 && *tptr <= T2 ) {
       if ( !startFlag ) { startFlag = true; peak = *pointer; peakTime = *tptr; }
       else if ( isMax ) {  if ( *pointer >= peak )  { peak = *pointer; peakTime = *tptr; } }
            else  if ( *pointer <= peak )  { peak = *pointer; peakTime = *tptr; }
     }
   }

//*******************************************************************************

#define PEAK_TIME_ID "\001"

PeakTrackObj::PeakTrackObj(TokenString &Param, const char *minmax, int i)
   {
   isMax = equal(minmax,"max") ? true : false;
   startFlag = false;

   long p = Param.token_index(minmax, i+1) - 1;

   ID = Param.getVarName(p, "min/max variable");

   timeID = StrCpy(PEAK_TIME_ID);
   if (p > 0) 
     if ( !Param.equal(p-1, CARR_RET_TOKEN) ) {
       delete [] timeID;
       timeID = Param.getVarName(p-1, "min/max time variable");
     }

   varPos = p + 2;
   varName = Param.StrCpy( varPos );
   
   T1 = Param.get_double(p + 3);
   T2 = Param.get_double(p + 4);
   double Tmax = get_sim_time(Param);

   if (T1 < 0) T1 = 0.0;
   if (T2 > Tmax) T2 = Tmax;

   if ( T1 >= T2 ) Param.errorMessage(p + 3, 
                   makeMessage("Zero time interval length in \"%s\" object: [ %g, %g ]", minmax, T1, T2) );
   if (VERBOSE)
       fprintf(stderr,"  # %s %s = %s of %s on [%g,%g] \n", timeID, ID, minmax, varName, T1, T2);
   }


//*******************************************************************************
//*******************************************************************************
/*
void PeakTrackArray::copy(const PeakTrackArray& PTA)  
{ 
    peak_num = PTA.peak_num;
    array = new PeakTrackObj *[peak_num];

    for (int i = 0; i < peak_num; i++)
      array[i] = new PeakTrackObj(*PTA.array[i]); 
}    
*/
//*******************************************************************************

PeakTrackArray::PeakTrackArray(TokenString &Param)
{
  int min_num = Param.token_count("min");
  int max_num = Param.token_count("max");

  peak_num = min_num + max_num;

  array = new PeakTrackObj *[peak_num];

  if (peak_num && VERBOSE) 
       fprintf(stderr,"\n### Setting up %d extremum tracking variable(s):\n", peak_num);

  int i;

  for (i = 0; i < max_num; i++)  array[i] = new PeakTrackObj(Param, "max", i);
  for (i = 0; i < min_num; i++)  array[i + max_num] = new PeakTrackObj(Param, "min", i);

}

//*******************************************************************************


void PeakTrackObj::set_pointer(TokenString &Param, class VarList *VL)
{

  if ( ! ( pointer = VL->ResolveID( varName, &tptr ) ) )
      Param.errorMessage(varPos, 0, "Cannot track min/max of an undefined variable");
  
}

//*******************************************************************************//*******************************************************************************


void PeakTrackArray::set_pointers(TokenString &Param, class VarList *VL)
{

  for (int i = 0; i < peak_num; i++ )  array[i]->set_pointer(Param, VL);
  
}

//*******************************************************************************
