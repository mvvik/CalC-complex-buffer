/*****************************************************************************
 *
 *                     Calcium Calculator (CalC)
 *              Copyright (C) 2001-2019 Victor Matveev
 *
 *                               peak.h
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

#ifndef CALC_PEAK_H_included
#define CALC_PEAK_H_included

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************


class PeakTrackObj
{
protected:

  double T1, T2;

public:

 double *pointer;
 double peak, peakTime, peakStored[2], peakTimeStored[2];
 bool   startFlag, startFlagStored[2];
 double *tptr;
 char   *ID,  *timeID;
 char   *varName;
 long   varPos; //char   MinMax[4];
 bool   isMax; // true if tracking a maximum; false if tracking a minimum
    
 void set_pointer(TokenString &Param, class VarList *VL);

 PeakTrackObj(TokenString &Param, const char *minmax, int i);

 ~PeakTrackObj() { delete [] ID; delete [] timeID; delete [] varName; }

 void Evaluate();

 void saveState(int level = 0) 
   { peakStored[level] = peak;  peakTimeStored[level] = peakTime; startFlagStored[level] = startFlag; }

 void recoverState(int level = 0) 
   { peak = peakStored[level];  peakTime = peakTimeStored[level]; startFlag = startFlagStored[level]; }
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************

class PeakTrackArray
{
 public:

  PeakTrackObj **array;
  int          peak_num;

  PeakTrackArray()  { array = 0; peak_num = 0; }

  PeakTrackArray(TokenString &);

  void set_pointers(TokenString &TS, class VarList *VL);

  ~PeakTrackArray()  { for (int i=0; i<peak_num; i++) delete array[i];
                       if (peak_num) delete [] array; }  

  void Evaluate()  {
    for (int i = 0; i < peak_num; i++) array[i]->Evaluate();
  }

  void saveState(int level = 0) { for (int i = 0; i < peak_num; i++) array[i]->saveState(level); }
  void recoverState(int level = 0) { for (int i = 0; i < peak_num; i++) array[i]->recoverState(level); }
};

//**************************************************************************

#endif
