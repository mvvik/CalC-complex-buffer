/**************************************************************************
 *
 *                      Calcium Calculator (CalC)
 *               Copyright (C) 2001-2019 Victor Matveev
 *
 *                               loop.h
 *
 *  LoopObj variable contains data for the (nested) "for" loop script statement
 *
 **************************************************************************
 
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

#ifndef CALC_LOOP_H_included
#define CALC_LOOP_H_included

#define MAX_FOR_STEPS 1000000

class LoopObj
{
public:

  char loopVarString[ 2048 ];

  PlotArray *plots;
  VectorObj *result;

  int       num;
  long      steps;
  long      count;
  double    fCount;
  bool      started;
  VectorObj var, var0, dvar;
  int       *index, *limit;
  char      **ids, **trackIDs, *tp;

  LoopObj(TokenString &, VectorObj *);
  ~LoopObj();
 
  int  operator()() { return num; }
  void step();
  void draw();
};

int  getTrackVarNum(TokenString &);
long getTrackVar(TokenString &TS, int n, int * = 0, int * = 0);

#endif
