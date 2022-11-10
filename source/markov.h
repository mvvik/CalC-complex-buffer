/*****************************************************************************
 *
 *                    Calcium Calculator (CalC)
 *             Copyright (C) 2001-2019 Victor Matveev
 *
 *                              markov.h
 *
 *  Markov variable object (discrete state / continuous time)
 *
 *  Array of markov process objects - class "MarkovArray"
 *  MarkovArray is a member of the "KineticObj" defined in gate.h/gate.cpp
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

#ifndef CALC_MARKOV_H_included
#define CALC_MARKOV_H_included

#define MY_SEED 293847944

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************


class MarkovObj
{
protected:

  double storedTime[2];
  int    storedState[2];

  struct TermStruct *transitions;  /* the transition generator matrix */

public:

  static double *errorTolerance;  /* bound to m_ODEaccuracy in SimulationObj constructor */
  int    States; /* number of discrete states */
  int    state;  /* current state of the markov variable (an integer) */
  double dstate; /* same, but float value */
  double currentTime;
  char   *ID;

  MarkovObj() { ID = 0; transitions = 0; }

  MarkovObj(TokenString &Param, int ind);

 ~MarkovObj() { if (transitions) delete [] transitions; if (ID) delete [] ID; }

 double Evaluate(double t);

 void   saveState(int level=0) { storedState[level] = state; storedTime[level] = currentTime; }

 void   recoverState(int level=0) 
 { state = storedState[level]; dstate = double(state); currentTime = storedTime[level]; }

 void   set_matrix(TokenString &Param, class VarList *VL);
};

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************

class MarkovArray
{
 public:

  MarkovObj  **array;
  int        number;

  MarkovArray(TokenString &);

  ~MarkovArray()                   { for (int i = 0; i < number; i++) delete array[i]; delete [] array; }

  void Evaluate(double t)          { for (int i = 0; i < number; i++) array[i]->Evaluate(t); }
  void saveState   (int level = 0) { for (int i = 0; i < number; i++) array[i]->saveState(level);    }
  void recoverState(int level = 0) { for (int i = 0; i < number; i++) array[i]->recoverState(level); }

  void set_matrices(TokenString &Param, class VarList *VL);
};

//**************************************************************************

#endif
