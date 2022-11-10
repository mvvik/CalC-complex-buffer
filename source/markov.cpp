/*****************************************************************************
 *
 *                     Calcium Calculator (CalC)
 *             Copyright (C) 2001-2019 Victor Matveev
 *
 *                              markov.cpp
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

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "syntax.h"
#include "markov.h"

double *MarkovObj::errorTolerance = 0;

//**************************************************************************

double MarkovObj::Evaluate(double Time)
   {
     int    steps, index;
     double Tstep, deltaT, maxProb = 0.0, prob;

     if ( (deltaT = Time - currentTime) <= 0) return Time;
     index = state * (States - 1);

     for (int j=0; j<States; j++) {
       if (j == state) continue;
       if ( (prob = ::Evaluate(transitions[index])) > maxProb )  maxProb = prob;
       index++;
     }

     if (maxProb <= 0.0) return (currentTime = Time);
     
     Tstep = sqrt(*errorTolerance / maxProb);
     steps = ( Tstep >= deltaT ? 1 : int(deltaT / Tstep) );

     for (int subdivisions = 0; subdivisions < steps; subdivisions++) {
       currentTime += Tstep;
       double random = double(rand()) / (double(RAND_MAX) + 1.0);
       index = state * (States - 1);
       double sum = 0;
       for (int j=0; j<States; j++) {
         if (j == state) continue;
         sum += Tstep * ::Evaluate(transitions[index]);
         if (sum > random) { dstate = double(state = j); return Evaluate(currentTime); } 
         index++;
       }
     }
  
     return 0;   
   }

//*******************************************************************************

MarkovObj::MarkovObj(TokenString &Param, int i)
   {
   long pos = Param.token_index("markov", i+1) - 1;

   ID = Param.getVarName(pos, "markov variable");

   States = Param.get_int(pos+2);
   state  = Param.get_int(pos+3);
   dstate = double(state);
   currentTime = 0.0;

   if (VERBOSE)
       fprintf(stderr,"  # %s = markov variable with %d states; %s(0)=%d \n", ID, States, ID, state );

   }



//*******************************************************************************

MarkovArray::MarkovArray(TokenString &Param)
{

  number = Param.token_count("markov");

  array = new MarkovObj *[number];

  if (number && VERBOSE) 
       fprintf(stderr,"\n### Setting up %d Markov variable(s):\n", number);

  for (int i = 0; i < number; i++)  array[i] = new MarkovObj(Param, i);
}

//*******************************************************************************


void MarkovObj::set_matrix(TokenString &Param, class VarList *VL)
{
  char   *element = 0;
  int    index = 0;
  double value, *ptr;


  transitions = new TermStruct[States*(States-1)];

  for (int i=0; i<States; i++)
    for (int j=0; j<States; j++) {
      if (i == j) continue;
      element = makeMessage("%s.%d.%d",ID,i,j);
      if ( Param.isConst(element, &value) ) {
        transitions[index].type = NUMBER_TYPE;
        transitions[index].val  = value;
      }
      else if ( (ptr = VL->ResolveID(element)) ) {
        transitions[index].type = POINTER_TYPE;
        transitions[index].ptr = ptr;
      }
      else  {   // is zero if undefined
        transitions[index].type = NUMBER_TYPE;
        transitions[index].val  = 0.0;
      }
      delete [] element;       
      index ++;
    }
}


//*******************************************************************************


void MarkovArray::set_matrices(TokenString &Param, class VarList *VL)
{

  for (int i = 0; i < number; i++ )  array[i]->set_matrix(Param, VL);
  
}

//*******************************************************************************
