/*****************************************************************************
 *
 *                        Calcium Calculator (CalC)
 *                 Copyright (C) 2001-2022 Victor Matveev
 *
 *                                loop.cpp
 *
 *  LoopObj variable contains data for the (nested) "for" loop script statement
 *
 ******************************************************************************
 
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

 ******************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "PlatformSpecific.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "vector.h"
#include "syntax.h"
#include "box.h"
#include "field.h"
#include "table.h"
#include "peak.h"
#include "interpol.h"
#include "markov.h"
#include "gate.h"
#include "fplot.h"
#include "loop.h"

extern char* globalLabelX;

//**********************************************************************************************

int getTrackVarNum(TokenString &TS) {
  int num = 0;
  for (int i=0; i<TS.token_count("track"); i++) {
    num += TS.tokens_to_eol( TS.token_index("track", i+1) + 1);
  }
  return num;
}

//**********************************************************************************************

long getTrackVar(TokenString &TS, int n, int *iGraph, int *iSets) {  // n = 0..track var number - 1
  int num = 0, m;
  for (int i=0; i<TS.token_count("track"); i++) {
    if (iGraph) *iGraph = i;
    long p = TS.token_index("track", i+1) + 1; 
    num += ( m = TS.tokens_to_eol( p ) );
    if (iSets) *iSets =  m;
    if (num >= n + 1) return p + m - ( num - n );
  }
  throw makeMessage("Could not find track variable number %d", n);
}

//**********************************************************************************************

LoopObj::LoopObj(TokenString &TS, VectorObj *res) : num( TS.token_count(LOOP_TOKEN) ), 
						    var(num), var0(num), dvar(num),
                            index(new int[num]), limit(new int[num]), 
                            ids(new char *[num]),  tp(new char[num])
  {
  count   = 0;
  fCount  = 0.0;
  started = false;
  steps   = 1;

  if (!num) {
      plots = new PlotArray(0);
      return;   
  }                                // ***  Below this line: only if "for" loop present ****

  char prefix[256], filename[256];
  double v0, v1, dv = 1.0;
  int i;

  strcpy(prefix, "");
  result   = res;
  trackIDs = new char *[res->size + 1];

  for (i = 0; i < res->size; i++) trackIDs[i] = TS.StrCpy( getTrackVar(TS, i) );
  
  for (i = 0; i < num; i++) 
  try    {       
      ids[i] = TS.StrCpy( TS.token_index(LOOP_TOKEN, i+1) + 1 ); 
      v0 = ExpressionObj(TS, TS.token3_index(LOOP_TOKEN, ids[i],ASSIGN_TOKEN) + 1, 
                             "Error in the \"for\" loop: bad loop index initialization").Evaluate();
      v1 = ExpressionObj(TS, TS.token_index("to", i+1) + 1, 
                             "Error in the \"for\" loop: bad \"to\" expression").Evaluate();
      dv = ExpressionObj(TS, TS.token_index("step", i+1) + 1, 
                             "Error in the \"for\" loop: bad \"step\" expression").Evaluate();

      if (  ( fabs(v0 - int(v0+0.5)) < 1.0e-6 ) && 
            ( fabs(v1 - int(v1+0.5)) < 1.0e-6 ) &&  
            ( fabs(dv - int(dv+0.5)) < 1.0e-6 )  )   tp[i] = 'i'; else tp[i] = 'd';

      var.elem[i] = var0.elem[i] = v0;
      dvar.elem[i] = dv;
   
      int n = int( (v1 - v0) / dv + 0.5 ) + 1;
      if (n > MAX_FOR_STEPS) n = MAX_FOR_STEPS;
      limit[i] = n;
      index[i] = 0;
      steps *= n;    

      if (VERBOSE) fprintf(stderr, "\n### For loop: for %s = %g to %g step %g\n", ids[i], v0, v1, dv);
  } 
  catch (char* str) { TS.errorMessage(TS.token_index(LOOP_TOKEN, i + 1), str); }
  catch (int ERR)   { TS.errorMessage(TS.token_index(LOOP_TOKEN, i + 1), 0,
                                      makeMessage("ERROR %d: Bad \"for\" loop statement", ERR));
                    }

  PlotObj::UPDATE_STEPS = steps;

  if ( TS.token2_count("plot.method","xmgr") || TS.token3_count("plot.method","=","xmgr") ) {
     plots = new PlotArray(result->size + TS.token_count("track"));
     plots->method = METHOD_XMGR;
     XmgrPlot::init(TS, TS.token_count("track"));
  }
#ifndef _NO_GLUT_
  else if (TS.token2_count("plot.method", "gl") || TS.token3_count("plot.method", "=", "gl")) {
      plots = new PlotArray(2 * result->size);
      plots->method = METHOD_GL;
      int graphs = result->size;
      GlPlotObj::init(TS, graphs, 1, 1);
      plots->gl_on = 1;
  }
#endif
  else {
     plots = new PlotArray(result->size);
     plots->method = METHOD_MUTE;
  }
 
  double* xPtr, xMax; // the pointer to x-label variable for all "track" graphs, and its max value

  if (num == 1) {
      globalLabelX = StrCpy(ids[0]);
      xPtr = var();
      xMax = v0 + (limit[0] - 1) * fabs(dvar[0]);
  }
  else {
      globalLabelX = StrCpy("iteration");
      xPtr = &fCount;
      xMax = double(steps);
  }

  long pos;
  if (TS.token_count("plot.print", &pos)) TS.line_string(pos + 1, prefix);

  int xmgrGraph, xmgrGraph0 = -1, graphNum = 0, xmgrSets;

  for (i = 0; i < result->size; i++) 
    {
    if (VERBOSE) fprintf(stderr, "    tracking variable #%d: %s\n", i+1, trackIDs[i]);
    snprintf(filename, 255, "%s%s", prefix, trackIDs[i]);

    if ( plots->method == METHOD_MUTE )
      plots->set_plot( new MutePointPlot( (*result)()+i, xPtr, 0, xMax, filename, trackIDs[i]) );
    else if ( plots->method == METHOD_XMGR ) {
      getTrackVar(TS, i, &xmgrGraph, &xmgrSets);
      if (xmgrGraph > xmgrGraph0) {
         plots->set_plot( new XmgrPointPlot(xmgrSets, 0, xMax) );
         graphNum = ( xmgrGraph0 =  xmgrGraph ) + i;
      }
      ((XmgrPointPlot *)(plots->array[graphNum]))->set_set( (*result)()+i, xPtr, trackIDs[i]);
      plots->set_plot( new MutePointPlot( (*result)()+i, xPtr, 0, xMax, filename, trackIDs[i]) );
    }
#ifndef _NO_GLUT_
    else if (plots->method == METHOD_GL) {
        plots->set_plot(new GlPointPlot((*result)() + i, xPtr, 0, xMax, trackIDs[i]));
        plots->set_plot(new MutePointPlot((*result)() + i, xPtr, 0, xMax, filename, trackIDs[i]));
    }
#endif // !_NO_GLUT_
  }

  }

//**************************************************************************************************
//**************************************************************************************************

  LoopObj::~LoopObj()
  {
   int i;
   if (!num) return;

   if (plots && plots->plot_num) delete plots;
   delete [] tp;
   delete [] index;
   delete [] limit;

   for (i=0; i<result->size; i++) delete trackIDs[i];
   for (i=0; i<num; i++)          delete ids[i];

   delete [] ids; delete [] trackIDs;
 
  }


//**************************************************************************************************

void LoopObj::step() {
  
  int i;
  char temp[1024];
  strcpy( loopVarString, "");

  if (!num) return;

  if ( !started ) started = true;
  else {
    count ++;    // increment cycle vars 
    fCount += 1;
    for (i = num - 1; i >= 0; i--) {
      if ( ++index[i] < limit[i] ) break;
      else index[i] = 0; 
    }
  }

  for (i = 0; i < num; i++)
    {       
    var[i] = var0[i] + index[i] * dvar[i];
    if (tp[i] == 'i')
      snprintf(temp, 1023, "%s = %d", ids[i], int(var[i] + 0.5) );
    else
      snprintf(temp, 1023, "%s = %g", ids[i], var[i] );

    strcat( loopVarString, temp );
    if ( i < num - 1 ) strcat( loopVarString, "; " );
    }
 
  if (VERBOSE > 0) {
    fprintf(stderr,  "\n    Next step in the \"for\" loop: %s \n", loopVarString);
    fprintf(stderr,  "========================================================\n");
  }
}

//**************************************************************************************************

  void LoopObj::draw()
    {
    int i;
    if (!num) return;
    plots->draw_all();  

    if (VERBOSE > 0 && result->size) {
      fprintf(stderr,"\n========================================================\n ");
      fprintf(stderr,"    %s : ", loopVarString);
      for (i = 0; i < result->size; i++) fprintf(stderr," %s=%g", trackIDs[i], (*result)[i]);
      fprintf(stderr,"\n========================================================\n");
    }
    fflush(stdout);
    }

//**************************************************************************************************
