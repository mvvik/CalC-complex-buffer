/*****************************************************************************
 *
 *                      Calcium Calculator (CalC)
 *               Copyright (C) 2001-2022 Victor Matveev
 *
 *                                gate.cpp
 *
 *  Variable of the class KineticObj contains and updates (using 5th order 
 *  adaptive Runge-Kutta) all the scalar (non-field) time-dependent quantities:
 *
 *   o ODE variables (VectorObj *var)
 *   o Auxiliary user-defined variables (VectorObj *ident)
 *   o Array of formulas defining the latter two (ExpressionObj  **formulas)
 *   o Array of imported time-dependent "table" file data (TableArray *tables)
 *   o Extremum tracking min/max variables (PeakTrackArray *maxima)
 *   o Array of interpolated local and global concentration variables (InterpolArray *locations)
 *       
 *  KineticObj variable is a member of the main SimulationObj object, defined
 *  in simulation.h / simulation.cpp
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

#include "PlatformSpecific.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>  // for compatibility with Visual C++
#include <stdarg.h>
#include "syntax.h"
#include "vector.h"
#include "box.h"
#include "grid.h"
#include "field.h"
#include "table.h"
#include "peak.h"
#include "markov.h"
#include "interpol.h"
#include "simulation.h"
#include "fplot.h"
#include "gate.h"

//**************************************************************************

  KineticObj::~KineticObj()  
    {
    int i;

    delete maxima;

    for (i=0; i < ident_num; i++) {
      delete ident_id[i];
      delete formulas[i];
      }

    for (i=0; i < var_num; i++)  {
      delete var_id[i];
      delete formulas[i+ident_num];
      }

    delete var;
    delete varOld[0];
    delete varOld[1];
    delete [] var_id;

    delete ident;
    delete [] ident_id;
    delete [] formulas;

    delete switches;
    delete locations;
    delete tables;
 }

//**************************************************************************

void KineticObj::init_ODE_AUX(TokenString *Param)
  {
    int i, pos;
    int var_check  = Param->token_count("dt");
    int var_num1   = Param->token_count(ODE_TOKEN);
    int var_num2   = Param->token3_count("/","dt","=");

    if ( var_check > var_num2 )    // Check for invalid use of reserved keyword "dt" 
      for (i=0; i<var_check; i++) {
	    pos = Param->token_index("dt", i+1);
	    if ( !Param->equal(pos-1,"/") || !Param->equal(pos+1,ASSIGN_TOKEN) ) 
	       Param->errorMessage(pos, 0, "Reserved word \"dt\" can only be used in ODE definitions (df/dt = expression)");
        }                          // End check

    var_num   = var_num1 + var_num2;
    ident_num = Param->token_count(VAR_TOKEN);

    var       = new VectorObj(var_num);
    varOld[0] = new VectorObj(var_num);
    varOld[1] = new VectorObj(var_num);
    var_id    = new char *[var_num];

    ident     = new VectorObj(ident_num);
    ident_id  = new char *[ident_num];

    formulas  = new ExpressionObj *[var_num + ident_num];
  }

//**************************************************************************


void KineticObj::sortAUX(class TokenString *Param, long *exprStart)
{
   int  j;
   long pos, p0, p1;

  long count = ident_num * ident_num * 2;  // max number of variable swaps
  int ii;

  for (ii=ident_num-1; ii > 0; ii--)  { // Sorting the identities according to inter-dependence
     j = 0; 
     while (j < ii)    {
        p0 = exprStart[j];
        p1 = Param->lastInLine(p0);
        for (pos = p0; pos <= p1; pos ++)
     	  if (equal(ident_id[ii], (*Param)[pos]) )   {
	         char *temp    = ident_id[ii];
             long exprPos  = exprStart[ii];
             ident_id[ii]  = ident_id[j];
             ident_id[j]   = temp;
             exprStart[ii] = exprStart[j];
             exprStart[j]  = exprPos;
             if ( --count <= 0 ) 
				 Param->errorMessage(p0, 0, "Recursive time-dependent variable definitions");
             break;
	     }
        if ( pos > p1 ) j++;   // made it through the whole loop without problem
        else j=0;              // otherwise, start all over 
    }
  }  

} 
//**************************************************************************


void KineticObj::setNames_ODE_IC(class TokenString *Param, long *exprStart)
{
  int  i;
  long pos;
  int var_num1   = Param->token_count(ODE_TOKEN);

  for (i = 0; i < var_num; i++) {
    if (i < var_num1) {
      pos = Param->token_index(ODE_TOKEN, i + 1);
      var_id[i] = Param->getVarName( pos - 1, "dynamic ODE variable");
    }
    else {
      pos = Param->token3_index("/","dt","=", i + 1 - var_num1);
      var_id[i] = Param->getVarName( pos - 3, "dynamic ODE variable", 1);
    }
    exprStart[ident_num + i] = pos + 1;
    (*var)[i] = 0.0;
    if ( (Param->token3_count(var_id[i], "(", "0", &pos )) )  {
      if ( Param->equal(pos + 1,")")  && Param->equal(pos + 2, "=") ) 
        (*var)[i] = ExpressionObj(*Param, pos + 3, "Bad expression for the initial condition").Evaluate();
      else Param->errorMessage(pos + 1, 0, "Bad expression for the initial condition"); 
	}
  }
}
//**************************************************************************


void KineticObj::setNames_AUX_Sort(class TokenString *Param, long *exprStart)
{
  int   i;
  long  pos;

  for (i = 0; i < ident_num; i++) {
    pos = Param->token_index(VAR_TOKEN, i+1);
    ident_id[i]  = Param->getVarName( pos - 1, "time-dependent variable");
    exprStart[i] = pos + 1;
  }

  sortAUX(Param, exprStart);
}
//**************************************************************************

void KineticObj::setAUXexpressions(class SimulationObj *Sim, long *exprStart)
{
  int   i;
  TokenString *Param = Sim->Params;

  if (VERBOSE && ident_num)
    fprintf(stderr, "\n### Setting up %d time-dependent (auxiliary) variable(s):\n", ident_num);

  for (i = 0; i < ident_num; i++)  {
    if (VERBOSE) fprintf(stderr,"  #%d: %s %s ", i+1, ident_id[i], VAR_TOKEN );
    formulas[i] = new ExpressionObj(*Param, exprStart[i], "Bad auxiliary variable definition", Sim, 0, &Time, "t");
    ident->elem[i] = formulas[i]->Evaluate();
    if (VERBOSE) { formulas[i]->print(stderr);
                   fprintf(stderr,"; %s(0) = %g\n",  ident_id[i], ident->elem[i] ); }
    }
}
//**************************************************************************

void KineticObj::setODEexpressions(class SimulationObj *Sim, long *exprStart)
{
  int   i;
  TokenString *Param = Sim->Params;

  if (VERBOSE && var_num) fprintf(stderr, "\n### Setting up %d ODE(s):\n",var_num);

  for (i = 0; i < var_num; i++)  {
    if (VERBOSE) fprintf(stderr,"  #%d: d%s/dt = ", i+1, var_id[i] );
    formulas[i + ident_num] = 
        new ExpressionObj(*Param, exprStart[i + ident_num], "Bad ODE expression", Sim, 0, &Time, "t");

    if (VERBOSE) { formulas[i + ident_num]->print(stderr);
                   fprintf(stderr,"; %s(0) = %g\n", var_id[i], (*var)[i] ); } 
    formulas[i + ident_num] -> Evaluate(); // make it crash here if expression is not valid
    }
}

//**************************************************************************


KineticObj::KineticObj(class SimulationObj *Sim)
{
   TokenString *Param = Sim->Params;
 
  if (Sim->Ca) {
    Param->addToken("_Charge"); Param->addToken("`");
    Param->addToken("="); Param->addToken("_ICa"); Param->addToken(";"); 
	Ca = Sim->Ca;
  }
  else Ca = 0;

  init_ODE_AUX(Param);                 
  long *exprStart = new long[var_num + ident_num]; 

  Time = 0.0;

  locations = new InterpolArray(*Param, Sim->Ca, Sim->Buffers);
  tables    = new TableArray(*Param);
  setNames_ODE_IC  (Param, exprStart);
  setNames_AUX_Sort(Param, exprStart);
  maxima    = new PeakTrackArray(*Param);
  switches  = new MarkovArray(*Param);    
 
  Sim->Gates = this;

  maxima -> set_pointers(*Param, Sim);
  setAUXexpressions(Sim, exprStart);
  setODEexpressions(Sim, exprStart);

  switches -> set_matrices(*Param, Sim);
  switches -> Evaluate(0.0);

  Evaluate();
  Equilibrate(Param);
};

//**************************************************************************

void KineticObj::Equilibrate(TokenString *Param) {

  if ( Param->token_count("equilibrate") && var_num)
    {
      eq_flag = 1; 

      double eps = KINETIC_DEFAULT_EPS, eq_time = 1.0;
      long maxSteps = 40000;

      Param->trail_pars("equilibrate", 1, 'l', &maxSteps, 'd', &eps, 'd', &eq_time);

      if (VERBOSE) fprintf(stderr, 
        "\n### Equilibrating the ODE(s) (steps = %ld, precision = %g, initial step = %g ms)",
         maxSteps, eps, eq_time);

      RunStatusString  *rss = 0;
      if (VERBOSE > 2)  rss = new RunStatusString(30, "", 0, 0);
      else if (VERBOSE) rss = new RunStatusString(1, "", 0, 0);
      eq_time = RungeKuttaAdaptive(eq_time, eps, 0, rss, 0.0, 0.0, maxSteps);
      if (VERBOSE) fprintf(stderr, "\n  # ODE(s) equilibrated after %g ms \n", eq_time);
      Time = 0.0;
    }

  eq_flag = 0;
  Evaluate();

  if (VERBOSE) {
	  fprintf(stderr, "\n     >> ");
	  for (int i = 0; i < ident_num; i++) fprintf(stderr," %s(0)=%g ",  ident_id[i], ident->elem[i] );
	  fprintf(stderr, "\n     >> ");
	  for (int i = 0; i < var_num; i++)   fprintf(stderr," %s(0)=%g ", var_id[i], (*var)[i] );
	  fprintf(stderr, "\n");
  }
};

//**************************************************************************

void KineticObj::Evaluate()
 { 
 VectorObj dvar(var_num);
 int i;

 if (!eq_flag)
   {
   locations->Evaluate(Time);
   tables->Evaluate(Time);
   }

 for (i = 0; i < ident_num; i++)    ident->elem[i] = formulas[i]->Evaluate();
 maxima -> Evaluate();
 for (i = 0; i < ident_num; i++)    ident->elem[i] = formulas[i]->Evaluate();  // In case aux depend on max/min
 //if (Ca) Ca->evaluateCurrents(); // Fix the bootstrap in the order of object declaration vs. evaluation
 }

//**************************************************************************

void KineticObj::print()
 { 
 VectorObj dvar(var_num);
 int i;

 fprintf(stderr, "\n\n Variables:  "); for (i = 0; i < var_num;   i++)  fprintf(stderr, "%s=%g ", this->var_id[i],   this->var->elem[i]);
 fprintf(stderr, "\n\n Identities: "); for (i = 0; i < ident_num; i++)  fprintf(stderr, "%s=%g ", this->ident_id[i], this->ident->elem[i]);
 }

//**************************************************************************

VectorObj KineticObj::derivative()
 { 
 VectorObj dvar(var_num);

 Evaluate();
 if (Ca && !eq_flag) Ca->evaluateCurrents();

 for (int i = 0; i < var_num;   i++) dvar[i] = formulas[i+ident_num]->Evaluate();

 return dvar;
 }


//**************************************************************************
//      Cash-Karp embedded adaptive step-size Runge-Kutta scheme 
//                 (see Numerical Recipes, p. 710)
//**************************************************************************

double KineticObj::RungeKuttaAdaptive(double T, double eps, class PlotArray *plots, class RunStatusString *rss, 
                                      double ceiling, double dt, long maxSteps) 
{
  double Time0 = Time, eqTime = 0.0, dthold;
  long   count = 0;
  double NORMeps = eps, norm = 0.0;

  if ( eq_flag )  { 
	dt = T; T = 1e20;  
	NORMeps = 10*eps;   // accuracy for the norm equilibrium condition should be less 
                        // than the accuracy of the difference method
  } else {
    if ( T <= 0 ) return Time;
    else if (dt <= 0) dt = T;    // dt is the initial time step, T is the integration time, ceiling is the 
    T += Time;                   //                                                  total simulation time
    if ( ceiling )  {
                    if ( Time >= ceiling )  return Time;
                    else if ( T > ceiling )  T = ceiling;
                    }
	if ( !var_num ) { // if there are no ODEs to integrate, update all the time variables
					Time = T;
					Evaluate();
					return Time;
					}
    }   // end if ( equilibrating )

  static double safety = 0.95;
  static const double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;

  static const double b21 = 0.2;
  static const double b31 = 3.0 / 40.0, b32 = 9.0/40.0;
  static const double b41 = 0.3, b42 = -0.9, b43 = 1.2;
  static const double b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0;
  static const double b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, b64 = 44275.0/110592.0, 
					  b65 = 253.0/4096.0;

  static const double c1  = 37.0 / 378.0,     c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0;
  static const double cc1 = 2825.0 / 27648.0, cc3 = 18575.0 / 48384.0, cc4 = 13525.0 / 55296.0, cc5 = 277.0/14336.0,
                      cc6 = 0.25;

 VectorObj v1(var_num), dv(var_num); 
 VectorObj k1(var_num), k2(var_num), k3(var_num), k4(var_num), k5(var_num), k6(var_num);

 double delta, L1norm;

 do
   {
   if ( Time + dt > T )  dt = T - Time;
   saveState(1);

   k1 = derivative();
   *var = *varOld[1] + dt * (b21 * k1 );
   Time = TimeOld[1] + a2 * dt;
   k2 = derivative();
   *var = *varOld[1] + dt * (b31 * k1 + b32 * k2 );
   Time = TimeOld[1] + a3 * dt;
   k3 = derivative();
   *var = *varOld[1] + dt * (b41 * k1 + b42 * k2 + b43 * k3 ); 
   Time = TimeOld[1] + a4 * dt;
   k4 = derivative();
   *var = *varOld[1] + dt * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4 ); 
   Time = TimeOld[1] + a5 * dt;
   k5 = derivative();
   *var = *varOld[1] + dt * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5 ); 
   Time = TimeOld[1] + a6 * dt;
   k6 = derivative();   

   v1 = c1  * k1 +  c3 * k3 +  c4 * k4            +  c6 * k6;
   dv = cc1 * k1 + cc3 * k3 + cc4 * k4 + cc5 * k5 + cc6 * k6;

   dv -= v1;
   L1norm = dv.L_infty_norm();
   delta = dt * L1norm; 

   if (VERBOSE > 8) { dv.print();  fprintf(stderr,"\n"); }

   if (delta > eps || !_finite(L1norm) )   // accuracy condition not satisfied
     {
     recoverState(1);
     if (VERBOSE > 7) fprintf(stderr, "# RK step: delta=%.3e > eps=%.3e time=%g dt=%.3g->", delta, eps, Time, dt);
     dthold = 0.1 * dt;
     if ( !_finite(L1norm) ) dt /= 5.0; 
     else dt *= safety * pow(delta / eps, -0.25);
     if (dt < dthold) dt = dthold;
	 while (Time + dt <= Time)  dt *= 2.0;
     // if ( dt < Time * DT_MIN ) 
     //  throw makeMessage("Runge-Kutta time step reduced below %g, at time %g\n", DT_MIN * Time, Time);
     if (VERBOSE > 7) fprintf(stderr, "%.3g\n", dt);
     else if (rss) rss->update('\"');

     continue;
     }
   else   // accuracy eps is satisfied
     {
     Time = TimeOld[1] + dt;
     *var = *varOld[1] + v1 * dt;
     Evaluate();   // ensure that all variables are evaluated at the current time 
     //switches->Evaluate(Time);
     if (plots) plots->draw_all();
     if (VERBOSE > 7) fprintf(stderr, "# RK step: delta=%.3e < eps=%.3e time=%g dt=%.3g->", delta, eps, Time, dt);
     dthold = 5.0 * dt;
     if (delta > 0) dt *= safety * pow(delta / eps, -0.2); 
     if (dt > dthold) dt = dthold;
     if (VERBOSE > 7) fprintf(stderr, "%g\n", dt); 
     if (rss) rss->update('\'');
     }

   if ( eq_flag ) { 
       eqTime += Time; Time = Time0;
      if ( (norm = derivative().norm() ) < NORMeps || ++count > maxSteps ) break; }
   }
 while ( (T - Time) / (Time + T + 1e-12) > 1e-12 );

 switches->Evaluate(Time);

 if (count > maxSteps)  {
     fprintf(stderr,"\n@@@ KineticObj equilibration: more than %ld iterations (norm = %.3g, step = %.3g ms) %c\n", maxSteps, norm, dt, 7);
     fprintf(stderr,  "@@@ Change the equilibration step or equilibration accuracy\n");
     fprintf(stderr,  "@@@ Equilibration may not work if ODEs depend on time explicitly\n");
     fprintf(stderr,  "    or if the ODE system is stiff\n");
 }

 if (VERBOSE > 7) fprintf(stderr, "## RK run completed: dt=%.3g\n", dt);

 if (eq_flag) return eqTime;
 return Time;
}

//*****************************************************************************

double *KineticObj::ResolveID(const char *svar, double **t) {

  int j;

  // if ( isFloat(svar) ) throw makeMessage("Expected a variable name, not a constant \"%s\"",svar);

    if (VERBOSE > 8) fprintf(stderr, ">> ResolveID: Resolving %s\n", svar);

    if (t) *t = &Time;

    if (var) 
        for (j = 0; j < var->get_size(); j++) 
          if (equal(svar, var_id[j])) 
            return (*var)() + j;

    if (ident)
        for (j = 0; j < ident->get_size(); j++)
          if (equal(svar, ident_id[j])) 
            return  (*ident)() + j;

    for (j = 0; j < tables->table_num; j++)
	 if (equal(svar, tables->array[j]->ID ) )
           return &(tables->array[j]->value);

    for (j = 0; j < maxima->peak_num; j++)
	 if (equal(svar, maxima->array[j]->ID ) ) 
           return &(maxima->array[j]->peak);
         else  if (equal(svar, maxima->array[j]->timeID ) ) 
           return &(maxima->array[j]->peakTime);

    for (j = 0; j < locations->interpol_num; j++)
		if (equal(svar, locations->array[j]->ID ) ) {
			if (t) *t = locations->array[j]->tptr;
			return &(locations->array[j]->result);
		}

    for (j = 0; j < switches->number; j++) 
		if (equal(svar, switches->array[j]->ID ) )  {
			if (t) *t = &(switches->array[j]->currentTime);
	        return &(switches->array[j]->dstate);
		}

   return 0;
}

//**************************************************************************


const char *KineticObj::ResolvePtr(double *ptr) {

   int j;

      if (var) 
        for (j = 0; j < var->get_size(); j++)
			 if ( ptr == ((*var)() + j) ) return var_id[j];
          
      if (ident)
        for (j = 0; j < ident->get_size(); j++)
			if ( ptr == ((*ident)() + j) ) return ident_id[j];

      for (j = 0; j < tables->table_num; j++)
		 if ( ptr == &(tables->array[j]->value) ) return tables->array[j]->ID;

      for (j = 0; j < maxima->peak_num; j++)
		if ( ptr == &(maxima->array[j]->peak) ) return maxima->array[j]->ID;
         else if ( ptr == &(maxima->array[j]->peakTime) ) return maxima->array[j]->timeID;

      for (j = 0; j < locations->interpol_num; j++)
		 if ( ptr == &(locations->array[j]->result) ) return locations->array[j]->ID;

      for (j = 0; j < switches->number; j++)
		if ( ptr == &(switches->array[j]->dstate) ) return switches->array[j]->ID;

    return 0;
}

//***************************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

