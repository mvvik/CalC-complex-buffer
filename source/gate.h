/*****************************************************************************
 *
 *                      Calcium Calculator (CalC)
 *               Copyright (C) 2001-2019 Victor Matveev
 *
 *                                gate.h
 *
 *  Variable of the class KineticObj contains and updates (using 5th order 
 *  adaptive Runge-Kutta) all the scalar (non-field) time-dependent quantities:
 *
 *   ODE variables (VectorObj *var)
 *   Auxiliary user-defined variables (VectorObj *ident)
 *   Array of formulas defining the latter two (ExpressionObj  **formulas)
 *   Array of imported time-dependent "table" file data (TableArray *tables)
 *   Extremum tracking min/max variables (PeakTrackArray *maxima)
 *   Array of interpolated local and global concentration variables (InterpolArray *locations)
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

#ifndef CALC_GATE_H_included
#define CALC_GATE_H_included

#define  KINETIC_DEFAULT_EPS 1.0e-6
#define  DT_MIN 1.0e-16   // smallest time step (multiplied by current time)

//*************************************************************************************

class KineticObj : public VarList
{

 public:

  VectorObj       *var, *varOld[2];
  VectorObj       *ident;
  int             var_num;
  int             ident_num;
  char            **var_id;
  char            **ident_id;
  ExpressionObj   **formulas;
  FieldObj        *Ca;

  TableArray      *tables;
  PeakTrackArray  *maxima;
  InterpolArray   *locations;
  MarkovArray     *switches;
  
  double Time, TimeOld[2];

  //  KineticObj(int n, int m) : var_num(n), ident_num(m), Time(0.0) { initialize(); }

  KineticObj(class SimulationObj *);

  void init_ODE_AUX(TokenString *);

  void sortAUX(TokenString *, long *exrStart);

  void setNames_ODE_IC  (class TokenString *Param, long *exprStart);
  void setNames_AUX_Sort(class TokenString *Param, long *exprStart);

  void setAUXexpressions(class SimulationObj *Sim, long *exprStart);
  void setODEexpressions(class SimulationObj *Sim, long *exprStart);

  void Evaluate();
  void Equilibrate(TokenString *);
  void print();

  VectorObj derivative();                      // returns dt satisfying accuracy eps
  double    RungeKuttaAdaptive(double T, double eps, class PlotArray *plots = 0, class RunStatusString *rss=0,
                               double total = 0.0, double dt = 0.0, long max_steps = 40000);   

  char      eq_flag;

  ~KineticObj();

  void saveState(int level = 0)    { *varOld[level] = *var; TimeOld[level] = Time; 
     switches->saveState(level); tables->saveState(level); maxima->saveState(level); locations->saveState(); }

  void recoverState(int level = 0) { *var = *varOld[level]; Time = TimeOld[level]; 
     switches->recoverState(level); tables->recoverState(level); maxima->recoverState(level); locations->recoverState(); }

  double     *ResolveID(const char *, double **t=0);
  const char *ResolvePtr(double *ptr);
};


//*************************************************************************************


#endif
