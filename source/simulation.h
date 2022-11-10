/*****************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *                Copyright (C) 2001-2019 Victor Matveev
 *
 *                              simulation.h
 *
 *  Defines the main (highest-level) simulation compound object SimulationObj, 
 *  which contains the concentration fields (FieldObj "Ca" object and the 
 *  BufferArray object) and the remaining scalar time-dependent variables 
 *  combined in the KineticObj variable "Gates". Also links the grid object 
 *  (GridObj *Grid), the geometry object (RegionObj *Synapse), the plot array
 *  variable (PlotArray *Plots), and the parsed script (TokenString *Params).
 *
 *  Includes as members the highest-level adaptive "Adaptive(double)" and 
 *  non-adaptive "FixedTimeStep(double)" time-update routines.
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

#ifndef CALC_SIMULATION_H_included
#define CALC_SIMULATION_H_included

typedef void (*CaMethod)       (FieldObj &,    VectorObj &,   BufferArray &, BufferArray &, double);
typedef void (*BufMethod)      (BufferArray &, BufferArray &, VectorObj &,   VectorObj &,   double); 

void Buf1DstepNew(BufferObj &Buf, VectorObj &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
void Buf2DstepNew(BufferObj &Buf, VectorObj &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
void Buf3DstepNew(BufferObj &Buf, VectorObj &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
 
void Buf1DstepCoop(BufferArray &Buf, BufferArray &BubNew, VectorObj   &Ca,  VectorObj   &CaNew,  double dt);
void Buf2DstepCoop(BufferArray &Buf, BufferArray &BubNew, VectorObj   &Ca,  VectorObj   &CaNew,  double dt);
void Buf3DstepCoop(BufferArray &Buf, BufferArray &BubNew, VectorObj   &Ca,  VectorObj   &CaNew,  double dt);
 
void Ca1DstepCoop(FieldObj  &Ca,  VectorObj   &CaNew,  BufferArray &Buf, BufferArray &BufNew, double dt);
void Ca2DstepCoop(FieldObj  &Ca,  VectorObj   &CaNew,  BufferArray &Buf, BufferArray &BufNew, double dt);
void Ca3DstepCoop(FieldObj  &Ca,  VectorObj   &CaNew,  BufferArray &Buf, BufferArray &BufNew, double dt);

//*******************************************************************************************

class SimulationObj : public VarList {

 protected:

  double       m_dt0, m_accuracy, m_dtStretch, m_ODEaccuracy, m_dtMax;
  char         ERROR_FLAG;

 public:

  double  totalSimTime;

  //double  ***React;  // 2D array of pointers to reaction rates stored in Gates 
  double  **DiffArray; // array of different diffusibility fields
  int     DiffNum;     // number of distinct diffusibility fields

  CaMethod  CaStep;
  BufMethod BufStep;

  class TokenString *Params;
  class BCarrayObj  *BCArray;
  class RegionObj   *Synapse;
  class GridObj     *Grid;
  class FieldObj    *Ca;
  class BufferArray *Buffers;
  class KineticObj  *Gates;
  class PlotArray   *Plots;
  class VectorObj   *kuptake;
  
  void initialize() { Params = 0; BCArray = 0; Synapse = 0; Grid = 0; Ca = 0; Buffers = 0; Gates = 0;
                      DiffArray = 0; DiffNum = 0; DiffArray = 0; Plots = 0; kuptake = 0; ERROR_FLAG = 1; };

  SimulationObj()  { initialize(); }
  SimulationObj(TokenString &TS);
 ~SimulationObj();

  void initTortuosity();
  void killTortuosity();

  void Export(const char *filename);
  void Import(const char *filename);
  
  double         *ResolveID   (const char *, double **t=0);
  const char     *ResolvePtr  (double *ptr);
  class FieldObj *ResolveField(const char *);

  virtual void Run();
  void FixedTimeStep(double T, int n);
  long Adaptive(double T);
};

//*******************************************************************************************

class ODESimulationObj : public SimulationObj {

 public:

   ODESimulationObj(TokenString &TS);
  ~ODESimulationObj() { }

  void Run();                
};

//*******************************************************************************************
#endif
