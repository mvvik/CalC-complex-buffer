/*****************************************************************************
 *
 *                         Calcium Calculator (CalC)
 *                 Copyright (C) 2001-2019 Victor Matveev
 *
 *                             simulation.cpp
 *
 *  Defines the main highest-level simulation compound object SimulationObj, 
 *  which contains the concentration fields (FieldObj "Ca" object and the 
 *  BufferArray object) and the remaining scalar time-dependent variables 
 *  combined in the KineticObj variable "Gates". Also links the grid object 
 *  (GridObj *Grid), the geometry object (RegionObj *Synapse), the plot array
 *  variable (PlotArray *Plots), and the parsed script (TokenString *Params).
 *
 *  Includes as members the highest-level adaptive "Adaptive(double)" and 
 *  non-adaptive "FixedTimeStep(double)" time-update routines. 
 *
 *  "Ca" on one side, and the buffers plus KineticObj variables on the other
 *  side, are updated on temporal grids staggered by dt/2
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

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>  // for compatibility with Visual C++
#include <string.h>
#include <stdarg.h>
#include "vector.h"
#include "syntax.h"
#include "box.h"
#include "grid.h"
#include "field.h"
#include "table.h"
#include "peak.h"
#include "markov.h"
#include "interpol.h"
#include "simulation.h"
#include "gate.h"
#include "fplot.h"

extern int    Number_Of_Iterations_Per_PDE_Step;   
extern double CHARGE_LOSS;

extern double get_sim_time(TokenString &TS);
extern void   setFieldByFunction(TokenString &TS, long p, double *Diff, const char *errStr);
extern const  char *getMethod(CaMethod &, BufMethod &);
extern void   getRun( TokenString &TS, int i, bool *adaptive, double *time=0, double *accuracy=0,
					  double *dtMax=0, double *dt=0, double *dtStretch=0, double *ODEaccuracy=0);

//*******************************************************************************************

SimulationObj::SimulationObj(TokenString &TS) {

  initialize();

  long pos;
  char temp[512];

  totalSimTime = get_sim_time(TS); // make sure simulation time is computed before Plot Array is constructed!
  if (VERBOSE) fprintf(stderr, "\n ### Total Simulation Time = %g \n", totalSimTime);

  strcpy(LABEL_DIM1,""); strcpy(LABEL_DIM2,""); strcpy(LABEL_DIM3,"");

  if      (TS.assert("geometry","spherical") )      GEOMETRY = SPHERICAL; 
  else if (TS.assert("geometry","conical") )        GEOMETRY = CONICAL; 
  else if (TS.assert("geometry","cylindrical") )    GEOMETRY = CYLINDRICAL;
  else if (TS.assert("geometry","polar") )          GEOMETRY = POLAR;
  else if (TS.assert("geometry","cartesian.1D") )   GEOMETRY = CARTESIAN1D;
  else if (TS.assert("geometry","cartesian.2D") )   GEOMETRY = CARTESIAN2D;
  else if (TS.assert("geometry","spherical.3D") )   GEOMETRY = SPHERICAL3D;
  else if (TS.assert("geometry","cylindrical.3D") ) GEOMETRY = CYLINDRICAL3D;
  else if (TS.assert("geometry","cartesian.3D") ) { }
  else if (TS.assert("geometry","disc.1D") || TS.assert("geometry","disc") || TS.assert("geometry","disk")) GEOMETRY = DISC;
  else if (TS.token_count("geometry"))  TS.errorMessage( TS.token_index("geometry") + 1, 0, "Unknown geometry specifier");
                                               else GEOMETRY = CARTESIAN3D;

  DIMENSIONALITY =  GEOMETRY & DIMENSION_MASK;
  GEOMETRY1 = (GEOMETRY >> 2) & DIMENSION_MASK;
  GEOMETRY2 = (GEOMETRY >> 4) & DIMENSION_MASK;
  GEOMETRY3 = (GEOMETRY >> 6) & DIMENSION_MASK;

  switch (GEOMETRY1) { case 0: strcpy(LABEL_DIM1,"x"); break;
                       case 1: strcpy(LABEL_DIM1,"r"); break;
                       case 2: strcpy(LABEL_DIM1,"r"); break;
                      }
  if (DIMENSIONALITY > 1) 
  switch (GEOMETRY2) { case 0: strcpy(LABEL_DIM2,"y");     break;
                       case 1: strcpy(LABEL_DIM2,"z");     break;
                       case 2: strcpy(LABEL_DIM2,"phi");   break;
                       case 3: strcpy(LABEL_DIM2,"theta"); break;
                      }
  if (DIMENSIONALITY > 2) 
  switch (GEOMETRY3) { case 0: strcpy(LABEL_DIM3,"z");   break;
                       case 1: strcpy(LABEL_DIM3,"phi"); break;
                      }
  
  if (VERBOSE) fprintf(stderr, "\n ### %dD Geometry, variables: %s %s %s \n", DIMENSIONALITY, LABEL_DIM1, LABEL_DIM2, LABEL_DIM3);

  Params  = &TS;
  BCArray = new BCarrayObj(TS);
  Synapse = new RegionObj(TS);
  Grid    = new GridObj( *Synapse, TS);
  FieldObj::setStaticData(*Synapse, *Grid, *BCArray); 
  Synapse->bindToGrid( *Grid );
  Synapse->computeFormulas(TS);
  Ca      = new FieldObj(TS);
  Ca->adjust_sources(*Ca);  // Increase current if some of it falls outsise the boundary
 kuptake = new VectorObj(Grid->Size);
  *kuptake = 0.0;
  Ca->kuptake = this->kuptake;

  if ( TS.token_count("uptake", &pos) )  {
    if (VERBOSE) fprintf(stderr, "\n### Distributed Ca2+ uptake (clearance) defined \n");
    setFieldByFunction(TS, pos + 1, kuptake->elem, "Bad uptake definition expression");
  }
  // if (!Ca->source_num) throw StrCpy("No calcium channels specified");
  Buffers = new BufferArray(TS);
 
  initTortuosity();
  Gates   = new KineticObj(this);
  MarkovObj::errorTolerance = &m_ODEaccuracy;
  if ( TS.token_count("Import", &pos) ) Import( TS.line_string( pos + 1, temp ) );
  Plots   = new PlotArray(*this);
  Plots->draw_all();

  // have to initialize run parameters before the check loop, to check the bounds on parameter values
  // also have to re-initialize after the loop, to restore the default settings

  m_dt0 = 1.0e-3; m_dtMax = 0.1; m_dtStretch = 1.03;
  m_ODEaccuracy = 2e-5;
  m_accuracy    = 4e-4;
  
  for (int i=1; i <= TS.token_count("Run"); i++) // make sure all "run" and "current" statements are Ok
  try { 
     bool adaptive; double t;
     getRun( *Params, i, &adaptive, &t, &m_accuracy, &m_dtMax, &m_dt0, &m_dtStretch, &m_ODEaccuracy);   
     getMethod( CaStep, BufStep );
          Ca->getCurrents(*Params, i, this, &Gates->Time);
     Buffers->getCurrents(*Params, i, this, &Gates->Time);
          Ca->killCurrents();
     Buffers->killCurrents();
  } catch (char *str) {
     TS.errorMessage( TS.token_index("Run", i), str, "Cannot process this \"Run\"/\"current\" statement");
  }

  m_dt0 = 1.0e-3; m_dtMax = 0.1;
  m_ODEaccuracy = 2e-5;
  m_accuracy    = 4e-4; 
  m_dtStretch = 1.03;

  Params->get_param    ("adaptive.dt0",      &m_dt0);         Params->get_param("adaptive.dtMax",     &m_dtMax);  
  Params->get_param    ("adaptive.accuracy", &m_accuracy);    Params->get_param("adaptive.dtStretch", &m_dtStretch); 
  Params->get_param    ("ODE.accuracy",      &m_ODEaccuracy); Params->get_param("adaptive.dtMax",     &m_dtMax);         

 }

//*******************************************************************************************

SimulationObj::~SimulationObj() {

  if (ERROR_FLAG)     return;  // Avoid segmentation faults if the constructor did not finish allowcating memory
  if (Plots)          delete Plots;
  if (Gates)          delete Gates;
  if (DiffArray)      killTortuosity();
  FieldObj::kill_tridiag();      
  if (Buffers)        delete Buffers;
  if (Ca)             delete Ca;
  if (Grid)           delete Grid;
  if (Synapse)        delete Synapse;
  if (BCArray)        delete BCArray;
  if (kuptake)        delete kuptake;
}

//********************************************************************************************
//********************************************************************************************


void SimulationObj::killTortuosity() {
  if (DiffArray) {
     for (int ti = 0; ti < DiffNum; ti ++) delete [] DiffArray[ti]; 
     delete [] DiffArray;
  }
}

//*****************************************************************************

void SimulationObj::initTortuosity() {

  bool undefined = false;
  int  bi;
  double *Diff = 0;
  long l = 0;

  DiffNum = 0;
  if ( Ca->tortDefined ) DiffNum++; else undefined = true;

  for (bi = 0; bi < Buffers->buf_num; bi++) 
    if ( Buffers->array[bi]->tortDefined ) DiffNum++;
    else undefined = true;

  if ( undefined ) DiffNum ++;

  DiffArray = new double *[DiffNum];
  int ind = 0;

  if ( undefined ) {   // at least one field has undefined tortuosity
    Diff = new double[ FieldObj::Size ];
    DiffArray[ind++] = Diff; 
    if ( !Params->token_count("tortuosity") && !Params->token_count("all.tortuosity") ) 
      for (l = 0; l < FieldObj::Size; l++) Diff[l] = 1.0;
    else {
      if (VERBOSE) fprintf(stderr,"\n### Global (default) tortuosity function specified:\n");
      long p = 0;
      if ( !Params->token_count("tortuosity", &p) ) Params->token_count("all.tortuosity", &p);
      if (VERBOSE) fprintf(stderr,"\n#### Tortuosity function = ");
      setFieldByFunction(*Params, p+1, Diff, "Bad tortuosity expression");
    }
  }

  if ( Ca->tortDefined ) DiffArray[ind++] = Ca->Diff; else Ca->Diff = Diff;
  for (bi = 0; bi < Buffers->buf_num; bi++) 
    if ( Buffers->array[bi]->tortDefined ) DiffArray[ind++] = Buffers->array[bi]->Diff; 
       else Buffers->array[bi]->Diff = Diff;
}
//*****************************************************************************

double *SimulationObj::ResolveID(const char *svar, double **t) {

  if (VERBOSE > 8) fprintf(stderr, ">> ResolveID: Resolving %s\n", svar);

    if (Ca) {
     if ( equal(svar, "_ICa") )  
       { if (t) *t = &Ca->Time; return &(Ca->ICa); }
     if ( equal(svar, "Charge.loss") )  
       { if (t) *t = &Ca->Time; return &(CHARGE_LOSS); }
    }     
	if (Gates) return Gates->ResolveID(svar, t);

    return 0;
}

//**************************************************************************


const char *SimulationObj::ResolvePtr(double *ptr) {

    if (ptr == &CHARGE_LOSS) return "Charge.loss";

    if (Ca) {
       if ( ptr == Ca->elem )  return Ca->ID; 
       else if ( ptr == &(Ca->ICa) )  return "_ICa";
       else if ( ptr == &CHARGE_LOSS )  return "Charge.loss";
    } 

    if (Buffers)
      for (int j = 0; j < Buffers->buf_num; j++)
	if ( ptr == Buffers->array[j]->elem )  return  Buffers->array[j]->ID; 

    if (Gates) return Gates->ResolvePtr(ptr);
    return 0;
    //fprintf(stderr,"ResolvePtr: cannot resolve pointer %X\n", ptr );
    //exit(1);
}

//***************************************************************************

FieldObj *SimulationObj::ResolveField(const char *label) {
      FieldObj *field = 0;
      if (Ca) {  if ( equal(label, Ca->ID) ) field = Ca; }
      if (Buffers) for (int i=0; i < Buffers->buf_num; i++)
        if ( equal(label, (Buffers->array[i])->ID) ) field = Buffers->array[i];
      return field;
}


//**************************************************************************************************


void SimulationObj::Export(const char *filename)
{
  FILE *f;
  char endtag[7];
  strcpy(endtag,"endtag");
  long wrote, toWrite;

  if (VERBOSE) 
    fprintf(stderr,"\n#### Dumping all state variables into file %s\n", filename);

  f = fopenAssure(filename, "wb", "Simulation state export", "");

  if (Ca) {
	toWrite = Ca->size;
    fwrite( (void *)&toWrite, sizeof(long), 1, f);
    wrote = fwrite( (void *)(Ca->elem), sizeof(double), toWrite, f);
	if ( wrote != toWrite ) 
	  throw makeMessage("Failure to export Ca profile to file %s: wrote only %ld of total %ld elements",
	                    filename, wrote, toWrite);
  }

  if (Buffers) 
	for (int i=0; i< Buffers->buf_num; i++) {
      toWrite = Buffers->array[i]->size;
      wrote = fwrite( (void *)(Buffers->array[i]->elem), sizeof(double), toWrite, f);
  	  if ( wrote != Buffers->array[i]->size ) 
 	    throw makeMessage("Failure to export %s profile to file %s: wrote only %ld of total %ld elements",
	                      Buffers->array[i]->ID, filename, wrote, toWrite);
	}

  if (Gates) {
     int m = Gates->var_num;
     fwrite( &m, sizeof(int), 1, f );
     fwrite( (void *)(Gates->var->elem), sizeof(double), Gates->var_num, f);
  }

  fwrite( (void *)endtag, sizeof(char), strlen(endtag)+1, f);
  fclose(f);
  fflush(f);
}


//**************************************************************************************************

void SimulationObj::Import(const char *filename)
{
  FILE *f;
  char endtag[7];
 
  if (VERBOSE) 
    fprintf(stderr,"\n#### Importing all state variables from file %s\n", filename);

  f = fopenAssure(filename, "rb", "Simulation state import", "");

  long n, Size, imported;
  
  if (Ca) {
    Size = Ca->size;
    fread( (void *)&n, sizeof(long), 1, f);
    if (n != Size) 
      throw makeMessage("Can't import the concentration fields: wrong number of elements (%ld vs %ld)\n", n, Size);
	imported = fread( (void *)(Ca->elem), sizeof(double), Size, f );
	if ( imported != Size ) 
		throw makeMessage("Failure to import Ca profile from file %s: read only %ld of the total %ld elements",
		                   filename, imported, Size);
  }

  if (Buffers) 
    for (int i=0; i< Buffers->buf_num; i++) {
	  Size = Ca->size;
      imported = fread( (void *)(Buffers->array[i]->elem), sizeof(double), Size, f);
  	  if ( imported != Size ) 
	 	 throw makeMessage("Failure to import %s profile from file %s: read only %ld of the total %ld elements",
		                   Buffers->array[i]->ID, filename, imported, Size);
    }

  if (Gates) {
    int m = 0;
    fread( &m, sizeof(int), 1, f);
    if (m != Gates->var_num)  
      throw makeMessage("Can't import the dynamic variables: wrong number of ODE variables (%d vs %d)\n", m, Gates->var_num);
    fread( (void *)(Gates->var->elem), sizeof(double), Gates->var_num, f);
    Gates->locations->reset();
    Gates->Evaluate();
  }

  fread( (void *)endtag, sizeof(char), strlen("endtag")+1, f );
  if ( !equal(endtag,"endtag") )  
    throw makeMessage("in Import: file size mismatch (endtag=%s)\n", endtag);

  fclose(f);
}


//**************************************************************************
//          H I N E S / N O N - A D A P T I V E    M E T H O D                 
//**************************************************************************

void SimulationObj::FixedTimeStep(double T, int n)
{
	BufferArray  BufNew(*Buffers);
	FieldObj     CaNew(*Ca);
					   
	double Time0 = Ca->Time;
	double dt = (T - Time0) / double(n);
	
	RunStatusString *rss = 0;
	if (VERBOSE > 4)      rss = new RunStatusString(30,"time",&(Ca->Time),T);						       
	else if (VERBOSE > 1) rss = new RunStatusString(1, "time",&(Ca->Time),T);						       
	
	for (int i = 1; i <= n; i++) {
		
		Buffers->setTime( Ca->Time = Time0 + (i - 1)*dt );
		Ca->evaluateCurrents();

		for (int iter = 1; iter <= Number_Of_Iterations_Per_PDE_Step; iter++) {	
			BufStep(*Buffers, BufNew, *Ca,      CaNew,  dt);
			CaStep (*Ca,      CaNew,  *Buffers, BufNew, dt);
		}
		*Ca      = CaNew;
		*Buffers = BufNew;

		Buffers->setTime( Ca->Time = Time0 + i * dt );
		Gates->RungeKuttaAdaptive(dt, m_ODEaccuracy, 0, rss, T);
		
		if (VERBOSE > 1)  rss->update('.');
		Plots->draw_all(); 
	}
	if (rss) delete rss;
}


//**************************************************************************
//         A D A P T I V E   N O N - S T A G G E R E D   M E T H O D                 
//**************************************************************************

long SimulationObj::Adaptive(double T)  {
  
  FieldObj    oldCa(*Ca), CaNew(*Ca), CaNew2(*Ca);
  BufferArray BufNew(*Buffers), oldBuf(*Buffers);
  double      one, two, error, old_dt = m_dt0, dt = m_dt0;
  int         rvalue = 0, flag = 0;
  long        between_checks = 0, total_steps = 0, since_last_divide = 0, total_backsteps = 0;
  int         errorODE = 0;

  RunStatusString *ODEstatus = 0, *status = 0;
  if (VERBOSE > 2)    status = new RunStatusString(40,"time", &(Ca->Time), T, "dt", &dt);
  else if (VERBOSE)   status = new RunStatusString(1, "time", &(Ca->Time), T);	
  if (VERBOSE > 7) ODEstatus = status;

  Ca->evaluateCurrents(); 
  Gates->saveState(); // store variables to recover

  while(1) {  // ******************************************************************

    unsigned long i=0;
    
    while (++i <= between_checks && !rvalue) { // loop over n="between_checks" iterations
      dt *= m_dtStretch;
	  if (dt > m_dtMax) dt = m_dtMax;
	  if (Ca->Time + dt >= T) { m_dt0 = dt / 2.0; dt = T - Ca->Time; rvalue = 1; }
      Ca->evaluateCurrents();
	  Plots->draw_all();
      
      one   = 1e32; 
	  two   = Ca->checkerBoardNorm(DIMENSIONALITY, Ca->xsize, Ca->xysize);
      error = 1; 
	  int iter = 0;

	  //for (int iter = 1; iter <= Number_Of_Iterations_Per_PDE_Step; iter++) {	
	  while ( error > m_accuracy  ||  iter < Number_Of_Iterations_Per_PDE_Step) {
			BufStep(*Buffers, BufNew, *Ca,      CaNew,  dt);
			CaStep (*Ca,      CaNew,  *Buffers, BufNew, dt);
			one = two;
			two = CaNew.checkerBoardNorm(DIMENSIONALITY, Ca->xsize, Ca->xysize);
			error = (two + one) == 0.0 ? 0 : fabs( 2 * (two - one) / (1 + two + one) );
			iter ++;  total_steps ++;  
			if (iter > 1000) {
				if (VERBOSE) fprintf(stderr, "Max PDE iteration number exceeded: rel. error = %g \n", error);
                break;
			}
			//fprintf(stderr,"iter=%d:  one=%g two=%g error=%g\n", iter, one, two, error);   
	  }
	  *Ca = CaNew;  *Buffers = BufNew;
      //fprintf(stderr,"\n");

	  Buffers->setTime( CaNew.Time = ( Ca->Time += dt ) );
	  try { Gates->RungeKuttaAdaptive(dt, m_ODEaccuracy, 0, ODEstatus, T); }
	  catch (char *errorMsg) { fprintf(stderr, "%s", errorMsg); errorODE = 1; break; }
	  Plots->draw_all();
      since_last_divide ++; 
      if (status) status->update('.');
      if (VERBOSE > 6) fprintf(stderr, "\n > Adaptive step #%ld (of %ld): T=%g(%g) dt=%.3g", i, between_checks, Ca->Time, Gates->Time, dt);  
	  if (Ca->Time + dt >= T)  break;    // force stability check at end of run time 

	}  // **********************  for-loop between checks *******************************
    
	if (rvalue) { 
		Plots->draw_all(); 
		m_dt0 = dt / 2.0; 
		if (status) delete status; 
		return total_steps; // procedure return point 
	} 

    CaStep (*Ca,      CaNew2,   *Buffers, *Buffers, 0.5*dt);  // *********** Check Accuracy
	BufStep(*Buffers, BufNew,   *Ca,      CaNew2,   dt);
    CaStep (CaNew2,   CaNew,    *Buffers, BufNew,   0.5*dt);
	one = CaNew.checkerBoardNorm(DIMENSIONALITY, Ca->xsize, Ca->xysize);

	CaStep(*Ca,    CaNew,  *Buffers, *Buffers, dt);
	two = CaNew.checkerBoardNorm(DIMENSIONALITY, Ca->xsize, Ca->xysize);

    error = (two + one) == 0.0 ? 0 : fabs( 2 * (two - one) / (1 + two + one) ); 

    if (VERBOSE > 6) fprintf(stderr,"\n > accuracy check after %ld steps: T=%g dt=%g norm=[%g, %g], error=%g ", since_last_divide, Ca->Time, dt, one, two, error);
    
    if ( _isnan(error) || error > 2 * m_accuracy || errorODE ) { //  *********** Max Error Exceeded
       if (status)  status->update('<'); 
       if (VERBOSE > 4) {
         fprintf(stderr,"\n*** max tolerance exceeded (err = %.2g%% > %.2g%%) \n", error*100, 2*m_accuracy*100);
         fprintf(stderr,"\n >  Stepping back by %.3gms (to time %.6g) and halving time step to %.3ems\n", Ca->Time-oldCa.Time, oldCa.Time, dt * 0.5);
       }
       *Ca = oldCa; *Buffers = oldBuf; BufNew = oldBuf;
	   CaNew.Time = Ca->Time;
	   Gates->recoverState(); Gates->Evaluate(); Ca->evaluateCurrents();
       dt = (old_dt *= 0.5);
	   if ( Ca->Time + dt <= Ca->Time )  globalError( makeMessage("*** Error: adaptive time step is too small (dt = %.2e)", dt) );
       Plots->draw_all();
       total_backsteps++;
       between_checks = since_last_divide = 0;
       }  
    else if ( error > m_accuracy ) {   //  *********** Halving the Time Step
       dt *= 0.5;
	   if (status)  status->update('|'); 
       if (VERBOSE > 5) fprintf(stderr, "\n# dt halved @ %ld steps: t=%gms, dt=%.2ems, error=%.2e \n", since_last_divide, Ca->Time, dt, error);
       since_last_divide = 0;
	   between_checks    = 2;
       } 
    else  {                 // *********** Accuracy Within Limits (error < eps), Time Step Unchanged
	   int maxSteps = 12 - int( 10 * error / m_accuracy );

	        if ( error < 0.15 * m_accuracy  )  {       if ( between_checks < 16       )  between_checks ++;         } 
	   else if ( error < 0.75 * m_accuracy  )  {       if ( between_checks < maxSteps )  between_checks ++;         } 
	   else                                    {       if ( between_checks > maxSteps )  between_checks = maxSteps; 
	                                              else if ( between_checks < 1 )         between_checks = 1;        }
	   
	   if (error > 0.8 * m_accuracy) dt *= (1.8 - error / m_accuracy);

       oldCa = *Ca;        oldBuf = *Buffers;   // store variables to recover 
       Gates->saveState(); old_dt = dt;         // when error is exceeded
	   if (Ca->Time >= T)  rvalue = 1;
       if (status)  status->update('+');
       if (VERBOSE > 6) 
         fprintf(stderr, "\n> Error within bounds: new time = %lg, one=%lg two=%lg error=%lg\n",Ca->Time,one,two,error);
       }
    } // while 1 ************************************************************
}


//**************************************************************************
//                                R U N                
//**************************************************************************

void SimulationObj::Run() {

  bool adaptive;
  double t, T = 0.0;
  long   nIterations;
  Plots->draw_all();
  
  for (int i = 1; i <= Params->token_count("Run"); i++)  {
      
     getRun( *Params, i, &adaptive, &t );
     Ca->getCurrents(*Params, i, (VarList *)this, &Gates->Time); 
     Buffers->getCurrents(*Params, i, (VarList *)this, &Gates->Time);
     T += t;
   
     char *method = StrCpy( getMethod( CaStep, BufStep ) );

       if ( adaptive ) {
          getRun( *Params, i, &adaptive, &t, &m_accuracy, &m_dtMax, &m_dt0, &m_dtStretch, &m_ODEaccuracy);
          if (VERBOSE) {
            fprintf(stderr,"\n\n#### Running adaptive %s scheme for %g ms, current = ", method, t);
            Ca->printCurrents(stderr); fflush(stderr); // make it crash here if expression for current is not valid
            fprintf(stderr, "\n > accuracy=%g, dt0=%gms, dtMax=%g, dtStretch=%g, ODEaccuracy=%g\n",
                            m_accuracy, m_dt0, m_dtMax, m_dtStretch, m_ODEaccuracy);
			fflush(stderr);
			}
          nIterations = Adaptive(T);
		}  else {  // non-adaptive method
          getRun( *Params, i, &adaptive, &t, &m_ODEaccuracy, &m_dt0);
          nIterations = int(t / m_dt0 + 0.5);
          if (VERBOSE) {
            fprintf(stderr,"\n#### Running non-adaptive %s scheme for %g ms, current = ", method, t);
            Ca->printCurrents(stderr); fflush(stderr);  // make it crash here if expression for current is not valid
			fprintf(stderr,"[I]\n > time step = %g ms, number of iterations = %ld\n", m_dt0, nIterations);
			fflush(stderr);
          }
          FixedTimeStep(T, int(nIterations) );
		}
        Plots->redraw_all();
        delete [] method;

        double exper = calcium_gain(Ca, Buffers);
        double theor = *(this->ResolveID("_Charge"));
        CHARGE_LOSS = ( exper == theor ) ? 0.0 : 100.0 * (1 - (exper + 1e-20) / (theor + 1e-20) );
        if (VERBOSE) 
			if (CHARGE_LOSS > 0)
				 fprintf(stderr," ## charge loss = %.4g%% (net charge = %g, integrated flux = %g), %ld iterations \n",  CHARGE_LOSS, exper, theor, nIterations ); 
			else fprintf(stderr," ## charge gain = %.4g%% (net charge = %g, integrated flux = %g), %ld iterations \n", -CHARGE_LOSS, exper, theor, nIterations ); 

		Ca->killCurrents();
  } // i loop over runs

  if (VERBOSE) fprintf(stderr, "\n### Simulation completed: total time = %g ms\n",T );
}
//**************************************************************************
//                           O D E   S O L V E R               
//**************************************************************************

ODESimulationObj::ODESimulationObj(TokenString &TS) : SimulationObj() {

     long pos;
     char temp[512];

     totalSimTime = get_sim_time(TS);
     Params       = &TS;
     Gates        = new KineticObj(this);
     if ( TS.token_count("Import", &pos) ) Import( Params->line_string( pos + 1, temp ) );
     Plots        = new PlotArray(*this);
     Plots->redraw_all();

     m_dt0 = 1.0;
     m_ODEaccuracy = 1e-6;
     Params->get_param("ODE.dt0",      &m_dt0);
     Params->get_param("ODE.accuracy", &m_ODEaccuracy);
  }                                 

//*******************************************************************************************

void ODESimulationObj::Run()  {
  double time = 0;
  bool   flag;

  RunStatusString *rss = 0;
  if (VERBOSE > 2)  rss = new RunStatusString(30, "time", &(Gates->Time), 0);
  else if (VERBOSE) rss = new RunStatusString(1,  "time", &(Gates->Time), 0);

  if (VERBOSE) fprintf(stderr,"\n\n#### Running the ODEs: total simulation time = %g\n", totalSimTime);

  for (int i = 1; i <= Params->token_count("Run"); i++) {
    
    getRun(*Params, i, &flag, &time, &m_ODEaccuracy, &m_dt0);

    if (VERBOSE) fprintf(stderr," ### ODE Run #%d: time = %g, accuracy = %g, initial time-step = %g\n", 
            i, time, m_ODEaccuracy, m_dt0);

    if (rss) rss->reset( Gates->Time + time );
    Gates->RungeKuttaAdaptive(time, m_ODEaccuracy, Plots, rss, totalSimTime, m_dt0);
    if (VERBOSE) fprintf(stderr,"\n");
    Plots->redraw_all();
  }

  if (rss) delete rss;
}
//*******************************************************************************************************
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//*******************************************************************************************************
