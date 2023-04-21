
/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2022 Victor Matveev
 *                    LBM/NIDDK/NIH and DMS/NJIT
 *
 *              Calcium Calculator.cpp / calc.cpp
 *
 *   The main() routine and the highest-order ADI engine procedures;
 *   tortuosity function interpreter/initializer
 *
 ************************************************************************
 
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
#define GL_SILENCE_DEPRECATION

#include "PlatformSpecific.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "syntax.h"
#include "box.h"
#include "grid.h"
#include "vector.h"
#include "field.h"
#include "table.h"
#include "peak.h"
#include "markov.h"
#include "interpol.h"
#include "simulation.h"
#include "fplot.h"
#include "gate.h"
#include "loop.h"
#include "time.h"

double  CHARGE_LOSS;
int     Number_Of_Iterations_Per_PDE_Step = 3;

int     DIMENSIONALITY = 3;
int     VERBOSE        = 3;
int     GEOMETRY       = CARTESIAN3D;
int     GEOMETRY1      = 0;
int     GEOMETRY2      = 0;
int     GEOMETRY3      = 0;
char    LABEL_DIM1[2];
char    LABEL_DIM2[6];
char    LABEL_DIM3[4];

PlotArray* GluPlotArray = NULL;
char*      globalLabelX = NULL;
char*      versionStr = StrCpy("7.10.7");
char       scriptFileName[1024];

#define EXTRA_PARAM_STRING  "; pA=5.182134 ; pi=4 atan(1) ; "

//********************************************************************************************

void getRun( TokenString &TS, int i, bool *adaptive=0, double *time=0, 
             double *accuracy=0, double *dtMax=0, double *dt=0, double *dtStretch=0, double *ODEaccuracy=0);

const char *getMethod(CaMethod &, BufMethod &);
double get_sim_time(TokenString &params);

void setFieldByFunction(TokenString &TS, long p, double *Diff, const char *errStr);
void setFieldByFunction(TokenString &TS, long p, bool   *Diff, const char *errStr);

//**************************************************************************************************

void Buf1DstepNew(BufferObj &Buf, VectorObj &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
void Buf2DstepNew(BufferObj &Buf, VectorObj &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
void Buf3DstepNew(BufferObj &Buf, VectorObj &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);

void Buf1DstepCoop(BufferArray &Buf, BufferArray &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
void Buf2DstepCoop(BufferArray &Buf, BufferArray &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);
void Buf3DstepCoop(BufferArray &Buf, BufferArray &BubNew, VectorObj &Ca, VectorObj &CaNew, double dt);

void Ca1DstepCoop(FieldObj &Ca, VectorObj &CaNew, BufferArray &Buf, BufferArray &BufNew, double dt);
void Ca2DstepCoop(FieldObj &Ca, VectorObj &CaNew, BufferArray &Buf, BufferArray &BufNew, double dt);
void Ca3DstepCoop(FieldObj &Ca, VectorObj &CaNew, BufferArray &Buf, BufferArray &BufNew, double dt);

//**************************************************************************************************

void header() {
     fprintf(stderr,"\n******************************************************************");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*  Calcium Calculator (CalC)  *  version 7.10.7  *  Apr 5, 2023  *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*                  Victor Matveev, 2001-2023                     *");
     fprintf(stderr,"\n*   CalC is distributed under GPLv3: see attached license file   *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*  Dept of Math Sciences, New Jersey Institute of Technology     *");
     fprintf(stderr,"\n*                     and LBM, NIDDK, NIH                        *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n*  Supported in part by NSF DMS0417416, DMS0817703, DMS1517085   *");
     fprintf(stderr,"\n*                                                                *");
     fprintf(stderr,"\n****************************************************************\n");
}


 //**************************************************************************************************
 //                                           M A I N
 //**************************************************************************************************


 int main(int argc, char **argv) {

 char fname[1024];

 try {

	 if (argc >1) strcpy(scriptFileName, argv[1]);
	 else {
		 FILE* f = fopen("DefaultScript.txt", "r");
		 if (f) {
			 fclose(f);
			 strcpy(scriptFileName, "DefaultScript.txt");
		 } else {
			 header();
			 fprintf(stderr, "\n\n Enter the CalC script file name: ");
			 fflush(stderr);
			 scanf("%s", fname);
			 size_t i = strlen(argv[0]);
			 while (argv[0][i] != '/' && argv[0][i] != '\\' && i > 0) i--;
			 strncpy(scriptFileName, argv[0], i + 1);
			 scriptFileName[i + 1] = 0;
			 fprintf(stderr, "\n Full path = %s \n", scriptFileName);
			 strcat(scriptFileName, fname);
			 fprintf(stderr, " File name = %s \n\n", scriptFileName);
		 }
	 }

	 TokenString* TS;
	 TS = new TokenString(scriptFileName, EXTRA_PARAM_STRING, "", argc, argv);
	 TS->get_int_param("verbose", &VERBOSE);
	 if (VERBOSE) header();
#ifndef _NO_GLUT_
	 if (TS->Assert("plot.method", "gl"))  glutInit(&argc, argv);
#endif

	 VectorObj result(getTrackVarNum(*TS));
	 LoopObj   vary(*TS, &result);
	 long seed = long(time(NULL));
	 seed = TS->get_long_param("seed");
	 srand(seed);
	 delete TS;

	 // ++++++++++++++++++++++++++++ MAIN LOOP ++++++++++++++++++++++++++++

	 for (int step = 0; step < vary.steps; step++) {
		 vary.step();
		 TS = new TokenString(scriptFileName, EXTRA_PARAM_STRING, vary.loopVarString, argc, argv);
		 if (TS->token_count("exit"))     { if (VERBOSE) fprintf(stderr, "\n > Exit command: breaking execution <\n ");     delete TS; break; }
		 if (TS->token_count("continue")) { if (VERBOSE) fprintf(stderr, "\n > Continue command: breaking execution <\n "); delete TS; continue; }
		 TS->get_int_param("verbose", &VERBOSE);
		 SimulationObj* Simulation;
		 double* trackPtr;

		 if (!TS->token_count("Run")) {                      //######## CALCULATOR MODE #########
			 if (VERBOSE) fprintf(stderr, "***** No simulation runs specified; running in calculator mode *****\n");
			 if (vary())
				 for (int i = 0; i < result.size; i++)  result[i] = TS->get_double(getTrackVar(*TS, i));
			 TS->printResults(0);
		 }
		 else {                                             //####### NOT CALCULATOR MODE ######
			 if (!(TS->token_count("Ca.D") || TS->token_count("grid")) || equal(TS->get_string("mode"), "ODE")) {
				 if (VERBOSE) fprintf(stderr, "***** Running in ODE-only mode *****\n");
				 Simulation = new ODESimulationObj(*TS);     //#######  ODE solver mode: ########
			 }
			 else Simulation = new SimulationObj(*TS);       //########  FULL PDE MODE  #########
			 TS->get_int_param("Number_Of_Iterations_Per_PDE_Step", &Number_Of_Iterations_Per_PDE_Step);
			 Simulation->Run();                              //#### RUN THE DIFFERENCE SCHEME ####
			 if (vary())
				 for (int i = 0; i < result.size; i++) {
					 if ((trackPtr = Simulation->ResolveID(vary.trackIDs[i]))) result[i] = *trackPtr;
					 else  TS->errorMessage(TS->token_index(vary.trackIDs[i]), 0, "Cannot track an undefined variable");
				 }
			 TS->printResults(Simulation);
#ifndef _NO_GLUT_
			 if (Simulation->Plots->gl_on && Simulation->Plots->plot_num)  glutMainLoop();
#endif
			 delete Simulation;
		 }

		 vary.draw();
#ifndef _NO_GLUT_
		 if (vary.plots->gl_on && step == vary.steps - 1) {
			 GluPlotArray = vary.plots;
			 glutMainLoop();
		 }
#endif
		 delete TS;

	 } // ****************************** Loop over steps

 }
 catch (char *str) { if ( !equal(str,"") ) {
                        fprintf(stderr, "\n\n*** Error: %s\n", str);
			fflush(stderr); 
                        delete [] str;  }
                   }

 catch (int i) { fprintf(stderr, "\n\n*** CalC: breaking execution on exception #%d ***\n", i); 
                 perror(" System error (if any): ");
                 fflush(stderr); }

 if (VERBOSE > 3) {  fprintf(stderr, "\n\n*** Enter any string to exit CalC ***\n"); 
		     fflush(stderr);
             char s[256] = "0";
		     scanf("%s", s); }

 return 0;
 }


 //**********************************************************************************************
 //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 //**********************************************************************************************

//**************************************************************************
//                C A L C I U M   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) Ca*  = (1 + Ax/2 + Ay + Az) Ca + dt H(Ca,B)
//  (1 - Ay/2) Ca** = -Ay/2 Ca + Ca*
//  (1 - Az/2) Ca^  = -Az/2 Ca + Ca** + dt / 2 * { H(Ca^,B^) - H(Ca,B) }
//
//  H(Ca,B) = Sum_B{ - kplus * Ca * B + kminus * (B_total - B) }
//
//  Ca = Ca(t_n), Ca^ = Ca(t_{n+1}), B = B(t_{n+1/2})
//**************************************************************************

void Ca3DstepCoop(FieldObj &Ca, VectorObj &CaNew, BufferArray &Buf, BufferArray &BufNew, double dt)
{
	double    dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
	double    nu = dtHalf * Ca.getD();	
	int       BN = Buf.buf_num;
	int       NC = Buf.nonCoopNum;

	VectorObj LinOld(FieldObj::Size);
	VectorObj LinNew(FieldObj::Size);
	VectorObj Temp(FieldObj::Size);

	CaNew  =  Ca  + dt * Ca.bgr * *Ca.kuptake;
	LinOld = -dtHalf * *Ca.kuptake; 
	LinNew = LinOld; 
	Temp = 0.0;
	
	for (int b = 0; b < NC; b++) {
		LinOld -= (dtHalf * Buf.array[b]->kplus ->Evaluate() ) *  *Buf.array[b];
		LinNew -= (dtHalf * Buf.array[b]->kplus ->Evaluate() ) *  *BufNew.array[b];
		CaNew  +=     (dt * Buf.array[b]->kminus->Evaluate() ) * (*Buf.array[b]->total - *Buf.array[b]);
		Temp   -= (dtHalf * Buf.array[b]->kminus->Evaluate() ) * (*BufNew.array[b] - *Buf.array[b]); // Note the overall sign 
    }
	
	for (int b = NC; b < BN-1; b++) {
		if ( (b-NC) % 3 == 2 ) continue;  // skip double-bound buffer
		LinOld -= (dtHalf * Buf.array[b]->kplus ->Evaluate() ) *  *Buf.array[b];
		LinNew -= (dtHalf * Buf.array[b]->kplus ->Evaluate() ) *  *BufNew.array[b];
		CaNew  +=    ( dt * Buf.array[b]->kminus->Evaluate() ) *  *Buf.array[b+1];
		Temp   += (dtHalf * Buf.array[b]->kminus->Evaluate() ) * (*BufNew.array[b+1] - *Buf.array[b+1]);
    }
	
	CaNew += 2 * (Ca % LinOld);
	Ca.add_sources(CaNew, dt);
    Ca.Run3Dz(-nu,   2*nu,  2*nu,    nu,           0,  CaNew.elem);
	Ca.Run3Dx(-nu,    -nu,     0.,   0.,           0,  CaNew.elem);
    CaNew += Temp - LinOld % Ca;
    Ca.Run3Dy(-nu,     0.,   -nu,    0., LinNew.elem, CaNew.elem);
}

//**************************************************************************
//                  B U F F E R   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) B*  = (1 + Ax/2 + Ay + Az) B + dt H(Ca,B)
//  (1 - Ay/2) B** = -Ay/2 B + B*
//  (1 - Az/2) B^  = -Az/2 B + B** + dt / 2 * { H(Ca^,B^) - H(Ca,B) }
//
//  H(Ca,B) = B [-kplus * Ca - kminus] + kminus * B_total
//
//  Ca = Ca(t_n), Ca^ = Ca(t_{n+1})
//**************************************************************************

void Buf3DstepNew(BufferObj &Buf, VectorObj &BufNew, VectorObj &Ca, VectorObj &CaNew, double dt)
{
    double nu     = 0.5 * dt * Buf.getD();
	double kplus  = 0.5 * dt * Buf.kplus->Evaluate();
	double kminus = 0.5 * dt * Buf.kminus->Evaluate();

	VectorObj    LinNew( -kplus * CaNew - kminus );
	VectorObj    LinOld( -kplus * Ca    - kminus );
	BufNew = Buf % ( 1.0 + 2 * LinOld ) + (2 * kminus) * *Buf.total;
	
	Buf.Run3Dz(-nu, 2*nu,  2*nu,  nu, 0, BufNew.elem);
	Buf.Run3Dx(-nu, -nu,   0.,    0., 0, BufNew.elem);
	
	BufNew -= (LinOld % Buf);
    Buf.Run3Dy(-nu, 0.,   -nu,    0.,  LinNew.elem, BufNew.elem);
}

//**************************************************************************
//         H(Ca,B) = B [-kplus * Ca ] + kminus * B[b+1]
//**************************************************************************

void Buf3DstepCoop(BufferArray &Buf, BufferArray &BufNew, VectorObj &Ca, VectorObj &CaNew, double dt)
{
	double  dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
	int     BN = Buf.buf_num;
	int     NC = Buf.nonCoopNum;

	for (int bn = 0; bn < NC; bn++)
		Buf3DstepNew(*Buf.array[bn], *BufNew.array[bn], Ca, CaNew, dt);

	for (int UB = NC; UB < BN; UB += 3) {

		int SB = UB + 1; // Singly-bound buffer index
		int DB = SB + 1; // Doubly-bound buffer index

		double nuUB = dtHalf * Buf.array[UB]->getD();
		double nuSB = dtHalf * Buf.array[SB]->getD();
		double nuDB = dtHalf * Buf.array[DB]->getD();
		
		double kp1  = dtHalf * Buf.array[UB]->kplus->Evaluate();
		double km1  = dtHalf * Buf.array[UB]->kminus->Evaluate();
		double kp2  = dtHalf * Buf.array[SB]->kplus->Evaluate();
		double km2  = dtHalf * Buf.array[SB]->kminus->Evaluate();
		
		VectorObj    LinNew( - kp1 * CaNew ); // ********** Unbound buffer block **************
		VectorObj    LinOld( - kp1 * Ca    );

		*BufNew.array[UB] = *Buf.array[UB] % (1.0 + 2 * LinOld) + (2*km1) * (*Buf.array[SB]);
		
		Buf.array[UB]->Run3Dx( -nuUB, nuUB, 2*nuUB, 2*nuUB,          0, BufNew.array[UB]->elem);
		Buf.array[UB]->Run3Dy( -nuUB,   0.,  -nuUB,     0.,          0, BufNew.array[UB]->elem);

		*BufNew.array[UB] += km1 * (*BufNew.array[SB] - *Buf.array[SB]) - (LinOld % *Buf.array[UB]);

		Buf.array[UB]->Run3Dz( -nuUB,   0.,     0.,  -nuUB, LinNew.elem, BufNew.array[UB]->elem);

		LinNew = -kp2 * CaNew - km1;   // ************* Single-bound buffer block **************
		LinOld = -kp2 * Ca    - km1;

		*BufNew.array[SB] = *Buf.array[SB] % (1.0 + 2 * LinOld) + (2*km2) * (*Buf.array[DB]) + (2*kp1) * (Ca % *Buf.array[UB]);
		
		Buf.array[SB]->Run3Dx(-nuSB, nuSB, 2*nuSB, 2*nuSB,           0, BufNew.array[SB]->elem);
		Buf.array[SB]->Run3Dy(-nuSB,   0.,  -nuSB,     0.,           0, BufNew.array[SB]->elem);

		*BufNew.array[SB] += km2 * (*BufNew.array[DB] - *Buf.array[DB])
						   + kp1 * (CaNew % *BufNew.array[UB] - Ca % *Buf.array[UB]) - (LinOld % *Buf.array[SB]);

		Buf.array[SB]->Run3Dz(-nuSB,   0.,     0.,  -nuSB,  LinNew.elem, BufNew.array[SB]->elem);

		LinOld = LinNew = -km2;           // ************* Double-bound buffer block **************
		
		*BufNew.array[DB] = *Buf.array[DB] % (1.0 + 2 * LinOld) + (2*kp2) * Ca % (*Buf.array[SB]);
		
		Buf.array[DB]->Run3Dx(-nuDB, nuDB, 2*nuDB, 2*nuDB,           0, BufNew.array[DB]->elem);
		Buf.array[DB]->Run3Dy(-nuDB,   0.,  -nuDB,     0.,           0, BufNew.array[DB]->elem);

		*BufNew.array[DB] += kp2 * (CaNew % *BufNew.array[SB] - Ca % *Buf.array[SB]) - (LinOld % *Buf.array[DB]);

		Buf.array[DB]->Run3Dz(-nuDB,   0.,     0.,  -nuDB, LinNew.elem, BufNew.array[DB]->elem);
	}
}


//**************************************************************************
//                C A L C I U M   2 D   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) Ca*  = (1 + Ay/2) Ca  + dt/2 H(Ca,B)
//  (1 - Ay/2) Ca^  = (1 + Ax/2) Ca* + dt/2 H(Ca^,B^)
//
//  H(Ca,B) = Sum_B{ - kplus * Ca * B + kminus * (B_total - B) }
//
//  Ca = Ca(t_n), Ca^ = Ca(t_{n+1}), B = B(t_{n+1/2})
//**************************************************************************

void Ca2DstepCoop(FieldObj &Ca, VectorObj &CaNew, BufferArray &Buf, BufferArray &BufNew, double dt)
{
	double    dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
	int       BN = Buf.buf_num;
	int       NC = Buf.nonCoopNum;
	VectorObj LinOld(FieldObj::Size);
	VectorObj LinNew(FieldObj::Size);
	FieldObj  Temp(Ca);
	double    nu = Ca.getD() * dtHalf;
	
	CaNew  = Ca  + dtHalf * Ca.bgr * *Ca.kuptake;
	LinOld = -dtHalf * *Ca.kuptake; 
	LinNew = LinOld; 
	
	for (int b = 0; b < NC; b++) {
		LinOld -= (dtHalf *  Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];   
		LinNew -= (dtHalf *  Buf.array[b]->kplus->Evaluate()) * *BufNew.array[b]; 
		CaNew  +=  dtHalf *(*Buf.array[b]->total - *Buf.array[b]) * Buf.array[b]->kminus->Evaluate();
    }

	for (int b = NC; b < BN-1; b++) {
		if ( (b-NC) % 3 == 2 ) continue;  // skip double-bound buffer
		LinOld -= (dtHalf *  Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];   
		LinNew -= (dtHalf *  Buf.array[b]->kplus->Evaluate()) * *BufNew.array[b]; 
		CaNew  +=  dtHalf *(*Buf.array[b+1]) * Buf.array[b]->kminus->Evaluate();
    }
	
	CaNew += Ca % LinOld;
	Ca.add_sources(CaNew, dtHalf);
    
	Ca.Run2Dx(-nu,   0.0,   nu,     0,  CaNew.elem);
	Temp = CaNew;   // IMPORTANT! cf. 3D methods: Ca* in second split operator, not Ca
    CaNew += dtHalf * Ca.bgr * *Ca.kuptake ; 

	for (int b = 0;  b < NC;   b++) CaNew += dtHalf * (*BufNew.array[b]->total - *BufNew.array[b]) * BufNew.array[b]->kminus->Evaluate();
    for (int b = NC; b < BN-1; b++) 
	  	  if ( (b-NC) % 3 != 2 )    CaNew += dtHalf * (*BufNew.array[b+1]) * Buf.array[b]->kminus->Evaluate();
	
    Ca.add_sources(CaNew, dtHalf);
    Temp.Run2Dy(-nu,    nu,     0.0,     LinNew.elem,   CaNew.elem);
}

//**************************************************************************

//**************************************************************************
//               B U F F E R   2 D   D - G   S T E P
//**************************************************************************
//
//  (1 - Ax/2) B*  = (1 + Ay/2) B  + dt/2 H(Ca,B)
//  (1 - Ay/2) B^  = (1 + Ax/2) B* + dt/2 H(Ca^,B^)
//
//  H(Ca,B) = - kplus * Ca * B + kminus * (B_total - B)
//
//  B = B(t_n), B^ = B(t_{n+1}), Ca = Ca(t_{n+1/2})
//**************************************************************************

void Buf2DstepNew(BufferObj &Buf, VectorObj &BufNew, VectorObj &Ca, VectorObj &CaNew, double dt)
{ 
  double   dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
  FieldObj Temp(Buf);
  double   nu = Buf.getD() * dtHalf;

  double kplus  = dtHalf * Buf.kplus->Evaluate();
  double kminus = dtHalf * Buf.kminus->Evaluate();

  VectorObj  LinNew( - (kplus * CaNew + kminus) );
  BufNew =  Buf % (1 - (kplus * Ca + kminus)) + kminus * *Buf.total;

  Buf.Run2Dx (-nu,     0.0,    nu,  0,  BufNew.elem);

  Temp = BufNew;
  BufNew += kminus * *Buf.total; // IMPORTANT! cf. 3D methods.
  Temp.Run2Dy(-nu, nu,     0.0,     LinNew.elem,  BufNew.elem);
}
//**************************************************************************

void Buf2DstepCoop(BufferArray &Buf, BufferArray &BufNew, VectorObj &Ca, VectorObj &CaNew, double dt)
{
	int BN = Buf.buf_num;
	int NC = Buf.nonCoopNum;
    double dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
	if (BN == 0) return;

	FieldObj  Temp(*Buf.array[0]);
 
	for (int bn = 0; bn < NC; bn++)
		Buf2DstepNew(*Buf.array[bn], *BufNew.array[bn], Ca, CaNew, dt);

	for (int UB = NC; UB < BN; UB += 3) {

		int SB = UB + 1; // Singly-bound buffer index
		int DB = SB + 1; // Doubly-bound buffer index
	
		double nuUB = Buf.array[UB]->getD() * dtHalf;
		double nuSB = Buf.array[SB]->getD() * dtHalf;
		double nuDB = Buf.array[DB]->getD() * dtHalf;
		
		double kp1 = dtHalf * Buf.array[UB]->kplus->Evaluate();
		double km1 = dtHalf * Buf.array[UB]->kminus->Evaluate();
		double kp2 = dtHalf * Buf.array[SB]->kplus->Evaluate();
		double km2 = dtHalf * Buf.array[SB]->kminus->Evaluate();
		
		// ********************************* Unbound Buffer:

		*BufNew.array[UB] = *Buf.array[UB] % (1.0 - kp1 * Ca) + km1 * (*Buf.array[SB]);

		Buf.array[UB]->Run2Dx( -nuUB, 0.0, nuUB, 0, BufNew.array[UB]->elem);

		Temp = *BufNew.array[UB];
		VectorObj    LinNew( -kp1 * CaNew );    
		*BufNew.array[UB] +=  km1 * *BufNew.array[SB];

		Temp.Run2Dy(-nuUB, nuUB, 0.0, LinNew.elem,  BufNew.array[UB]->elem);

		// ********************************* Singly-Bound Buffer:

		*BufNew.array[SB] = *Buf.array[SB] % (1.0 - km1 - kp2*Ca) + (km2 * *Buf.array[DB] + kp1 * Ca % *Buf.array[UB]);

		Buf.array[SB]->Run2Dx( -nuSB, 0.0, nuSB, 0, BufNew.array[SB]->elem);
		Temp = *BufNew.array[SB];
		LinNew = -km1 - kp2*CaNew;    
		*BufNew.array[SB] += km2 * *BufNew.array[DB] + kp1 * CaNew % *BufNew.array[UB];

		Temp.Run2Dy(-nuSB, nuSB, 0.0, LinNew.elem, BufNew.array[SB]->elem);

		// ********************************* Double-Bound Buffer:

		*BufNew.array[DB] = *Buf.array[DB] * (1.0 - km2) + kp2 * Ca % *Buf.array[SB];

		Buf.array[DB]->Run2Dx( -nuDB,  0.0,   nuDB, 0,   BufNew.array[DB]->elem);
		Temp = *BufNew.array[DB];
		LinNew = -km2;
	
		*BufNew.array[DB] += kp2 * CaNew % *BufNew.array[SB];

		Temp.Run2Dy(-nuDB, nuDB, 0.0, LinNew.elem, BufNew.array[DB]->elem);
	}
}

/**************************************************************************
     C A L C I U M   C R A N K - N I C H O L S O N   1 D   S T E P
***************************************************************************

  (1 - Ax/2) Ca^ = (1 + Ax/2) Ca + dt/2 ( H(Ca, B) + H(Ca^,B^) ) + dt sources(t_{n+1/2})

  (1 - Ax/2 - LinNew ) Ca^ = Ax/2 Ca + (1 + LinOld ) Ca
                         + dt { Sum_B{ kminus * (B_total - <B>) + sources(t_{n+1/2}) }

  Lin = - dt/2 Sum_B{ kplus * B }
****************************************************************************
  (1 - Ax/2 - LinNew ) Ca^ = Ax/2 Ca + (1 + LinOld ) Ca
                         + dt { Sum_B{ kminus * <B*> + sources(t_{n+1/2}) }
****************************************************************************/

void Ca1DstepCoop(FieldObj &Ca, VectorObj &CaNew, BufferArray &Buf, BufferArray &BufNew, double dt)
{
	double    dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
	//double    tStore = Ca.Time;
    int       BN = Buf.buf_num;
	int       NC = Buf.nonCoopNum;   
	VectorObj LinOld(FieldObj::Size);
	VectorObj LinNew(FieldObj::Size);
	double    nu = Ca.getD() * dtHalf;

	CaNew  = Ca  + dt * Ca.bgr * *Ca.kuptake;
	LinOld =  -dtHalf * *Ca.kuptake; 
	LinNew = LinOld; 
	
	for (int b = 0; b < NC; b++) {
		LinOld -= (dtHalf * Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];
		LinNew -= (dtHalf * Buf.array[b]->kplus->Evaluate()) * *BufNew.array[b];
		CaNew  +=  dtHalf * (2 * *Buf.array[b]->total - *Buf.array[b] - *BufNew.array[b]) * Buf.array[b]->kminus->Evaluate();
    }

	for (int b = NC; b < BN-1; b++) {
		if ( (b-NC) % 3 == 2 ) continue;  // skip double-bound buffer
		LinOld -= (dtHalf * Buf.array[b]->kplus->Evaluate()) * *Buf.array[b];
		LinNew -= (dtHalf * Buf.array[b]->kplus->Evaluate()) * *BufNew.array[b];
		CaNew  +=  dtHalf * (*Buf.array[b+1] + *BufNew.array[b+1]) * Buf.array[b]->kminus->Evaluate();
    }
	
	CaNew += Ca % LinOld;
	Ca.add_sources(CaNew, dt);
    Ca.Run1D(-nu, nu, LinNew.elem, CaNew.elem);
}

/**************************************************************************
       B U F F E R    C R A N K - N I C H O L S O N   1 D   S T E P
**************************************************************************

  (1 - Ax/2) B^ = (1 + Ax/2) B + dt/2 ( H(Ca, B) + H(Ca^, B^) )

  ( 1 - Ax/2 - LinNew ) B^ = Ax/2 B + (1 + LinOld ) B + dt kminus B_total

  Lin = - dt/2 (kplus Ca + kminus)

  H(Ca,B) = - kplus * Ca * B + kminus * (B_total - B)

  B = B(t_n), B^ = B(t_{n+1}), Ca = Ca(t_{n+1/2})
****************************************************************************/
/**************************************************************************
  ( 1 - Ax/2 - Lin ) B^ = Ax/2 B + (1 + Lin ) B + dt kminus B_total

  Lin = - dt/2 (kplus Ca + kminus)
****************************************************************************/

void Buf1DstepNew(BufferObj &Buf, VectorObj &BufNew, VectorObj &Ca, VectorObj &CaNew, double dt)
{
	double dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
    double nu     = dtHalf * Buf.getD();
	double kplus  = dtHalf * Buf.kplus->Evaluate();
	double kminus = dtHalf * Buf.kminus->Evaluate();

	VectorObj    LinNew( -kplus * CaNew - kminus );
	VectorObj    LinOld( -kplus * Ca    - kminus );
	BufNew = Buf % ( 1.0 + LinOld ) + (2 * kminus) * *Buf.total;
	
    Buf.Run1D( -nu, nu, LinNew.elem, BufNew.elem);
}

//**************************************************************************

void Buf1DstepCoop(BufferArray &Buf, BufferArray &BufNew, VectorObj &Ca, VectorObj &CaNew, double dt)
{ 
	double dtHalf = 0.5 * dt;  // NOTE FACTOR OF 1/2 IN ALL FLUXES!
    int BN = Buf.buf_num;
	int NC = Buf.nonCoopNum;

	for (int bn = 0; bn < NC; bn++)
		Buf1DstepNew(*Buf.array[bn], *BufNew.array[bn], Ca, CaNew, dt);

	for (int UB = NC; UB < BN; UB += 3) {

		int SB = UB + 1; // Singly-bound buffer index
		int DB = SB + 1; // Doubly-bound buffer index

		double nuUB = Buf.array[UB]->getD() * dtHalf, nuSB = Buf.array[SB]->getD() * dtHalf, nuDB = Buf.array[DB]->getD() * dtHalf;
		
		double kp1 = dtHalf * Buf.array[UB]->kplus->Evaluate(), km1 = dtHalf * Buf.array[UB]->kminus->Evaluate();
		double kp2 = dtHalf * Buf.array[SB]->kplus->Evaluate(), km2 = dtHalf * Buf.array[SB]->kminus->Evaluate();
		
		VectorObj    LinNew( - kp1 * CaNew );  // ********* Unbound buffer block ******************
		VectorObj    LinOld( - kp1 * Ca    );

		*BufNew.array[UB] = *Buf.array[UB] % (1.0 + LinOld) + km1 * (*Buf.array[SB] + *BufNew.array[SB]);
		Buf.array[UB]->Run1D(-nuUB, nuUB,  LinNew.elem, BufNew.array[UB]->elem);

		LinNew = - kp2 * CaNew - km1 ;      // ********* Single-bound buffer block **************
		LinOld = - kp2 * Ca    - km1 ;

		*BufNew.array[SB] = *Buf.array[SB] % (1.0 + LinOld) + km2 * (*Buf.array[DB] + *BufNew.array[DB])
						  + (0.5 * kp1) * (Ca + CaNew) % (*Buf.array[UB] + *BufNew.array[UB]);
		
		Buf.array[SB]->Run1D(-nuSB, nuSB,  LinNew.elem, BufNew.array[SB]->elem);

		LinOld = LinNew = -km2;              // ********* Double-bound buffer block **************
		
		*BufNew.array[DB] = *Buf.array[DB] % (1.0 + LinOld) + (0.5 * kp2) * (Ca + CaNew) % (*Buf.array[SB] + *BufNew.array[SB]);
		Buf.array[DB]->Run1D(-nuDB, nuDB,  LinNew.elem, BufNew.array[DB]->elem);
	}
}
//**************************************************************************

const char *getMethod( CaMethod &CaStep, BufMethod &BufStep) {

 switch (DIMENSIONALITY)  { 

 case 1: CaStep = &Ca1DstepCoop; BufStep = &Buf1DstepCoop; 
	     return makeMessage("1D (%s) Crank-Nicholson implicit", LABEL_DIM1);
         break;

 case 2: CaStep = &Ca2DstepCoop; BufStep = &Buf2DstepCoop; 
	     return  makeMessage("2D (%s,%s) Crank-Nicholson ADI", LABEL_DIM1, LABEL_DIM2);
	     break;

 case 3: CaStep = &Ca3DstepCoop; BufStep = &Buf3DstepCoop;       
	     return  makeMessage("3D (%s,%s,%s) Douglas-Gunn ADI", LABEL_DIM1, LABEL_DIM2, LABEL_DIM3);
		 break;
		 
 default: throw makeMessage("Cannot interpret Dimensionality=%d", DIMENSIONALITY);
		  return 0;		  

 }
 }

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


void getRun( TokenString &TS, int i, bool *adaptive, double *time, double *accuracy,
			 double *dtMax, double *dt, double *dtStretch, double *ODEaccuracy)
{
  char firstArg[MAX_TOKEN_LENGTH];

  if (adaptive) *adaptive = false;             // non-adaptive is the default
  long pos = TS.token_index("Run", i);

  TS.trail_pars(pos, 's', firstArg, 'E');
  pos ++;
  int pmax = TS.tokens_to_eol(pos);

  if ( pmax > 0 && !TS.isConst(firstArg) ) { 
    pos++; 
    if (adaptive) {
		if ( equal(firstArg,"adaptive") ) *adaptive = true;
	    else if ( !equal(firstArg, "nonadaptive") ) TS.errorMessage(pos, 0, "Undefined run method");
	}
  }

  if (time) {
    *time = TS.get_double(pos);
    if (*time <= 0) TS.errorMessage(pos, makeMessage("Simulation time has to be positive (T = %lf)", *time) );
  } 

  if ( !isLineEnd(TS[pos + 1]) )  {
	if (*adaptive) 
		TS.trail_pars(pos, 'd', accuracy, 'd', dtMax, 'd', dt, 'd', dtStretch, 'd', ODEaccuracy, 'E');
	else 
		TS.trail_pars(pos, 'd', dtMax, 'd', accuracy, 'E');
  }
  if (dt) if (*dt <= 0) 
    TS.errorMessage(pos+1, makeMessage("Time step must be positive (dt = %lf)", *dt) );
  if (accuracy) if (*accuracy <= 0 || *accuracy > 1.0) // this argument is actually ODEaccuracy for non-adaptive run
    TS.errorMessage(pos+2, makeMessage("Accuracy must be between 0 and 1 (accuracy = %lf)", *accuracy) );
  if (dtMax) if (*dtMax <= 0 ) 
    TS.errorMessage(pos+1, makeMessage("Max time step must be positive (dtMax = %lf)", *dtMax) );
  if (dtStretch) if (*dtStretch <= 1 || *dtStretch >= 2) 
    TS.errorMessage(pos+4, makeMessage("Time step stretch must be between 1 and 2 (dtStretch = %lf)", *dtStretch) );
  if (ODEaccuracy) if (*ODEaccuracy <= 0 || *ODEaccuracy > 1.0) 
    TS.errorMessage(pos+5, makeMessage("ODE accuracy must be between 0 and 1 (ODEaccuracy = %lf)", *ODEaccuracy) );
}


//**************************************************************************

double get_sim_time(TokenString &params)
  {
  double tt, t = 0.0;
  bool adaptive;

  for (int i = 1; i<= params.token_count("Run"); i++) {
    getRun(params, i, &adaptive, &tt); 
    t += tt;
  }
  return t;
  }

//**************************************************************************


void setFieldByFunction(TokenString &TS, long p, double *Field, const char *errStr) {

  int xsize = FieldObj::Grid->xsize;  double *xcoord = FieldObj::Grid->xcoord;
  int ysize = FieldObj::Grid->ysize;  double *ycoord = FieldObj::Grid->ycoord;
  int zsize = FieldObj::Grid->zsize;  double *zcoord = FieldObj::Grid->zcoord;

  int  ix, iy, iz;
  long l = 0;
  double x, y, z;
  ExpressionObj *T = 0;

  switch (DIMENSIONALITY) {

  case 1:
    T = new ExpressionObj(TS, p, makeMessage("%s (function of %s)", errStr, LABEL_DIM1), 0, 0, &x, LABEL_DIM1); 
    for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate();
    }
    break;  

  case 2:
    T = new ExpressionObj(TS, p, 
            makeMessage("%s (function of %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2), 0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2); 
    for (iy = 0; iy < ysize; iy++) {
        y = ycoord[iy];
        for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate();
	}
    }
    break;  
 
  case 3:
    T = new ExpressionObj(TS, p, 
         makeMessage("%s (function of %s, %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2, LABEL_DIM3), 
         0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2, &z, LABEL_DIM3); 
    for (iz = 0; iz < zsize; iz++) {
      z = zcoord[iz];
      for (iy = 0; iy < ysize; iy++) {
        y = ycoord[iy];
        for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate();
	}
      }
    }
    break;  
 }

 if (T && VERBOSE) T->print(stderr); 
 if (T) delete T;
 return;
}

//**************************************************************************


void setFieldByFunction(TokenString &TS, long p, bool *Field, const char *errStr) {

  int xsize = FieldObj::Grid->xsize;  double *xcoord = FieldObj::Grid->xcoord;
  int ysize = FieldObj::Grid->ysize;  double *ycoord = FieldObj::Grid->ycoord;
  int zsize = FieldObj::Grid->zsize;  double *zcoord = FieldObj::Grid->zcoord;

  int  ix, iy, iz;
  long l = 0;
  double x, y, z;
  ExpressionObj *T = 0;

  switch (DIMENSIONALITY) {

  case 1:
    T = new ExpressionObj(TS, p, makeMessage("%s (function of %s)", errStr, LABEL_DIM1), 0, 0, &x, LABEL_DIM1); 
    for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate() > 0;
    }
    break;  

  case 2:
    T = new ExpressionObj(TS, p, 
            makeMessage("%s (function of %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2), 0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2); 
    for (iy = 0; iy < ysize; iy++) {
        y = ycoord[iy];
        for (ix = 0; ix < xsize; ix++)  {
          x = xcoord[ix];
          Field[l++] = T->Evaluate() > 0;
	}
    }
    break;  
 
  case 3:
    T = new ExpressionObj(TS, p, 
         makeMessage("%s (function of %s, %s, %s)", errStr, LABEL_DIM1, LABEL_DIM2, LABEL_DIM3), 
         0, 0, &x, LABEL_DIM1, &y, LABEL_DIM2, &z, LABEL_DIM3); 

	for (iz = 0; iz < zsize; iz++) 
	{
		z = zcoord[iz]; 
		for (iy = 0; iy < ysize; iy++) 
		{
			y = ycoord[iy]; 
			for (ix = 0; ix < xsize; ix++)  
			{
				x = xcoord[ix];
				double val = T->Evaluate();
				Field[l++] = (val > 0);
			}
		}
	}
	break;  
 }

  if (T && VERBOSE) { T->print(stderr); fprintf(stderr, "\n"); }
 if (T) delete T;
 return;
}

//**************************************************************************
