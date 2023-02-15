/**************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *                Copyright (C) 2001-2021 Victor Matveev
 *
 *                               field.cpp
 *
 *  Definitions of the main FieldObj concentration field class (derived from
 *  the VectorObj class), the daughter BufferObj buffer concentration field 
 *  class, and the compound buffer array class BufferArray.
 *
 *  FieldObj methods include all low-level ADI difference schemes routines.
 *  Higher-level numeric engine routines are defined in "calc.cpp"       
 *       
 *  FieldObj and BufferArray variables are member of the main SimulationObj
 *  object, defined in simulation.h / simulation.cpp
 *
 *************************************************************************
 
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
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include "vector.h"
#include "syntax.h"
#include "box.h"
#include "grid.h"
#include "field.h"

extern void setFieldByFunction(TokenString &, long, double *, const char *errStr);                                     

const double pi = atan(1.0) * 4;

double (*integrate)(double, double, double, int) = &gaussian;

//****************************************************************************
//                      INITIALIZING STATIC MEMBERS:
//****************************************************************************

class RegionObj      *FieldObj::Region  = 0;
class GridObj        *FieldObj::Grid    = 0;
class BCarrayObj     *FieldObj::BCarray = 0;

char   FieldObj::SAME_CURRENT = 1;

int    FieldObj::xsize   = 0,    FieldObj::ysize   = 0,  FieldObj::zsize  = 0;
long   FieldObj::Size    = 0,    FieldObj::xysize  = 0;

double *FieldObj::xgrid  = 0,   *FieldObj::ygrid   = 0, *FieldObj::zgrid  = 0;
double *FieldObj::xcoord = 0,   *FieldObj::ycoord  = 0, *FieldObj::zcoord = 0;

double *FieldObj::dxplus = 0,   *FieldObj::dxminus = 0, *FieldObj::dvx = 0;
double *FieldObj::dyplus = 0,   *FieldObj::dyminus = 0, *FieldObj::dvy = 0;
double *FieldObj::dzplus = 0,   *FieldObj::dzminus = 0, *FieldObj::dvz = 0;
  
double *FieldObj::bc_deriv = 0, *FieldObj::bc_coef = 0, *FieldObj::bc_pump  = 0, *FieldObj::bc_lin = 0;
double *FieldObj::bc_Kn    = 0, *FieldObj::bc_pow  = 0, *FieldObj::bc_const = 0;

char  **FieldObj::bc_id       = 0;
int     FieldObj::bc_type_num = 0;

double *FieldObj::right  = 0, *FieldObj::diag   = 0, *FieldObj::sup    = 0,  *FieldObj::sub = 0;

struct TermStruct *Currents   = 0;
class  ExpressionObj *Current = 0;

//*************************************************************************
//*                T R I D I A G O N A L   S Y S T E M
//*************************************************************************

void FieldObj::setStaticData(RegionObj &RO, GridObj &GO, BCarrayObj &BCO) {

  Region = &RO;

  Grid    = &GO;
  Size    = GO.Size;     xysize  = GO.xysize; 
  xsize   = GO.xsize;    ysize   = GO.ysize;    zsize   = GO.zsize;

  xgrid   = GO.xgrid;    ygrid   = GO.ygrid;    zgrid   = GO.zgrid; 
  xcoord  = GO.xcoord;   ycoord  = GO.ycoord;   zcoord  = GO.zcoord; 

  dxplus  = GO.dxplus;   dyplus  = GO.dyplus;   dzplus  = GO.dzplus; 
  dxminus = GO.dxminus;  dyminus = GO.dyminus;  dzminus = GO.dzminus;
  dvx     = GO.dvx;      dvy     = GO.dvy;      dvz     = GO.dvz;  
  
  BCarray  = &BCO;          bc_type_num = BCO.bc_type_num;  
  bc_id    = BCO.bc_id;     bc_const    = BCO.bc_const;
  bc_deriv = BCO.bc_deriv;  bc_coef     = BCO.bc_coef;  
  bc_lin   = BCO.bc_lin;    bc_pump     = BCO.bc_pump;
  bc_Kn    = BCO.bc_Kn;     bc_pow      = BCO.bc_pow;  

  int max = ( xsize > ysize ? xsize : ysize );
  if (zsize > max) max = zsize;
  init_tridiag(max + 2);
} 


void FieldObj::init_tridiag(int n)
{
diag = new double[n];  right = new double[n];
sup  = new double[n];  sub   = new double[n];
}

void FieldObj::kill_tridiag()
{
  if (diag) {
    delete [] diag; delete [] right;
    delete [] sup;  delete [] sub;
  }
}

//*******************   S O L V E R :  **************************

void FieldObj::tridiag(int n)
{
int    i;
double sb, sp, dg, rt;

sp = ( sup[1]   /= (dg = diag[1]) );
rt = ( right[1] /= dg );

for (i=2; i<=n; i++)
  {
  sb = sub[i];
  rt = right[i] = (right[i] - sb * rt) / (dg = diag[i] - sb * sp);
  sp = (sup[i] /= dg);
  }

for (i=n-1; i>=1; i--) right[i] -= sup[i] * right[i+1];
}

//***************************************************************************

void FieldObj::getCurrents(TokenString &Params, int simID, class VarList *VL, double *tptr) 
{
  long pos;
  double res, *ptr;
  int  si, c0, c1;
  
  if (equal(ID,"Ca")) {
     c0 = Params.token_count("current" );
     c1 = Params.token_count("currents");
  } else {
     c0 = Params.token_count(ID, "current" );
     c1 = Params.token_count(ID, "currents");
  } 

  if (source_num == 0) return;
  if (c0 && c1) globalError(makeMessage("May not mix \"current\" and \"currents\" commands in the same file"));

  FieldObj::SAME_CURRENT = (c0 > 0);

  if ( FieldObj::SAME_CURRENT ) {

     if (equal(ID,"Ca")) pos = Params.token_index   ( "current", simID ) + 1;
	                else pos = Params.token_index(ID, "current", simID ) + 1;
     if ( Params.equal(pos, VAR_TOKEN) || Params.equal(pos, ASSIGN_TOKEN) ) pos++;
     Current = new ExpressionObj(Params, pos, "Bad expression for the current", VL, 0, tptr, "t");
     Current->Evaluate();

  } else {
 
	 Current = 0;
     if (equal(ID,"Ca")) pos = Params.token_index   ( "currents", simID ) + 1;
	                else pos = Params.token_index(ID, "currents", simID ) + 1;
     if ( Params.equal(pos, VAR_TOKEN) || Params.equal(pos, ASSIGN_TOKEN) ) pos++;

     for (si = 0; si < source_num; si++) 
       if ( Params.isConst(Params[pos + si], &res) ) { 
         Currents[si].type = NUMBER_TYPE;
         Currents[si].val  = res; 
       }
       else if ( ( ptr = VL->ResolveID(Params[pos + si]) ) ) {
         Currents[si].type = POINTER_TYPE;
         Currents[si].ptr  = ptr;
       }
       else Params.errorMessage(pos + si, makeMessage("Bad specifier for current #%d: expecting a number or time-dependent variable",si+1));
    }  
}


//***************************************************************************

void FieldObj::evaluateCurrents() 
{
  if (!source_num) return;
  ICa = 0;

  for (int i=0; i<source_num; i++)
    ICa += FieldObj::SAME_CURRENT ? Current->Evaluate() : Evaluate(Currents[i]);
}

//***************************************************************************

void FieldObj::printCurrents(FILE *out) 
{
  if (!source_num) return;

  if (FieldObj::SAME_CURRENT && Current )  fprintf(out, "%g ", Current->Evaluate() );
  else {
    fprintf(out,"{ ");
    for (int i=0; i<source_num; i++) if (Currents) fprintf(out, "%g ", Evaluate(Currents[i]) );
    fprintf(out,"}");
  }
  fprintf(out,"\n");
}

//***************************************************************************

void FieldObj::killCurrents() 
{
  if (!source_num) return;

  if (FieldObj::SAME_CURRENT) { if (Current) delete Current; Current = 0; }
}

//****************************************************************************
//*          C L A S S   F I E L D :   I M P L E M E N T A T I O N
//****************************************************************************


FieldObj::FieldObj(TokenString &params, const char *id) : VectorObj(Size)
{
  int i;

  integrate  = params.Assert("current.shape","square") ? &square : &gaussian;
  fieldObstNum = params.token_count(id, "obstacle");
  source_num = params.token_count(id, "source");

  //params.checkName(id, "concentration field");
  allocate(id);

  kuptake  = 0; // Bound in simulation constructor;
  Time     = 0.0;  
  *this    = bgr = 0.0;
  D        = 0.2;

  params.get2_param(id, "D", &D);   
  if (VERBOSE)
    fprintf(stderr, "\n### Concentration field %s: diffusion coefficient = %g um^2 / ms\n", ID, D);

  tortDefined    = false;
  long tortIndex = 0;                             // this blocks sets tortuosity
  Diff           = 0;
  params.token_count(id, "tortuosity", &tortIndex);

  if ( tortIndex ) {
    tortDefined = true;
    Diff = new double[ Size ];
    if (VERBOSE) fprintf(stderr,"\n#### Tortuosity function = ");
    setFieldByFunction( params, tortIndex + 1, Diff, "Bad tortuosity expression" ); 
    if (VERBOSE) fprintf(stderr,"\n");
  }

  if (params.token_count(id, "bgr")) {
    bgr = params.get2_param(id, "bgr");
    if (bgr < 0) 
       globalError( makeMessage("Background %s concentration cannot be negative ([%s]bgr = %g)", ID, ID, bgr) );
    *this = bgr;
    if (VERBOSE)  fprintf(stderr, "    Background [%s] = %g uM \n", ID, bgr);
  }

  long pos;
  if (params.token_count(id, "import", &pos))
    { 
    char *filename = params.get_string(pos + 1);
    if (VERBOSE) fprintf(stderr, "  # Importing %s from file %s\n", ID, filename);
    import(filename);
    delete filename;
    }

//******************************************************************
  
  fieldObstArray = fieldObstNum ? new VolumeObjClass[fieldObstNum] : 0;

  if (fieldObstNum) 
	  for (i = 0; i < fieldObstNum; i++)  
	  {
      pos = params.token_index(ID, "obstacle", i + 1);
	  fieldObstArray[i] = VolumeObjClass(params, pos, 1);
	  }
 
  //****************************************************************

  int defs  = params.token_count(id, "bc");
  int boxes = Region->get_box_num() + fieldObstNum;

  char bc1[80], bc2[80], bc3[80], bc4[80], bc5[80], bc6[80];

  for (i = 1; i<= boxes; i++)
    try { 

    if ( i > defs )   // no b.c. specification for box i:
      {
      if (VERBOSE) fprintf(stderr,"    Box #%d: not specified; assuming all boundaries are zero flux\n", i);
      set_bc(i-1, "Noflux");     
      }
    else if ( params.equal( params.token_index(id, "bc", i) + 1, "all") ) 
      {
      params.trail_pars(id, "bc", i, 0,0, 's', bc1, 'E');
      if (VERBOSE) fprintf(stderr,"    Box #%d: all boundaries are %s\n", i, bc1);
      set_bc(i-1, bc1);
      }
    else if (DIMENSIONALITY == 1)
      {
      params.trail_pars(id, "bc", i, 's', bc1, 's', bc2, 'E');
      if (VERBOSE) 
        fprintf(stderr,"    Box #%d: boundary cond'ns: xmin:%s xmax:%s\n", i, bc1, bc2);
      set_bc(i-1, bc1, bc2, bc2, bc2, bc2, bc2);
      }
    else if (DIMENSIONALITY == 2)
      {
      params.trail_pars(id, "bc", i, 's', bc1, 's', bc2, 's', bc3, 's', bc4,  'E');
      if (VERBOSE) 
        fprintf(stderr,"    Box #%d: boundary cond'ns: xmin:%s xmax:%s ymin:%s ymax:%s\n", i, bc1, bc2, bc3, bc4);
      set_bc(i-1, bc1, bc2, bc3, bc4, bc4, bc4);
      }
    else  // DIMENSIONALITY == 3
      {
      params.trail_pars(id, "bc", i,
          's', bc1, 's', bc2, 's', bc3, 's', bc4, 's', bc5, 's', bc6, 'E');
      if (VERBOSE) 
        fprintf(stderr,"    Box #%d: boundary cond'ns: xmin:%s xmax:%s ymin:%s ymax:%s zmin:%s zmax:%s\n", 
                i, bc1, bc2, bc3, bc4, bc5, bc6);
      set_bc(i-1, bc1, bc2, bc3, bc4, bc5, bc6);
      }
    }
  catch(char *str) { params.errorMessage( params.token_index(id, "bc", i) + 1, str); }

  init_boundary();

  //******************************************************************

 if (source_num)
  {
  double x, y, z, sigmax, sigmay, sigmaz;
  if (VERBOSE) fprintf(stderr, "    %s has %d sources:\n", ID, source_num);
  for (i = 1; i <= source_num; i++)
	  try {
		x = y = z = 0;
		sigmax = 0.0; sigmay = sigmaz = -1;

		switch (DIMENSIONALITY) {
		case 1:
		   params.trail_pars(id, "source", i, 'd', &x, 'd', &sigmax);
		   break;
		case 2: 
		   params.trail_pars(id, "source", i, 'd', &x, 'd', &y, 'd', &sigmax, 'd', &sigmay);
		   break;
		case 3:
		   params.trail_pars(id, "source", i, 'd', &x, 'd', &y, 'd', &z, 'd', &sigmax, 
							  'd', &sigmay,  'd', &sigmaz);
		   break;
		}

		if (sigmay < 0) sigmay = sigmax;  //  by default widths in all directions are equal
		if (sigmaz < 0) sigmaz = sigmay;  //

		set_source(i-1, x, y, z, sigmax, sigmay, sigmaz); 
		} 
		catch (char *str) { 
		  params.errorMessage( params.token_index(id, "source", i) + 1, str, "Cannot initialize this source"); 
		}
  }

}

//***************************************************************************

FieldObj::FieldObj(const FieldObj &f) : VectorObj(f)
{
  int nx, ny, nz;
  source_num = f.source_num;
  fieldObstNum = f.fieldObstNum;
  allocate( f.ID );
 
  if (source_num)
    {
    Current = f.Current;
    for(int i = 0; i < source_num; i++)  
      {
      Currents[i]       = f.Currents[i];
      nx = source_nx[i] = f.source_nx[i];
      ny = source_ny[i] = f.source_ny[i];
      nz = source_nz[i] = f.source_nz[i];
      source_x0[i]      = f.source_x0[i];
      source_y0[i]      = f.source_y0[i];
      source_z0[i]      = f.source_z0[i];

      source_xsum[i] = new double[nx];
      source_ysum[i] = new double[ny];
      source_zsum[i] = new double[nz];

      int ii;
      for (ii = 0; ii < nx; ii++) source_xsum[i][ii] = f.source_xsum[i][ii];
      for (ii = 0; ii < ny; ii++) source_ysum[i][ii] = f.source_ysum[i][ii];
      for (ii = 0; ii < nz; ii++) source_zsum[i][ii] = f.source_zsum[i][ii];
      }
    }

  fieldObstArray = 0;

  if (fieldObstNum) {

    fieldObstArray = new VolumeObjClass[fieldObstNum]; 
    for (int i=0; i < fieldObstNum; i++)
      fieldObstArray[i].set( f.fieldObstArray[i] );
  }


  for(long i=0; i< size; i++) ptype[i] = f.ptype[i];
  for(int n=0; n < Region->get_surface_num() + fieldObstNum * 2 * DIMENSIONALITY; n++)  bccond[n] = f.bccond[n];

  ICa  = f.ICa;	
  D    = f.D;
  Time = f.Time;
  bgr  = f.bgr;
  Diff    = f.Diff;
  kuptake = f.kuptake;
}

//***************************************************************************

void FieldObj::allocate(const char *id)
{
ID = StrCpy(id);

 if ( !equal(ID,"Ca") ) source_num = 0;

  bccond = new int[Region -> get_surface_num() + 2 * DIMENSIONALITY * fieldObstNum];
  ptype  = new long[size];

  Currents     = new struct TermStruct[source_num];
  source_x0    = new int[source_num];
  source_nx    = new int[source_num];
  source_y0    = new int[source_num];
  source_ny    = new int[source_num];
  source_z0    = new int[source_num];
  source_nz    = new int[source_num];
  source_xsum  = new double *[source_num];
  source_ysum  = new double *[source_num];
  source_zsum  = new double *[source_num];
}

//****************************************************************************

FieldObj &FieldObj::operator=(const FieldObj &f)
{
 *((VectorObj *)this) = f;
 Time = f.Time;
 return *this;
}

//***************************************************************************

FieldObj &FieldObj::operator=(const VectorObj &v)
{
  *((VectorObj *)this) = v;

  //for (long i = 0; i < size; ++i) elem[i] = v[i]; 

return *this;
}

//***************************************************************************

FieldObj &FieldObj::operator=(const double val)
{
*((VectorObj *)this) = val;

//for (long i = 0; i < size; ++i) elem[i] = val; 
return *this;
}

//************************************************************************

FieldObj &FieldObj::operator=(double (*f0)(double, double, double) )
{
long   i = -1;

for (int iz = 0; iz < zsize; iz++)
  for (int iy = 0; iy < ysize; iy++)
    for (int ix = 0; ix < xsize; ix++)
      if ( ptype[++i] & VARY_MASK )
         elem[i] = (*f0)(xcoord[ix], ycoord[iy], zcoord[iz]);
 
return *this;
}


//************************************************************************************
//                   F I E L D   D E S T R U C T O R
//************************************************************************************

  FieldObj::~FieldObj()
    {
    //fprintf(stderr, "destroying %s\n",ID); fflush(stderr); 

    delete [] bccond; 
    delete [] ptype;

    for (int i = 0; i < source_num; i++)   {
      delete [] source_xsum[i]; delete [] source_ysum[i]; delete [] source_zsum[i];
      }
    
    if (fieldObstNum) delete [] fieldObstArray;
    
    delete [] source_xsum; delete [] source_ysum; delete [] source_zsum;

    delete [] source_x0;   delete [] source_nx;
    delete [] source_y0;   delete [] source_ny;
    delete [] source_z0;   delete [] source_nz; 

	delete [] ID; 

  }


//************************************************************************************
//                   F I E L D   I N I T I A L I Z A T I O N
//************************************************************************************


void  FieldObj::set_bc(int boxid, int bcx0, int bcx1, int bcy0, int bcy1, int bcz0, int bcz1)
{
  if (boxid >= Region->get_box_num() + fieldObstNum)  
    throw makeMessage("Cannot set boundary condition: box ID=%d > total number of volumes = %d", 
                      boxid, Region->get_box_num());

  if ( (bcx0 > bc_type_num) || (bcx1 > bc_type_num) || (bcy0 > bc_type_num) ||
       (bcy1 > bc_type_num) || (bcz0 > bc_type_num) || (bcz1 > bc_type_num) )  
    throw makeMessage("Bad boundary condition ID for geometry element #%d", boxid);
  
  bccond[boxid*2*DIMENSIONALITY    ] = bcx0;
  bccond[boxid*2*DIMENSIONALITY + 1] = bcx1; if (DIMENSIONALITY == 1) return;
  bccond[boxid*2*DIMENSIONALITY + 2] = bcy0;
  bccond[boxid*2*DIMENSIONALITY + 3] = bcy1; if (DIMENSIONALITY == 2) return;
  bccond[boxid*2*DIMENSIONALITY + 4] = bcz0;
  bccond[boxid*2*DIMENSIONALITY + 5] = bcz1;
}

//******************************************

void  FieldObj::set_bc(int boxid, const char *bcx0, const char *bcx1, const char *bcy0,
                                  const char *bcy1, const char *bcz0, const char *bcz1)
{
  set_bc( boxid, BCarray->bcidtoint(bcx0), BCarray->bcidtoint(bcx1), 
                 BCarray->bcidtoint(bcy0), BCarray->bcidtoint(bcy1), 
	             BCarray->bcidtoint(bcz0), BCarray->bcidtoint(bcz1) );
}

//******************************************

void  FieldObj::set_bc(int boxid, int bc) {
  set_bc( boxid, bc, bc, bc, bc, bc, bc);
}

//******************************************

void  FieldObj::set_bc(int boxid, const char *bcid) {
  set_bc( boxid, BCarray->bcidtoint(bcid) ); 
}

//************************************************************************************

signed long FieldObj::location_to_index(double x, double y, double z, bool errorMsg)
{
	long ind;

	ind = Grid->location_to_index(x, y, z);
	if (ptype[ind] & VARY_MASK)  return ind;

	int  ix, iy, iz;
	Grid->split(ind, ix, iy, iz);

	if (ix > 0) 
		if (ptype[ind-1] & VARY_MASK) return ind-1;
	if (ix < xsize-1) 
		if (ptype[ind+1] & VARY_MASK) return ind+1;
	if (DIMENSIONALITY == 1) { 
		if (errorMsg) throw makeMessage("location (%g) out of bounds", x);
		return -1;
	}

	if (iy > 0) 
		if (ptype[ind-xsize] & VARY_MASK) return ind-xsize;
	if (iy < ysize-1) 
		if (ptype[ind+xsize] & VARY_MASK) return ind+xsize;
	if (DIMENSIONALITY == 2) {
		if (errorMsg) throw makeMessage("location (%g, %g) out of bounds", x, y);
		return -1;
	}

	if (iz > 0) 
		if (ptype[ind-xysize] & VARY_MASK) return ind-xysize;
	if (iz < zsize-1) 
		if (ptype[ind+xysize] & VARY_MASK) return ind+xysize;
	if (DIMENSIONALITY == 3) {
		if (errorMsg) throw makeMessage("location (%g, %g, %g) out of bounds", x, y, z);
		return -1;
	}
	
	if (errorMsg) throw makeMessage("Can't parse location (%g, %g, %g)", x, y, z);
	return -1;
}

//*****************************************************************************
//                              S E T   S O U R C E
//************************************************************************************

void  FieldObj::set_source(int id, double x, double y, double z, double sigmax, double sigmay, double sigmaz)
{

  if (VERBOSE) {
                             fprintf(stderr, "    Setting source #%d at (%g",id+1, x);
     if (DIMENSIONALITY > 1) fprintf(stderr, ",%g",y);
     if (DIMENSIONALITY > 2) fprintf(stderr, ",%g",z);
                             fprintf(stderr, "), with half-width of (%g", sigmax);
     if (DIMENSIONALITY > 1) fprintf(stderr, ",%g",sigmay);
     if (DIMENSIONALITY > 2) fprintf(stderr, ",%g",sigmaz);
                             fprintf(stderr, ")");
                             fprintf(stderr, "\n");
  }

  if (id >= source_num)
    throw makeMessage("Bad source definition: source id=%d >= source_num=%d",id,source_num);

  long ind = long(location_to_index(x, y, z, 1));
  int ix, iy, iz;
  Grid->split(ind, ix, iy, iz);

  if (VERBOSE) {
                             fprintf(stderr, "    Closest node: (%d", ix);
     if (DIMENSIONALITY > 1) fprintf(stderr, ",%d",iy);
     if (DIMENSIONALITY > 2) fprintf(stderr, ",%d",iz);
                             fprintf(stderr, ") <--> (%.4g", xcoord[ix]);
     if (DIMENSIONALITY > 1) fprintf(stderr, ",%.4g",ycoord[iy]);
     if (DIMENSIONALITY > 2) fprintf(stderr, ",%.4g",zcoord[iz]);
                             fprintf(stderr, "), point_type=");  
							 print_type(ptype[ind]);
							 fprintf(stderr, "\n");
  }


  //*************************

  if (VERBOSE > 2)    { fprintf(stderr,"      Non-zero current on %s-nodes ",LABEL_DIM1); fflush(stderr); }
    set_source_xyz(id, x, ix, sigmax, xcoord, xgrid, dvx, xsize, source_x0, source_nx, source_xsum);

  if ( DIMENSIONALITY > 1) {
    if (VERBOSE > 2)  { fprintf(stderr,"      Non-zero current on %s-nodes ",LABEL_DIM2); fflush(stderr); }
    set_source_xyz(id, y, iy, sigmay, ycoord, ygrid, dvy, ysize, source_y0, source_ny, source_ysum);
  } else { 
    source_y0[id] = 0; source_ny[id] = 1;
    source_ysum[id] = new double[1]; source_ysum[id][0] = 1.0;
    }
  if ( DIMENSIONALITY > 2) {
    if (VERBOSE > 2)  { fprintf(stderr,"      Non-zero current on %s-nodes ",LABEL_DIM3); fflush(stderr); }
    set_source_xyz(id, z, iz, sigmaz, zcoord, zgrid, dvz, zsize, source_z0, source_nz, source_zsum);
    }
  else { 
    source_z0[id] = 0; source_nz[id] = 1;
    source_zsum[id] = new double[1]; source_zsum[id][0] = 1.0;
    }
}


//************************************************************************************
//                            S O U R C E   S P R E A D
//************************************************************************************

void  FieldObj::set_source_xyz(int id, double c0, int ix, double sigma, 
                               double *coord, double *grid, double *dvArray, int N,
                               int *source_coord0, int *interval_num, double **source_weights )
{
	  double *temp_weights = new double[N];
	  double check_sum = 0.0;
 
	  if (sigma < 0.25 * grid[ix])  sigma = 0.25 * grid[ix];  // if source is point-like, set half-width at 
                                                             // 40% of grid spacing
	  bool flag = false;    // to flag the entry into the region of non-zero current

	  for (int ii = 0; ii < N; ii++) 
			{
			double dx0 = 0.5 * grid[ii];
			double dx1 = 0.5 * grid[ii+1];

			double x0 = coord[ii] - c0 - dx0;
			double x1 = coord[ii] - c0 + dx1;

			double dv = dvArray[ii];  // default volume element

			     if (coord == xcoord && c0 == 0.0 && GEOMETRY1 >  0)            temp_weights[ii] = integrate( x0, x1, sigma, GEOMETRY1 );  // geom1 == r or geom1 == rho
			//else if (coord == xcoord && c0 >  0.0 && GEOMETRY1 >  0 && ii < 2)  temp_weights[ii] = 0.0;
			else if (coord == ycoord && c0 == 0.0 && GEOMETRY2 == 3)            temp_weights[ii] = integrate( x0, x1, sigma, 1 );  // geom2 == theta
			//else if (coord == ycoord && c0 >  0.0 && GEOMETRY2 == 3 && ii < 2)  temp_weights[ii] = 0.0;                            // geom2 == theta  
			else                                                             	temp_weights[ii] = integrate( x0, x1, sigma, 0 ); 

			check_sum += temp_weights[ii];
			temp_weights[ii] /= dv;
   
			if ( temp_weights[ii] > SOURCE_EPS )   // average between ii-1/2 and ii+1/2 > eps
			  {
			   if ( !flag )  
				 {
					 flag = true;
					 source_coord0[id] = ii;
					 interval_num[id] = 1;
				 }
			   else interval_num[id]++;
			   if (VERBOSE > 2) fprintf(stderr, "%d(%.5g):%.4g ", ii, coord[ii], temp_weights[ii]*dv);
			  }
		} // end for (ii = 0; ii <= N; ii++) 

	  //*************************

	  if (interval_num[id] <= 0) throw makeMessage(" Current source %d has no overlap with the domain", id);

	  source_weights[id] = new double[interval_num[id]];

	  for (int i = 0; i < interval_num[id]; i++) 
		source_weights[id][i] = temp_weights[i + source_coord0[id]];

	  if (fabs(check_sum - 1) > 1.0e-4 && VERBOSE > 2) 
		  fprintf(stderr," => Total = %.2e ", check_sum );

	  if (VERBOSE > 2) fprintf(stderr, "\n");

	  delete [] temp_weights;
}

//************************************************************************************
//                      O L D    S O U R C E   S P R E A D
//************************************************************************************
/*  
 
  Evaluating the integral of exp(-x^2/sigma^2) x^dim dx
  Normalization factor is evaluated by setting to 1 the intergal from 0 to inf:

  Int( exp(-x^2/sigma^2) x^dim dx, 0, inf ) = sigma^(dim+1) Int( exp(-y*y) y^dim dy, 0, inf) =
  = sigma^(dim+1) / 2 * Int( exp(-y*y) y^(dim-1) d(y^2), 0 , inf) =
  = sigma^(dim+1) / 2 * Gamma( (dim+1) / 2)
  
  So the normalization factor is 2 / (sigma^(dim+1) Gamma( (dim+1) / 2) ).
  For dim=0, integration is from -inf to inf, not from 0 to inf, so the factor is half that

  dim = 0: 1/(sigma sqrt(pi))
  dim = 1: 2/(sigma^2);
  dim = 2: 4/(sigma^3 sqrt(pi)) 

************************************************************************************/

double gaussian(double x0, double x1, double sigma, int dim)
{
  double eps = 1.0e-10;

  const double a[3] = { 1.0 / ( sigma * sqrt(pi) ),
                        2.0 / ( sigma * sigma    ),
			            4.0 / ( sigma * sigma * sigma * sqrt(pi) ) };
  double b = -1.0 / (sigma * sigma);

  long   N = 2;
  double resN = 0, res2N = 0, res4N = 0;
  double x, dx, sum;

  do
    {
    resN = res2N;
    N *= 2;
    dx = (x1 - x0) / double(N);
    sum = 0.0;
    for (int n = 1; n < N; n++)
      {
      x = x0 + dx * double(n);
      sum += ((n & 1) + 1) * exp( b * x * x ) * pow(x, dim);
      }
    res2N = dx / 3.0 * a[dim] * ( 2 * sum + exp( b * x0 * x0 ) * pow(x0, dim) + exp( b * x1 * x1 ) * pow(x1, dim) );
    res4N = (16.0 * res2N - resN) / 15.0;
    }
  while ( fabs(res4N - res2N) > eps * fabs(res4N + res2N) && fabs(res4N + res2N) > 1e-12);

  return res4N;
}

//************************************************************************************/

double square(double x0, double x1, double sigma, int dim)
{
  if (x0 >= sigma || x1 <= -sigma) return 0;

  double norm = 1.0 / pow(sigma, double(dim + 1) );
  if (dim == 0) norm /= 2.0;

  if (x0 >= 0) {
     if (x1 < sigma) return  norm * (pow(x1,    dim + 1) - pow(x0, dim + 1) );
                else return  norm * (pow(sigma, dim + 1) - pow(x0, dim + 1) );
  } else {  // must be 1D (dim=0) case; already know that x1 > -sigma
     if (x0 < -sigma) if (x1 < sigma) return norm * (x1 + sigma);
                      else return 1.0;
                 else if (x1 < sigma) return norm * (x1 - x0);
                      else return norm * (sigma - x0);
  }
 
}


//***************************************************************

void  FieldObj::add_sources(VectorObj &f, double dt)
{
int    i, ix, iy, iz;
long   xind, yind, zind;
double addx, addy, charge, ica;

  ICa = 0.0; 

  for (i=0; i < source_num; i++) {

    ica = SAME_CURRENT ? Current->Evaluate() : Evaluate(Currents[i]);
    if (fabs(ica) < MIN_CURRENT) continue;
	charge = ica * dt;
    ICa += ica;
    int nx = source_nx[i];
    int ny = source_ny[i];
    int nz = source_nz[i];
    xind = source_x0[i] + xsize * source_y0[i] + xysize * source_z0[i];

    for (ix = 0; ix < nx; ix++, xind ++) {
      addx = charge * source_xsum[i][ix];
      for (iy = 0, yind = xind; iy < ny; iy++, yind += xsize) {
		addy = addx * source_ysum[i][iy];
	    for (iz = 0, zind = yind; iz < nz; iz++, zind += xysize) //{ 
			   f.elem[zind] += addy * source_zsum[i][iz];
      }
    }
  }

}

//***************************************************************


void  FieldObj::adjust_sources(FieldObj f)
{
int    i, ix, iy, iz;
long   xind, yind, zind;
double addx, addy;
double newAverage, volume, source_adjust;

  ICa = 0.0;

  for (i=0; i < source_num; i++) {

    int   nx = source_nx[i];
    int   ny = source_ny[i];
    int   nz = source_nz[i];
	f = 0.0;

    xind = source_x0[i] + xsize * source_y0[i] + xysize * source_z0[i];

    for (ix = 0; ix < nx; ix++, xind ++) {
      addx = source_xsum[i][ix];
      for (iy = 0, yind = xind; iy < ny; iy++, yind += xsize) {
		addy = addx * source_ysum[i][iy];
	    for (iz = 0, zind = yind; iz < nz; iz++, zind += xysize)  
			   f.elem[zind] += addy * source_zsum[i][iz];
		}
	}

	newAverage = f.average(&volume);
	source_adjust = fabs( 1 / (volume * newAverage) );

	for (ix = 0; ix < nx; ix++)
		source_xsum[i][ix] = source_xsum[i][ix] * source_adjust;

	if (VERBOSE) fprintf(stderr, "#### Source #%d boundary spill-over adjustment factor = %g\n", i+1, source_adjust);
	if (source_adjust > 50) 
		fprintf(stderr,"\n *** WARNING: Check source #%d width or resolution: spill-over adjustment multiplier very large (%g)\n", i+1, source_adjust);
  }
}


//************************************************************************
//            B O U N D A R Y   I N I T I A L I Z A T I O N
//************************************************************************
/*
   RegionObj::point_type(int,int,int) bit fields are:
   bit # is:      ...8       |   7    |   6  |  5   |  4   |  3   |  2   |  1  
   meaning:    box id number | inside | zmax | zmin | ymax | ymin | xmax | xmin

  ***************************************************************

   ptype array bit fields are:

   bit # is:  ...10    |    9    |    8      |   7    |   6  |  5   |  4   |  3   |  2   |  1  
   meaning:  bc_type # | outside | dirichlet | inside | zmax | zmin | ymax | ymin | xmax | xmin
  
   bit # is: |   30...24    |   23...17    |   16...10    |
   meaning:  | z bc_type id | y bc_type id | x bc_type id |

  ***************************************************************/

void FieldObj::init_boundary()    {

	long      i;
	int       p, boxid;
	int       ix, iy, iz;

	for (i = 0; i < size; i++) ptype[i] =  _INSIDE_;   // default: point inside

	for (i = 0; i < size; i++)
	{
		Grid->split(i, ix, iy, iz);
		p = point_type(ix, iy, iz); 
		boxid = p >> 7;
		p = p & INSIDE_MASK;              // remove box id (keep bits 1 through 7)

		if (p == 0)                       // point is outside
		{
			ptype[i] = _OUTSIDE_; 
			elem[i] = 0.0;
		}
		else if (p & SURF_MASK)           // point on the surface
		{
			if (p & SURF_XMIN) ptype[i] |= (get_bc(boxid,0) << XSHIFT) + SURF_XMIN;
			if (p & SURF_XMAX) ptype[i] |= (get_bc(boxid,1) << XSHIFT) + SURF_XMAX;
			if (p & SURF_YMIN) ptype[i] |= (get_bc(boxid,2) << YSHIFT) + SURF_YMIN;
			if (p & SURF_YMAX) ptype[i] |= (get_bc(boxid,3) << YSHIFT) + SURF_YMAX;
			if (p & SURF_ZMIN) ptype[i] |= (get_bc(boxid,4) << ZSHIFT) + SURF_ZMIN;
			if (p & SURF_ZMAX) ptype[i] |= (get_bc(boxid,5) << ZSHIFT) + SURF_ZMAX;
		} // end if (p & SURF_MASK)
	} // end for loop


	bool flag = true;
	while(flag) {
		while (flag) {
			flag = false;

			if (VERBOSE > 3) fprintf(stderr, " ## Building geometry: second pass (detecting min-max boundaries)\n");
			for (i = 0; i < size; i++) {  // SECOND PASS: move min-max double boundary by one point inward, in all directions
				p = ptype[i];
				if (((p & SURF_XMIN) && (p & SURF_XMAX)) ||
					((p & SURF_YMIN) && (p & SURF_YMAX)) ||
					((p & SURF_ZMIN) && (p & SURF_ZMAX)) )   {

						Grid->split(i, ix, iy, iz);

						if ( (p & SURF_XMIN) && !(p & SURF_XMAX) && (ix+1 < Grid->xsize) ) if (ptype[i+1]      & _INSIDE_)
						{ ptype[i+1] |= ( SURF_XMIN | (((p >> XSHIFT) & 127 ) << XSHIFT) );  ptype[i] = _OUTSIDE_; flag = true;  }  

						if ( (p & SURF_XMAX) && !(p & SURF_XMIN) && (ix > 0) )             if (ptype[i-1]      & _INSIDE_)
						{ ptype[i-1] |= ( SURF_XMAX | (((p >> XSHIFT) & 127 ) << XSHIFT) );  ptype[i] = _OUTSIDE_; flag = true; }  

						if ( (p & SURF_YMIN) && !(p & SURF_YMAX) && (iy+1 < Grid->ysize) ) if (ptype[i+ysize]  & _INSIDE_)
						{ ptype[i+ysize] |= ( SURF_YMIN | (((p >> YSHIFT) & 127 ) << YSHIFT) );  ptype[i] = _OUTSIDE_; flag = true;  } 

						if ( (p & SURF_YMAX) && !(p & SURF_YMIN) && (iy > 0) )             if (ptype[i-ysize]  & _INSIDE_)
						{ ptype[i-ysize] |= ( SURF_YMAX | (((p >> YSHIFT) & 127 ) << YSHIFT) );  ptype[i] = _OUTSIDE_; flag = true; }  

						if ( (p & SURF_ZMIN) && !(p & SURF_ZMAX) && (iz+1 < Grid->zsize) ) if (ptype[i+xysize] & _INSIDE_)
						{ ptype[i+xysize] |= ( SURF_ZMIN | (((p >> ZSHIFT) & 127 ) << ZSHIFT) );  ptype[i] = _OUTSIDE_; flag = true;  }  

						if ( (p & SURF_ZMAX) && !(p & SURF_ZMIN) && (iz > 0) )             if (ptype[i-xysize] & _INSIDE_)
						{ ptype[i-xysize] |= ( SURF_ZMAX | (((p >> ZSHIFT) & 127 ) << ZSHIFT) );  ptype[i] = _OUTSIDE_; flag = true;  }

						ptype[i] = _OUTSIDE_;
				}
			} // end for loop
		}

		if (VERBOSE > 3) fprintf(stderr, " ## Building geometry: third and fourth pass (eliminating slit gaps)\n");
		for (i = 0; i < size; i++) {    // THIRD PASS: remove slit-gaps between boundaries and obstacles, marked as _OUTSIDE_ at third pass
			p = ptype[i];               // Note: bc type for filled-in gap boundaries = 0 (Noflux) or inherited from other boundaries
			Grid->split(i, ix, iy, iz);
			if (p & _OUTSIDE_) 
			{
				if (ix > 0)       if ( (ptype[i-1]      & VARY_MASK) && !(ptype[i-1]      & SURF_XMAX) ) { ptype[i-1]      |= SURF_XMAX; flag = true; }
				if (ix < xsize-1) if ( (ptype[i+1]      & VARY_MASK) && !(ptype[i+1]      & SURF_XMIN) ) { ptype[i+1]      |= SURF_XMIN; flag = true; }
				if (iy > 0)       if ( (ptype[i-xsize]  & VARY_MASK) && !(ptype[i-xsize]  & SURF_YMAX) ) { ptype[i-xsize]  |= SURF_YMAX; flag = true; }
				if (iy < ysize-1) if ( (ptype[i+xsize]  & VARY_MASK) && !(ptype[i+xsize]  & SURF_YMIN) ) { ptype[i+xsize]  |= SURF_YMIN; flag = true; }
				if (iz > 0)       if ( (ptype[i-xysize] & VARY_MASK) && !(ptype[i-xysize] & SURF_ZMAX) ) { ptype[i-xysize] |= SURF_ZMAX; flag = true; }
				if (iz < zsize-1) if ( (ptype[i+xysize] & VARY_MASK) && !(ptype[i+xysize] & SURF_ZMIN) ) { ptype[i+xysize] |= SURF_ZMIN; flag = true; }
			}
			                            // FOURTH PASS: remove false boundaries (max boundary immediately followed by min boundary)
			if (ix > 0 && ix < xsize-1)   if ( (ptype[i] & SURF_XMAX) && (ptype[i+1]      & SURF_XMIN) ) { ptype[i] ^= SURF_XMAX;  ptype[i+1]      ^= SURF_XMIN; flag = true; }
			if (iy > 0 && iy < ysize-1)   if ( (ptype[i] & SURF_YMAX) && (ptype[i+xsize]  & SURF_YMIN) ) { ptype[i] ^= SURF_YMAX;  ptype[i+xsize]  ^= SURF_YMIN; flag = true; }
			if (iz > 0 && iz < zsize-1)   if ( (ptype[i] & SURF_ZMAX) && (ptype[i+xysize] & SURF_ZMIN) ) { ptype[i] ^= SURF_ZMAX;  ptype[i+xysize] ^= SURF_ZMIN; flag = true; }

		} // end for loop
	}

	if (VERBOSE > 5 && xsize < 40) {    // Pretty-print the geometry using ASCII characters
		for (i = 0; i < size; i++)  {
			Grid->split(i, ix, iy, iz);
			p = ptype[i];
			if (ix == 0) { fprintf(stderr, "\n");
			if (iy == 0)   fprintf(stderr, "\n *** Geometry on surface z=%g:\n\n", zcoord[iz]); }
			if (p & 256)   fprintf(stderr, "    ");
			else if (p & SURF_MASK) fprintf(stderr, "%d%d%d ", p&3, (p&12) >> 2, (p&48) >> 4);
			else fprintf(stderr, " *  ");
		} 

		fprintf(stderr, "\n");
	}
}



//***********************************************************************************
//*******************************  3D AVERAGES  ************************************* 
//***********************************************************************************

double FieldObj::average3D(double *V)
{
double s = 0.0, v = 0.0, dv, vx, vxy;
int    ix, iy, iz;
long   i;

  for(ix = 0; ix < xsize; ix++)
  {
	  vx = dvx[ix];
	  for(iy = 0; iy < ysize; iy++) 
	  {
		vxy = vx * dvy[iy];
		for(iz = 0, i = ix + iy * xsize; iz < zsize; iz++, i += xysize) 
		  if ( ! (ptype[i] & _OUTSIDE_) ) {
			  v += ( dv = vxy * dvz[iz] );
			  s += ( dv * elem[i]       );  
		  }
	  }
  }
  if (V) *V = v;
  return s / v;
}

//***********************************************************************************
//*******************************  2D AVERAGES  ************************************* 
//***********************************************************************************

double FieldObj::average2D(double *V)
{
double s = 0.0, v = 0.0, dv, vx;
int    ix, iy;
long   i;

  for(ix = 0; ix < xsize; ix++) 
  {
    vx = dvx[ix];
    for(iy = 0, i = ix; iy < ysize; iy++, i += xsize) 
      if ( ! (ptype[i] & _OUTSIDE_) )
		{
		  v += ( dv = vx * dvy[iy] );
	      s += ( dv * elem[i]      );  
		}
  }
  
 if (V) *V = v;
 return s / v;
}

//***********************************************************************************
//*******************************  1D AVERAGES  ************************************* 
//***********************************************************************************

double FieldObj::average1D(double *V)
{
double s = 0.0, v = 0.0, dv;

for(int i = 0; i < xsize; i++)
  if ( ! (ptype[i] & _OUTSIDE_) ) {
    v += ( dv = dvx[i]  );
    s += ( dv * elem[i] );
  }

 if (V) *V = v;
 return s / v;
}


//***********************************************************************************
/*
double FieldObj::sum()
{
double s = 0.0;
for(long i=0; i<size; i++) if (ptype[i] & REGION_MASK) s += elem[i];
return s;
}
*/

//***********************************************************************************

/*
double FieldObj::L2norm()
{
double s = 0.0;
for(long i=0; i<size; i++) if (ptype[i] & REGION_MASK) s += elem[i] * elem[i];
return s / double(size);
}
*/

//***********************************************************************************
/*
void print3d(long *var)
{
for (int iz=0; iz < zsize; iz++, fprintf(stderr,("\n"))
 for (int iy=0; iy < ysize; iy++, fprintf(stderr,("\n"))
  for (int ix=0; ix < xsize; ix++) 
      fprintf(stderr,("%ld ",var[ix + iy * xsize + iz * xysize]);
}
*/  
//***********************************************************************************
/*
void FieldObj::print3d(int n)
{
        char s1[20], s2[20];
        strcpy(s1,"%.");
        sprintf(s2,"%d",n);
        strcat(s1,s2);
        strcat(s1,"f ");

for (int iz=0; iz<zsize; iz++, fprintf(stderr,"\n"))
 for (int iy=0; iy<ysize; iy++, fprintf(stderr,"\n"))
  for (int ix=0; ix<xsize; ix++) fprintf(stderr,s1,elem[ix+iy*xsize+iz*xysize]);

}
*/
//***********************************************************************************

void FieldObj::print_ptype()
{
for (int iz=0; iz < zsize; iz++, fprintf(stderr,"\n"))
 for (int iy=0; iy < ysize; iy++, fprintf(stderr,"\n"))
  for (int ix=0; ix < xsize; ix++) 
      print_ptype(ix + iy * xsize + iz * xysize);
}

//***************************************************************************************
//                              X   D I R E C T I O N
//***************************************************************************************
/*
 Find vector Result_new (to be stored in array Result) satisfying

 (1 + c * Dx^2 - Lin) Result_new = ( cx * Dx^2 + cy * Dy^2 + cz * Dz^2 ) This + Result 

Spherical:
                                                        1                         1
 (1 + c * Dr^2 - Lin) Result_new = ( cx * Dr^2 + cth * --- * Dth^2 + cphi * ------------- * Dphi^2 ) This + Result 
                                                       r^2                  sin(th)^2 r^2
Cylindrical:
                                                         1                  
 (1 + c * Dr^2 - Lin) Result_new = ( cx * Dr^2 + cphi * --- * Dphi^2 + cz * Dz^2 ) This + Result 
                                                        r^2                 

 | diag[1] sup[1]   0      0 :                    |                | right[1] |
 | sub[2]  diag[2] sup[2]  0 :            0       |                |          |
 |...........................:....................| x Result_new = |          |
 |                           : diag[n-1] sup[n-1] |                |          |
 |            0              : sub[n]    diag[n]  |                | right[n] |

 Operator Dx^2:   Dx^2 f[i] == xminus[i] f[ix-1] + xcenter[i] f[i] + xplus[i] f[ix+1]
                          i == (ix,iy,iz)

***************************************************************************************/

void FieldObj::Run3Dx(double c, double cx, double cy, double cz, double *Lin, double *Result)
{
int       ix = 0, iy = 0, iz = 0;
long      i = 0, istart = 0, bcvar;
int       n = 0;
double    xminus, xplus, xcenter, xcoef;
double    yminus, yplus, ycenter, ycoef;
double    zminus, zplus, zcenter, zcoef;
double    ri, invSine = 1.0, invR = 1.0;
const int angleY = GEOMETRY2 >> 1;  // == 1 if ycoord == phi or theta

if ( D <= 0.0 )  {
  for (i = 0; i < size; i ++)  if (ptype[i] & VARY_MASK) 
      Result[i] /= (1.0 -  (Lin ? Lin[i] : 0.0) );
  return;
}

if (GEOMETRY == SPHERICAL3D) 
	invSine = 1.0 / sin(ycoord[0]);    // "Phi" gradient : d/(r sin(theta) dphi)

//****************************************************************************************

do  { // loop over iy, iz:  while (iz < zsize)  NOTE: invSine computed at bottom of loop

  while ( (bcvar = ptype[i]) & VARY_MASK )  {  // loop over ix
		n++;
		if ( bcvar & SURF_XMIN ) istart = i; 
		
		if ( angleY )  invR = 1 / xcoord[ix];

		nablaX(i, ix, bcvar, xminus, xcenter, xplus, xcoef);
		nablaY(i, iy, bcvar, yminus, ycenter, yplus, ycoef, invR);
		nablaZ(i, iz, bcvar, zminus, zcenter, zplus, zcoef, invR*invSine);

		ri = Result[i] + (cx - c) * xcoef + cy * ycoef + cz * zcoef + (cx * xcenter + cy * ycenter + cz * zcenter) * elem[i];
		if ( ! (bcvar & SURF_XMIN) )  ri += cx * xminus * elem[i-1]; 
		if ( ! (bcvar & SURF_XMAX) )  ri += cx * xplus  * elem[i+1];     
		if ( ! (bcvar & SURF_YMIN) )  ri += cy * yminus * elem[i-xsize]; 
		if ( ! (bcvar & SURF_YMAX) )  ri += cy * yplus  * elem[i+xsize];    
		if ( ! (bcvar & SURF_ZMIN) )  ri += cz * zminus * elem[i-xysize]; 
		if ( ! (bcvar & SURF_ZMAX) )  ri += cz * zplus  * elem[i+xysize];    
 
		diag[n] = 1 - (Lin ? Lin[i] : 0.0) + c * xcenter;
		sup[n] = c * xplus;  sub[n] = c * xminus;
		right[n] = ri; 

		i ++;
		if (++ix == xsize) break;
  } // while bcvar > 0 and ix < xsize

  if (n)  { 
    tridiag(n);
    for (int ii = 1, ind = istart; ii<=n; ii++, ind ++) Result[ind] = right[ii];
    n = 0;
  }

  if (ix < xsize) 
    while ( !(ptype[i] & VARY_MASK) ) {
      i ++; 
      if ( (++ix) == xsize ) break; 
    }

  if (ix == xsize) { 
     ix = 0; i += (xysize - xsize);
     if (++iz == zsize) { iz = 0; iy++; i += (xsize - size); 
						  if (GEOMETRY3 && iy < ysize) invSine = 1.0 / sin(ycoord[iy]); } 
  }
  
} while (iy < ysize);   

}


//***************************************************************************************
//                               Y   D I R E C T I O N
//
//  (1 + c * Dy^2 - Lin) Result_new = ( cx * Dx^2 + cy * Dy^2 + cz * Dz^2 ) This + Result
/*
Spherical:
           1                                                 1                         1
 (1 + c * --- Dth^2 - Lin) Result_new = ( cx * Dr^2 + cth * --- * Dth^2 + cphi * ------------- * Dphi^2 ) This + Result 
          r^2                                               r^2                  sin(th)^2 r^2
Cylindrical:
           1                                                   1                  
 (1 + c * --- Dphi^2 - Lin) Result_new = ( cx * Dr^2 + cphi * --- * Dphi^2 + cz * Dz^2 ) This + Result 
          r^2                                                 r^2                 
*/
//***************************************************************************************

void FieldObj::Run3Dy(double c, double cx, double cy, double cz, double *Lin, double *Result)
{
int       ix = 0, iz = 0, iy = 0;
long      i = 0, istart = 0, bcvar;
int       n = 0;
double    xminus, xplus, xcenter, xcoef;
double    yminus, yplus, ycenter, ycoef;
double    zminus, zplus, zcenter, zcoef;
double    ri, invR=1.0, invSine = 1.0;
const int angleY = GEOMETRY2 >> 1;  // == 1 if ycoord == phi or theta

if ( D <= 0.0 )  {
  for (i = 0; i < size; i ++) if (ptype[i] & VARY_MASK) 
     Result[i] /= (1.0 -  (Lin ? Lin[i] : 0.0) );
  return;
}

if (angleY)  invR = 1.0 /  xcoord[0];

do  {  // loop over iz, ix: while (ix < xsize)  ** NOTE: invR updated at bottom of loop
	  
	  while ( (bcvar = ptype[i]) & VARY_MASK ) { // loop over iy
			n++;
			if ( bcvar & SURF_YMIN ) istart = i; 

			if (GEOMETRY == SPHERICAL3D) invSine = 1.0 / sin(ycoord[iy]); // "Phi" gradient : d/(r sin(theta) dphi)

			nablaX(i, ix, bcvar, xminus, xcenter, xplus, xcoef);
			nablaY(i, iy, bcvar, yminus, ycenter, yplus, ycoef, invR);
			nablaZ(i, iz, bcvar, zminus, zcenter, zplus, zcoef, invR*invSine);
    
			ri = Result[i] + (cy - c) * ycoef + cx * xcoef + cz * zcoef + (cx * xcenter + cy * ycenter + cz * zcenter) * elem[i];
			if ( ! (bcvar & SURF_XMIN) )  ri += cx * xminus * elem[i-1]; 
			if ( ! (bcvar & SURF_XMAX) )  ri += cx * xplus  * elem[i+1];     
			if ( ! (bcvar & SURF_YMIN) )  ri += cy * yminus * elem[i-xsize]; 
			if ( ! (bcvar & SURF_YMAX) )  ri += cy * yplus  * elem[i+xsize];    
			if ( ! (bcvar & SURF_ZMIN) )  ri += cz * zminus * elem[i-xysize]; 
			if ( ! (bcvar & SURF_ZMAX) )  ri += cz * zplus  * elem[i+xysize];    

			diag[n] = 1 - (Lin ? Lin[i] : 0.0) + c * ycenter;
			sup[n] = c * yplus;  sub[n] = c * yminus;
			right[n] = ri; 

			i += xsize;  
			if (++iy == ysize) break;
	  } // while bcvar > 0 and iy < ysize

	  if (n > 0)  {
		  tridiag(n);
		  for (int ii = 1, ind = istart; ii<=n; ii++, ind += xsize) Result[ind] = right[ii]; 
		  n = 0;
	  }

	  if (iy < ysize)
		  while ( !(ptype[i] & VARY_MASK) ) {  
			 i += xsize; if ( (++iy) == ysize ) break; 
	  }

	  if (iy == ysize)    {
		  iy = 0 ;
		  if (++iz == zsize) { i += (1 - size); iz = 0; ix++; 
							   if (angleY && ix < xsize) invR = 1.0 /  xcoord[ix]; }
	  }

} while (ix < xsize); 

}


//***************************************************************************************
//                              Z   D I R E C T I O N
//***************************************************************************************
//
// (1 + c * Dz^2 - Lin) Result_new = ( cx * Dx^2 + cy * Dy^2 + cz * Dz^2 ) This + Result 
/*
Spherical:
                1                                                     1                         1
 (1 + c * ------------ Dth^2 - Lin) Result_new = ( cx * Dr^2 + cth * --- * Dth^2 + cphi * ------------- * Dphi^2 ) This + Result 
          sin(th)^2 r^2                                              r^2                  sin(th)^2 r^2
Cylindrical:
                                                         1                  
 (1 + c * Dz^2 - Lin) Result_new = ( cx * Dr^2 + cphi * --- * Dphi^2 + cz * Dz^2 ) This + Result 
                                                        r^2                 
*/
//***************************************************************************************

void FieldObj::Run3Dz(double c, double cx, double cy, double cz, double *Lin, double *Result)

{
int       iz = 0, ix = 0, iy = 0;
long      i = 0, istart = 0, bcvar;
int       n = 0;
double    xminus, xplus, xcenter, xcoef;
double    yminus, yplus, ycenter, ycoef;
double    zminus, zplus, zcenter, zcoef;
double    ri, invR = 1.0, invSine = 1.0;
const int angleY = GEOMETRY2 >> 1;  // == 1 if ycoord == phi or theta

if ( D <= 0.0 )  {
  for (i = 0; i < size; i ++) if (ptype[i] & VARY_MASK) 
    Result[i] /= (1.0 - (Lin ? Lin[i] : 0.0) );
  return;
}

  if (GEOMETRY3) invSine = 1.0 / sin(ycoord[0]); // "Phi" gradient : d/(r sin(theta) dphi)

  do  {  // loop over ix, iy:  while (iy < ysize): ** NOTE: invSine updated at bottom of loop
	  if (angleY) invR = 1.0 /  xcoord[ix];

	  while ( (bcvar = ptype[i]) & VARY_MASK ) {  // loop over iz
			n++;

			if ( bcvar & SURF_ZMIN ) istart = i; 
			nablaX(i, ix, bcvar, xminus, xcenter, xplus, xcoef);
			nablaY(i, iy, bcvar, yminus, ycenter, yplus, ycoef, invR);
			nablaZ(i, iz, bcvar, zminus, zcenter, zplus, zcoef, invR*invSine);
    
			ri = Result[i] + (cz - c) * zcoef + cx * xcoef + cy * ycoef + (cx * xcenter + cy * ycenter + cz * zcenter) * elem[i];
			if ( ! (bcvar & SURF_XMIN) )  ri += cx * xminus * elem[i-1]; 
			if ( ! (bcvar & SURF_XMAX) )  ri += cx * xplus  * elem[i+1];     
			if ( ! (bcvar & SURF_YMIN) )  ri += cy * yminus * elem[i-xsize]; 
			if ( ! (bcvar & SURF_YMAX) )  ri += cy * yplus  * elem[i+xsize];    
			if ( ! (bcvar & SURF_ZMIN) )  ri += cz * zminus * elem[i-xysize]; 
			if ( ! (bcvar & SURF_ZMAX) )  ri += cz * zplus  * elem[i+xysize];    

			diag[n] = 1 - (Lin ? Lin[i] : 0.0) + c * zcenter;
			sup[n] = c * zplus;  sub[n] = c * zminus;
			right[n] = ri; 

			i += xysize; ++iz;
			if (iz == zsize) break;

	  } // while ptype > 0 and iz < zsize

	  if (n > 0)    {  
		tridiag(n);
		for (int ii = 1, ind = istart; ii<=n; ii++, ind += xysize) Result[ind] = right[ii]; 
		n = 0;
	  } 
 
	  if (iz < zsize)
		while ( !(ptype[i] & VARY_MASK) ) {
		   i += xysize; if ( (++iz) == zsize) break; 
		}

	  if (iz == zsize)  {
		  iz = 0; i += (1 - size);
		  if (++ix == xsize) { ix = 0; iy ++; 
							   if (GEOMETRY3 && iy < ysize) invSine = 1.0 / sin(ycoord[iy]); }
	  }

  } while (iy < ysize);

}


     
/***************************************************************************************
 
                     CYLINDRICAL GEOMETRY RUN, "r" DIRECTION
                                                            
 Solving  (1 + c Dr^2 - Lin) Result_new = ( cr Dr^2 + cz Dz^2 ) This + Result
                                                          
***************************************************************************************

                     CONICAL GEOMETRY RUN, "r" DIRECTION:
                                                           1 
 Solving  (1 + c Dr^2 - Lin) Result_new = ( cr Dr^2 + cth --- Dth^2 ) This + Result
                                                          r^2
***************************************************************************************

                     POLAR GEOMETRY RUN, "r" DIRECTION:
                                                            1 
 Solving  (1 + c Dr^2 - Lin) Result_new = ( cr Dr^2 + cphi --- Dphi^2 ) This + Result
                                                           r^2
***************************************************************************************/

void FieldObj::Run2Dx( double c, double cr, double cz, double *Lin, double *Result)
{
long      i = 0, istart = 0, bcvar;
int       ix = 0, iy = 0, n = 0;
double    rplus, rminus, rcenter, rcoef;
double    zplus, zminus, zcenter, zcoef;
const int angleY = GEOMETRY2 >> 1;  // == 1 for conical or polar geometry
double    ri, invR = 1.0;

if ( D <= 0.0 )  {
  for (i = 0; i < size; i ++)
    if (ptype[i] & VARY_MASK) Result[i] /= (1.0 - (Lin ? Lin[i] : 0.0) );
  return;
}

do  { // while iy < ysize
	  while ( (bcvar = ptype[i]) & VARY_MASK )   { // increment ix
			n++;
			if ( bcvar & SURF_XMIN )   istart = i; 
			if (angleY) invR = 1.0 / xcoord[ix];  // Scale factor for angular y-direction: d/(r dy)
			nablaX(i, ix, bcvar, rminus, rcenter, rplus, rcoef);     // calculating D_r^2
			nablaY(i, iy, bcvar, zminus, zcenter, zplus, zcoef, invR);

			ri = Result[i] + (cr - c) * rcoef + cz * zcoef + (cr * rcenter + cz * zcenter) * elem[i];
			if ( ! (bcvar & SURF_XMIN) )  ri += cr * rminus * elem[i-1]; 
			if ( ! (bcvar & SURF_XMAX) )  ri += cr * rplus  * elem[i+1];     
			if ( ! (bcvar & SURF_YMIN) )  ri += cz * zminus * elem[i-xsize]; 
			if ( ! (bcvar & SURF_YMAX) )  ri += cz * zplus  * elem[i+xsize];    
 
			diag[n] = 1 - (Lin ? Lin[i] : 0.0) + c * rcenter;
			sup[n] = c * rplus;  sub[n] = c * rminus;
			right[n] = ri; 
			++i;
			if (++ix == xsize) break;

	  } // while bcvar > 0 and i < xsize
 
	  if (n)  {
		tridiag(n);
		for (int ii = 1, ind = istart; ii<=n; ii++, ind ++) Result[ind] = right[ii];
		n = 0;
	  }

	  if (ix < xsize)
		while ( !(ptype[i] & VARY_MASK) )   {
		   ++i;
		   if ( ++ix == xsize) break; 
		}

	  if (ix == xsize) {
		ix = 0;
		iy ++;
	  }
} while (iy < ysize);

}

/***************************************************************************************
 
                     CYLINDRICAL GEOMETRY RUN, "z" DIRECTION
                                                            
 Solving  (1 + c Dz^2 - Lin) Result_new = ( cr Dr^2 + cz Dz^2 ) This + Result
                                                          
***************************************************************************************

                     CONICAL GEOMETRY RUN, "theta" DIRECTION:
                  1                                             1 
 Solving  (1 + c --- Dth^2 - Lin) Result_new = ( cr Dr^2 + cth --- Dth^2 ) This + Result
                 r^2                                           r^2

***************************************************************************************

                     POLAR GEOMETRY RUN, "phi" DIRECTION:
                  1                                               1 
 Solving  (1 + c --- Dphi^2 - Lin) Result_new = ( cr Dr^2 + cphi --- Dphi^2 ) This + Result
                 r^2                                             r^2
***************************************************************************************/

void FieldObj::Run2Dy( double c, double cr, double cz, double *Lin, double *Result)
{
long      i = 0, istart = 0, bcvar;
int       ix = 0, iy = 0, n = 0;
double    ri, invR = 1.0;
double    rplus, rminus, rcenter, rcoef;
double    zplus, zminus, zcenter, zcoef;
const int angleY = GEOMETRY2 >> 1;

if ( D <= 0.0 )  {
  for (i = 0; i < size; i ++)
    if (ptype[i] & VARY_MASK) Result[i] /= (1.0 - (Lin ? Lin[i] : 0.0) );
  return;
}

do  { // while ix < xsize
	  if (angleY) invR = 1.0 / xcoord[ix];  // Scale factor for angular y-direction: d/(r dy)

	  while ( (bcvar = ptype[i]) & VARY_MASK ) { // increment iy

			n++;
			if ( bcvar & SURF_YMIN )  istart = i; 

			nablaY(i, iy, bcvar, zminus, zcenter, zplus, zcoef, invR);
			nablaX(i, ix, bcvar, rminus, rcenter, rplus, rcoef);     // calculating D_r^2  
   
			ri = Result[i] + cr * rcoef + (cz - c) * zcoef + (cr * rcenter + cz * zcenter) * elem[i];
			if ( ! (bcvar & SURF_XMIN) ) ri += cr * rminus * elem[i-1];
			if ( ! (bcvar & SURF_XMAX) ) ri += cr * rplus  * elem[i+1];   
			if ( ! (bcvar & SURF_YMIN) ) ri += cz * zminus * elem[i-xsize];
			if ( ! (bcvar & SURF_YMAX) ) ri += cz * zplus  * elem[i+xsize];   

			diag[n]  = 1 - (Lin ? Lin[i] : 0.0) + c * zcenter;
			sup[n]   = c * zplus; sub[n] = c * zminus;
			right[n] = ri; 

			i += xsize;
			if (++iy == ysize) break;

		} // while bcvar > 0 and iy < ysize
 
	  if (n)
		{
		  //for  (int t = 1; t <= n; t++) fprintf(stderr,"%g %g %g %g\n",sub[t],diag[t],sup[t],right[t]);
		tridiag(n);
		for (int ii = 1, ind = istart; ii<=n; ii++, ind += xsize) Result[ind] = right[ii];
		n = 0;
		}

	  if (iy < ysize)
		while ( !(ptype[i] & VARY_MASK) )  
		   {
		   i += xsize;
		   if ( ++iy == ysize) break; 
		   }

	  if (iy == ysize) {
		iy = 0;
		i = ++ix;
	  }

  } while (ix < xsize);

}

/***************************************************************************************
                                    1 D    R U N
                                    ^^^^^^^^^^^^

 Find vector Result_new (to be stored in array Result) satisfying

 (1 + c * Dx^2 - Lin) Result_new = cc * Dx^2 This + Result 

 | diag[1] sup[1]   0      0 :                    |                | right[1] |
 | sub[2]  diag[2] sup[2]  0 :            0       |                |          |
 |...........................:....................| x Result_new = |          |
 |                           : diag[n-1] sup[n-1] |                |          |
 |            0              : sub[n]    diag[n]  |                | right[n] |

 Operator Dx^2:   Dx^2 f[i] == minus f[i-1] + center f[i] + plus f[i+1] + coef
 off the boundary: minus = dxminus[i], center = dxcenter[i], plus = dxplus[i]

 sub[i]   = c minus
 diag[i]  = 1 + c center - Lin[i]
 sup[i]   = c plus
 right[i] = Result[i] + cc ( elem[i] center + elem[i-1] minus + elem[i+1] plus ) + (cc - c) coef
                                       
***************************************************************************************/

void FieldObj::Run1D( double c, double cc, double *Lin, double *Result)
{
long     i = 0, istart = 0, bcvar;
int      n = 0;
double   ri, plus, minus, center, coef;

if ( D <= 0.0 )  {
  for (i = 0; i < size; i ++)
    if (ptype[i] & VARY_MASK) Result[i] /= (1.0 -  (Lin ? Lin[i] : 0.0) );
  return;
  }

do  {
  while ( (bcvar = ptype[i]) & VARY_MASK ) {
    n++;
    if ( bcvar & SURF_XMIN )   istart = i; 

    nablaX(i, i, bcvar, minus, center, plus, coef);

    ri = Result[i] + (cc - c) * coef + cc * center * elem[i];

    if ( ! (bcvar & SURF_XMIN) ) ri += cc * minus * elem[i-1];
    if ( ! (bcvar & SURF_XMAX) ) ri += cc * plus * elem[i+1]; 

    diag[n] = 1 - (Lin ? Lin[i] : 0.0) + c * center ;
    sup[n] = c * plus;  sub[n] = c * minus;
    right[n] = ri;
    if (++i == xsize) break;

    } // while bcvar > 0 and i < xsize
 
  if (n)  {
    tridiag(n);
    for (int ii = 1, ind = istart; ii<=n; ii++, ind ++) Result[ind] = right[ii];
    n = 0;
  }

  if (i < xsize)
    while ( !(ptype[i] & VARY_MASK) )   {
       if ( ++i == xsize) break; 
    }

  } while (i < xsize);

}

/***************************************************************************************
 
  Laplace f[i]= dplus[i] f[i+1] + dminus[i] f[i-1] - (dplus[i] + dminus[i]) f[i] 
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Boundary at i=0:  f[-1] = f[0] - dx (a - 0.5 b (f[0] + f[-1]) )  where a = bc_const = bc_coef + bgr b, 
                                                                        b = bc_lin + bc_pump * elem[i]^(bc_pow-1) / (bc_Kn + elem[i]^bc_pow)
                   f[-1] = f[0] - dx a + 0.5 b dx (f[0] + f[-1])  
                   f[-1](1 - b dx / 2) = f[0] (1 + b dx / 2) - a dx

 dplus[0] f[1] + dminus[0] (f[0]*(1 + b dx / 2) - a dx) / (1 - b dx / 2)  - (dplus[0] + dminus[0]) f[0] 
  
 dplus[0] f[1] +  (dminus[0]*(1 + b dx / 2)/(1 - b dx / 2) - dplus[0] - dminus[0]) f[0] - a dx dminus[0] / (1 - b dx / 2)
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
 dplus[0] f[1] +  (dminus[0]*b*dx / (1 - b dx / 2) - dplus[0]) f[0] - a dx dminus[0] / (1 - b dx / 2)
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Dirichlet boundary at i = 0: linear intepolation (averaging): f[-1] + f[0] = 2 (bc_coef + bgr)

       dplus[0] f[1] + dminus[0] (2 (bc_coef + bgr) - f[0])  - (dplus[0] + dminus[0]) f[0] 
  
 dplus[0] f[1] - (2 dminus[0] + dplus[0]) f[0] + 2 dminus[0] (bc_coef + bgr)
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 Boundary at i=N:  f[N+1] = f[N] - dx (a - 0.5 b (f[N+1] + f[N]) )  where a = bc_const = bc_coef + bgr b, 
                                                                          b = bc_lin + bc_pump * elem[i]^(bc_pow-1) / (bc_Kn + elem[i]^bc_pow)

 dminus[N] f[N-1] +  (dplus[N]*b*dx / (1 - b dx / 2) - dminus[N]) f[N] - a dx dplus[N] / (1 - b dx / 2)
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
****************************************************************************************************************************/


//**********************************************************************************************

void FieldObj::nablaX(long i, int ix, long bcvar, double &minus, double &center, double &plus, double &coef)
{
    int bc;
    double bdx2;
    double Di = Diff[i];
    double Dm = 0.5 * (Di + Diff[i - (ix ? 1 : 0)] );
    double Dp = 0.5 * (Di + Diff[i + (ix + 1 < xsize ? 1 : 0)] );

    center = - (plus = Dp * dxplus[ix]) - (minus = Dm * dxminus[ix]);
    coef = 0.0;

    if ( bcvar & SURF_XMIN ) {
      if ( bc_deriv[bc = (bcvar >> XSHIFT) & BC_ID_MASK ] == 0) {  // Dirichlet
	     coef = 2 * minus * (bc_coef[bc] + bgr);
		 center -= minus;
	  } else {
	    bdx2 = 0.5 * xgrid[ix] * (bc_lin[bc] + bc_pump[bc] * pow(elem[i], bc_pow[bc]-1) / (bc_Kn[bc] + pow(elem[i], bc_pow[bc]) ) );
        center = -plus + minus * 2 * bdx2 / (1 - bdx2);
		coef   = -minus * xgrid[ix] * bc_const[bc] / (1 - bdx2);
      }
      minus = 0;
    }
    else if ( bcvar & SURF_XMAX ) {
      if ( bc_deriv[bc = (bcvar >> XSHIFT) & BC_ID_MASK ] == 0) {  // Dirichlet
	     coef = 2 * plus * (bc_coef[bc] + bgr);
		 center -= plus;
	  } else {
 	    bdx2 = 0.5 * xgrid[ix+1] * (bc_lin[bc] + bc_pump[bc] * pow(elem[i], bc_pow[bc]-1) / (bc_Kn[bc] + pow(elem[i], bc_pow[bc]) ) );
		center = -minus + plus * 2 * bdx2 / (1 - bdx2);
		coef   = -plus * xgrid[ix+1] * bc_const[bc] / (1 - bdx2);
      }
      plus = 0;
    }
}

//******************* calculating D_y^2  *************************

void FieldObj::nablaY(long i, int iy, long bcvar, double &minus, double &center, double &plus, double &coef, double h)
{
    int bc;
    double bdx2;
    double Di = Diff[i];
    double Dm = 0.5 * (Di + Diff[i - (iy > 0 ? xsize : 0)] );
    double Dp = 0.5 * (Di + Diff[i + (iy + 1 < ysize ? xsize : 0)] );

    center = - (plus = Dp * h * h * dyplus[iy]) - (minus = Dm * h * h * dyminus[iy]);
    coef = 0.0;

     if ( bcvar & SURF_YMIN ) {
      if ( bc_deriv[bc = (bcvar >> YSHIFT) & BC_ID_MASK ] == 0) {  // Dirichlet
	     coef = 2 * minus * (bc_coef[bc] + bgr);
		 center -= minus;
	  } else {
	    bdx2 = 0.5 * h * ygrid[iy] * (bc_lin[bc] + bc_pump[bc] * pow(elem[i], bc_pow[bc]-1) / (bc_Kn[bc] + pow(elem[i], bc_pow[bc]) ) );
        center = -plus + minus * 2 * bdx2 / (1 - bdx2);
		coef   = -minus * h * ygrid[iy] * bc_const[bc] / (1 - bdx2);
      }
      minus = 0;
    }
    else if ( bcvar & SURF_YMAX ) {
      if ( bc_deriv[bc = (bcvar >> YSHIFT) & BC_ID_MASK ] == 0) {  // Dirichlet
	     coef = 2 * plus * (bc_coef[bc] + bgr);
		 center -= plus;
	  } else {
 	    bdx2 = 0.5 * h * ygrid[iy+1] * (bc_lin[bc] + bc_pump[bc] * pow(elem[i], bc_pow[bc]-1) / (bc_Kn[bc] + pow(elem[i], bc_pow[bc]) ) );
        center = -minus + plus * 2 * bdx2 / (1 - bdx2);
		coef   = -plus * h * ygrid[iy+1] * bc_const[bc] / (1 - bdx2);
      }
      plus = 0;
    }
	
}

//******************* calculating D_z^2  *************************

void FieldObj::nablaZ(long i, int iz, long bcvar, double &minus, double &center, double &plus, double &coef, double h)
{
    int bc;
    double bdx2;
    double Di = Diff[i];
    double Dm = 0.5 * (Di + Diff[i - (iz > 0 ? xysize : 0)] );
    double Dp = 0.5 * (Di + Diff[i + (iz + 1 < zsize ? xysize : 0)] );

    center = - (plus = Dp * h * h * dzplus[iz]) - (minus = Dm * h * h * dzminus[iz]);
    coef = 0.0;

    if ( bcvar & SURF_ZMIN ) {
      if ( bc_deriv[bc = (bcvar >> ZSHIFT) & BC_ID_MASK ] == 0) {  // Dirichlet
	     coef = 2 * minus * (bc_coef[bc] + bgr);
		 center -= minus;
	  } else {
 	    bdx2 = 0.5 * h * zgrid[iz] * (bc_lin[bc] + bc_pump[bc] * pow(elem[i], bc_pow[bc]-1) / (bc_Kn[bc] + pow(elem[i], bc_pow[bc]) ) );
        center = -plus + minus * 2 * bdx2 / (1 - bdx2);
		coef   = -minus * h * zgrid[iz] * bc_const[bc] / (1 - bdx2);
      }
      minus = 0;
    }
    else if ( bcvar & SURF_ZMAX ) {
      if ( bc_deriv[bc = (bcvar >> ZSHIFT) & BC_ID_MASK ] == 0)  {  // Dirichlet
	     coef = 2 * plus * (bc_coef[bc] + bgr);
		 center -= plus;
	  } else {
 	    bdx2 = 0.5 * h * zgrid[iz+1] * (bc_lin[bc] + bc_pump[bc] * pow(elem[i], bc_pow[bc]-1) / (bc_Kn[bc] + pow(elem[i], bc_pow[bc]) ) );
		center = -minus + plus * 2 * bdx2 / (1 - bdx2);
		coef   = -plus * h * zgrid[iz+1] * bc_const[bc] / (1 - bdx2);
      }
      plus = 0;
    }
}


//******************* calculating D_y^2  *************************

/*

 void FieldObj::cleanBoundaries()
    { 
    for (long i = 0; i < size; i++)  
      if ( !(ptype[i] & VARY_MASK) )
        elem[i] = (ptype[i] & DIRICHLET) ? bc_coef[ptype[i] >> XSHIFT] : 0;
    }
*/
//***************************************************************

void FieldObj::print_ptype(long i)
{
long x = ptype[i];
if (x & _OUTSIDE_) fprintf(stderr,"[OUT] ");
else
  {
  if (x & VARY_MASK) fprintf(stderr,"[IN");
  if (x & SURF_MASK)
    {
    if (x & SURF_XMIN) fprintf(stderr," %s XMIN", bc_id[(x >> XSHIFT) & BC_ID_MASK]);
    if (x & SURF_XMAX) fprintf(stderr," %s XMAX", bc_id[(x >> XSHIFT) & BC_ID_MASK]);
    if (x & SURF_YMIN) fprintf(stderr," %s YMIN", bc_id[(x >> YSHIFT) & BC_ID_MASK]);
    if (x & SURF_YMAX) fprintf(stderr," %s YMAX", bc_id[(x >> YSHIFT) & BC_ID_MASK]);
    if (x & SURF_ZMIN) fprintf(stderr," %s ZMIN", bc_id[(x >> ZSHIFT) & BC_ID_MASK]);
    if (x & SURF_ZMAX) fprintf(stderr," %s ZMAX", bc_id[(x >> ZSHIFT) & BC_ID_MASK]);
    }
   fprintf(stderr,"] ");
  }
}


//****************************************************************************
//*                   C L A S S   B U F F E R   A R R A Y
//****************************************************************************

  BufferObj :: BufferObj(TokenString &params, const char *id, const char coopFlag) : FieldObj(params, id)   {

    long pos;
    coopType = coopFlag;
    params.get2_param(id, "D", &D, " undefined coefficient of diffusion"); // Assure it is defined

	if (coopFlag < 3) {
		double KD;
		int    fplus=0, fminus=0, fKD=0;

		if ( params.token_count(id, "KD") ) {
		  fKD = 1; 
		  params.get2_param(id, "KD",    &KD);  
		}

		if (params.token_count(id, "kplus")) {
		  fplus = 1;
		  pos = params.token_index(id, "kplus", 1) + 1;
		  if ( params.equal(pos, VAR_TOKEN) || params.equal(pos, ASSIGN_TOKEN) ) pos++;
		  kplus = new ExpressionObj(params, pos, "Bad expression for binding rate", 0, 0, &Time, "t");
		} 

		if (params.token_count(id, "kminus")) {
		  fminus = 1;
		  pos = params.token_index(id, "kminus", 1) + 1;
		  if ( params.equal(pos, VAR_TOKEN) || params.equal(pos, ASSIGN_TOKEN) ) pos++;
		  kminus = new ExpressionObj(params, pos, "Bad expression for unbinding rate", 0, 0, &Time, "t");
		}

		if ( fplus == 0) {
		  if (fKD && fminus) {
			kplus = new ExpressionObj(1);
			kplus->term_array[0].type = NUMBER_TYPE;
			kplus->term_array[0].val  = kminus->Evaluate() / KD;
		  }
		  else throw makeMessage("Undefined Ca binding properties for buffer \"%s\"", ID);
		}
		if ( fminus == 0) { 
		  if (fKD && fplus) {
			kminus = new ExpressionObj(1);
			kminus->term_array[0].type = NUMBER_TYPE;
			kminus->term_array[0].val  = kplus->Evaluate() * KD;
		  }
		  else throw makeMessage("Undefined Ca binding properties for buffer \"%s\"", ID);
		}
		
		if ( fplus && fminus && !fKD) KD = kminus->Evaluate() / kplus->Evaluate();

		if ( KD < 0 || kminus->Evaluate() < 0 || kplus->Evaluate() < 0 ||  fabs( kminus->Evaluate() - KD * kplus->Evaluate() ) > 1e-6 ) 
			throw makeMessage("Negative or incompatible Ca binding characteristics for buffer \"%s\": kplus=%g, kminus=%g, KD=%g", id, kplus->Evaluate(), kminus->Evaluate(), KD);
	    
		if (VERBOSE) fprintf(stderr,"    Kinetics for buffer %s: KD=%g uM, kminus=%g /ms, kplus=%g /(ms uM)\n", ID, KD, kminus->Evaluate(), kplus->Evaluate());

	} // if (coopFlag < 3) (end of binding rate definition block)

    total = new VectorObj(Size);
    double Total = 0;
    long ii, pMembrane = 0;
    
	pos = 0;
	params.token_count(id, "total", &pos);
    params.token3_count("buffer", ID, "membrane", &pMembrane);
    
	if (VERBOSE) fprintf(stderr, "    Total [%s] = ", ID);

    if ( !pos && coopFlag < 2 ) throw makeMessage("Total concentration not defined for buffer \"%s\"", id);

    if ( (D != 0.0 || pMembrane) && pos ) {
	   Total = ExpressionObj(params, pos + 1, makeMessage("Bad total concentration definition for buffer \"%s\"", id) ).Evaluate();
	   if (VERBOSE) fprintf(stderr,"%g \n", Total);
	}

    if ( D != 0.0 && pos )  *total = Total;

       if ( pMembrane ) {
		*total = 0.0;
	    int  bid, boxNum = Region->get_box_num();
		int  C_not_N = params.token2_count("membrane.definition", "concentration");  // <<==== New definition here
		if (VERBOSE) {
			if (C_not_N) fprintf(stderr, "    Membrane-bound total buffer concentration is expressed in uM units \n");
		    else         fprintf(stderr, "    Membrane-bound total buffer concentration is defined as molecules per um^2 \n"); 
		}
	
        if (D > 0.0) 
           params.errorMessage(params.token_index(ID, "D") + 1, 0, "Diffusion coefficient must be zero for membrane-bound buffer");

		double  Depth = 0;
        params.trail_pars( pMembrane, 'd', &Depth);

		int count = params.tokens_to_eol( pMembrane );

		if (count <= 2) {
			if (Depth == 0.0)  setMembraneBound(Total, C_not_N, boxNum);
						else   setMembraneLayer(Total, Depth,   boxNum);
		} 
		else 
			for(int kk = 2; kk <= count - 1; kk++) {
				bid = params.get_int( pMembrane + kk ) - 1;
				if (bid >= boxNum || bid < 0) 
					params.errorMessage( pMembrane + kk,  makeMessage("Box ID = %d exceeds the number of available boxes = %d", bid+1, boxNum) );
				if (Depth <= 0.0)  setMembraneBound(Total, C_not_N, bid);
				         	else   setMembraneLayer(Total, Depth,   bid);
			}
      }
    else if (D == 0.0 && pos ) 
            setFieldByFunction(params, pos + 1, total->elem, makeMessage("Bad expression for total concentration of %s",ID) );

	if (coopFlag) return;

   if (!params.token_count(id, "import") ) {
      double calcium = 0;
      calcium = params.get_param("Ca.bgr", 0, "Undefined background [Ca]");
      for (ii=0 ; ii < Size; ii++) {
         (*this)[ii] = (*total)[ii] * kminus->Evaluate() / (kminus->Evaluate() + kplus->Evaluate() * calcium);
         if ((*total)[ii] < 0 || (*this)[ii] < 0  )
             params.errorMessage( pos+1, makeMessage("Negative concentration for buffer %s at node %ld", ID, ii) );
      }
    }

    bgr = average();
	if ( VERBOSE ) fprintf(stderr,"\n### Free background [%s] = %g uM\n", ID, bgr);

    for (long ind = 0; ind < Size; ind++) if ( !(ptype[ind] & VARY_MASK) ) (*total)[ind] = (*this)[ind] = 0.0; 
    }

//****************************************************************************

  void  BufferObj :: setMembraneLayer(double Total, double depth, int boxid) 
    {
    long    ii, jj, Incr = 1;
    int     ix, iy, iz, j = 0, incr = 1, bid, point;  // j was uninitialized in previous version!!
    double  *coord = 0;
	int     boxNum = Region->get_box_num();

	if (VERBOSE) {
		if (boxid < boxNum) fprintf(stderr,"    Buffer %s is membrane-localized with depth %g um to box %d    \n", ID, depth, boxid+1); 
		               else fprintf(stderr,"    Buffer %s is membrane-localized with depth %g um to all boxes \n", ID, depth); 
	}
	
    for (ii=0 ; ii < Size; ii++) 
     if ( ptype[ii] & SURF_MASK ) { 

       Grid->split(ii, ix, iy, iz);
       point = point_type(ix, iy, iz);
	   bid = point >> 7;

	   if (boxid == bid || boxid >= boxNum) {

		   if ( point & (SURF_XMAX | SURF_XMIN) ) 
			  { Incr = 1;      coord = xcoord; j = ix; incr = (point & SURF_XMIN) ? 1 : -1; }  // dx = xgrid[j = ix]; 
		   if ( point & (SURF_YMAX | SURF_YMIN) ) 
			  { Incr = xsize;  coord = ycoord; j = iy; incr = (point & SURF_YMIN) ? 1 : -1; }  // dx = ygrid[j = iy]; 
		   if ( point & (SURF_ZMAX | SURF_ZMIN) ) 
			  { Incr = xysize; coord = zcoord; j = iz; incr = (point & SURF_ZMIN) ? 1 : -1; }  // dx = zgrid[j = iz]; 

		   Incr *= incr;
		   jj    = ii;

		   //fprintf(stderr,"ii=%ld (%d,%d,%d) point=%ld Incr=%ld incr=%d dx=%g\n",ii,ix,iy,iz,point,Incr,incr,dx); 

		   double coordinate = coord[j];
		   do {
			 if (ptype[jj] & VARY_MASK )  
			 total->elem[jj] = Total; 
			 jj += Incr; j += incr; 
		   } while ( fabs(coord[j] - coordinate) <= depth && jj < Size );
	   }
   }    
 }

//****************************************************************************

  void  BufferObj :: setMembraneBound(double Total, int C_not_N, int boxid) 
    {
    long    ii;
    int     ix, iy, iz, point, bid;

    double  Avogadro    = 602.2;  // Units of concentration per area are:
                               // uM * um = 10^(-13) mol /cm^2 = 10^(-21) mol / um^2 = 602.2 1 / um^2
	double  AreaDensity = Total / Avogadro;
    int     boxNum      = Region->get_box_num();

	if (VERBOSE) {
		if (boxid < boxNum) fprintf(stderr,"    Buffer %s is membrane-bound to box %d    \n", ID, boxid+1);
		               else fprintf(stderr,"    Buffer %s is membrane-bound to all boxes \n", ID);
	}

    for (ii=0 ; ii < Size; ii++) 
     if ( ptype[ii] & SURF_MASK ) { 

       Grid->split(ii, ix, iy, iz);
       point = point_type(ix, iy, iz);
	   bid   = point >> 7;

	   if (boxid == bid || boxid >= boxNum) {
		   if ( C_not_N )
			   total->elem[ii] = Total;
		   else  {
					 if ( point & (SURF_XMAX | SURF_XMIN) ) total->elem[ii] += AreaDensity / xgrid[ix];
				else if ( point & (SURF_YMAX | SURF_YMIN) ) total->elem[ii] += AreaDensity / ygrid[iy];
				else if ( point & (SURF_ZMAX | SURF_ZMIN) ) total->elem[ii] += AreaDensity / zgrid[iz];
		   }
	   }
	 }
  }    
//**************************************************************************

BufferArray::BufferArray(TokenString &params) {
   
	ID          = StrCpy("buffers");
    int coopNum = params.token2_count("buffer", "cooperative");
	nonCoopNum  = params.token_count("buffer") - coopNum;
	buf_num     = nonCoopNum + 3 * coopNum;
    array       = new BufferObj *[buf_num];

    if ( !buf_num )  return;
    int ind = 0;

    for (int i=0; i < params.token_count("buffer"); i++) {
      long p = params.token_index("buffer", i+1) + 1;
	  if ( params.equal(p, "cooperative") )  continue;
	  char *id = params.getVarName( p, "buffer" );
      if  ( params.token2_count("buffer", id) > 1 )
           params.errorMessage(p+1, makeMessage("Buffer \"%s\" defined more than once", id) );
      array[ind++] = new BufferObj(params, id);
      delete [] id;
    }  
	
    int index = nonCoopNum, iUB, iSB, iDB;

	for (int i=0; i < coopNum; i++) {

      long p = params.token2_index("buffer", "cooperative", i+1) + 1;
	  char *idUB = params.StrCpy( p );
	  char *idSB = StrCpy("Ca",  idUB);
	  char *idDB = StrCpy("Ca2", idUB);
      
	  iDB = 1 + (iSB = 1 + (iUB = index) );
	  index += 3;

      array[iUB] = new BufferObj(params, idUB, 1);
	  array[iSB] = new BufferObj(params, idSB, 2);
	  array[iDB] = new BufferObj(params, idDB, 3);

	  double C0 = 0.0;
      params.get_param("Ca.bgr", &C0);
	  double K1 = array[iUB]->kminus->Evaluate() / array[iUB]->kplus->Evaluate();
      double K2 = array[iSB]->kminus->Evaluate() / array[iSB]->kplus->Evaluate();
		  
      int totalSB = params.token_count(idSB, "total");
	  int totalDB = params.token_count(idDB, "total");

	  if ( ( !totalSB && totalDB) || ( totalSB && !totalDB ) )
		   params.errorMessage(p+1, makeMessage("Total concentrations of \"%s\" and \"%s\" should be either both defined or both undeclared (equilibrium)", idSB, idDB) );
	  else if (!totalSB) {
		   if (!params.token_count(idUB, "import") ) *array[iUB] = *array[iUB]->total * (1 / ( 1 + C0/K1 + C0*C0/(K1*K2) ));
		   if (!params.token_count(idSB, "import") ) *array[iSB] = *array[iUB]->total * (C0 / ( K1 + C0 + C0*C0/K2 ) );
		   if (!params.token_count(idDB, "import") ) *array[iDB] = *array[iUB]->total - *array[iSB] - *array[iUB];
	  } else {
	 	  *array[iUB] = *array[iUB]->total;
		  *array[iSB] = *array[iSB]->total;
		  *array[iDB] = *array[iDB]->total;
	  }

	  array[iUB]->bgr = array[iUB]->average();
	  array[iSB]->bgr = array[iSB]->average();
      array[iDB]->bgr = array[iDB]->average();

	  if ( VERBOSE ) {
		  fprintf(stderr,"\n### Initial average concentration [%s] = %g uM\n", idUB, array[iUB]->bgr );
		  fprintf(stderr,  "### Initial average concentration [%s] = %g uM\n", idSB, array[iSB]->bgr );
		  fprintf(stderr,  "### Initial average concentration [%s] = %g uM\n", idDB, array[iDB]->bgr );
	  }
      for (long ind = 0; ind < array[0]->Size; ind++)
		  if ( !(array[0]->ptype[ind] & VARY_MASK) ) {
 		   array[iUB]->total->elem[ind] = array[iSB]->total->elem[ind] = array[iDB]->total->elem[ind] = 0.0;
           array[iUB]->elem[ind]        = array[iSB]->elem[ind]        = array[iDB]->elem[ind]        = 0.0;
		  }

	  delete idUB; delete idSB; delete idDB;
    }
}

//****************************************************************************

void BufferArray::allocateCopy(const BufferArray &bufs)
    {
    if (bufs.ID) ID = StrCpy( bufs.ID );

    if ( (buf_num = bufs.buf_num) > 0 )
      {
      nonCoopNum = bufs.nonCoopNum;
      array = new BufferObj *[buf_num];
      for (int i = 0; i < buf_num; i++) array[i] = new BufferObj(*bufs.array[i]);
      }
    }

//****************************************************************************

double calcium_gain( FieldObj *Ca, BufferArray *Bufs) {
  double    x =  0.0;
  if (Ca)   x =  Ca->gain(); 
  if (Bufs) x += Bufs->gain(); 
  return x;
}
//****************************************************************************
