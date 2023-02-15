/*****************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *                Copyright (C) 2001-2019 Victor Matveev
 *
 *                               interpol.cpp
 *
 *  Interpolation/concentration averaging object "InterpolObj" 
 *  (appear in script as "Ca[x,y,z]" or "Ca[]")
 *
 *  Array of interpolation objects class "InterpolArray"
 *  InterpolArray is a member of the "KineticObj" defined in gate.h/gate.cpp
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
#include <stdarg.h>
#include <math.h>
#include "syntax.h"
#include "vector.h"
#include "box.h"
#include "grid.h"
#include "field.h"
#include "interpol.h"

//**************************************************************************

double InterpolObj::Evaluate(double t)
{
     if (*tptr != newtime)   // field has been updated since last interpolation evaluation
     {
       oldtime = newtime;
       oldres = newres;
       newtime = *tptr;

	   if (average) newres = Field->average();
	   else {
		   newres = 0.0;
		   for (int i = 0; i < nPointers; i++) newres += *pointers[i] * factors[i];
       }
	 }
     if (newtime > oldtime) // temporal extrapolation using values at two different times  
       result = ( (newtime - t) * oldres + (t - oldtime) * newres ) / (newtime - oldtime);
     else 
		 result = newres;
     
     return result;
}

//**************************************************************************

void InterpolObj::reset()
   {
   if (average) return;

   result = 0.0;
   for (int i = 0; i < nPointers; i++) result += *pointers[i] * factors[i];

   oldres = newres = result;
   oldtime = newtime = *tptr;   
   }

//*******************************************************************************
//*******************************************************************************

InterpolObj::InterpolObj(TokenString &Param, long p0, FieldObj *FO)
   {
     double *xcoord = FieldObj::Grid->xcoord; int xsize = FieldObj::Grid->xsize;
     double *ycoord = FieldObj::Grid->ycoord; int ysize = FieldObj::Grid->ysize;
     double *zcoord = FieldObj::Grid->zcoord; int zsize = FieldObj::Grid->zsize;
     int     xysize = FieldObj::Grid->xysize;

   double  x=0.0, y=0.0, z=0.0;
   long    p = p0, pLast;        // Param[p0] = "[", Param[p0-1] = "FieldID"
   
   pointers = 0; factors = 0;

   average = false;   
   char xyz[128];
   
   Field = FO;
   double *field = (*FO)();
   tptr = &(FO->Time); 
 
   if ( Param.equal(p + 1, "]") )  {
			 strcpy(xyz, Param[p0-1]); strcat( xyz, "[]" );
			 average = true;
			 Evaluate(0.0);
			 pLast = p0;
   } else {
		   if ( DIMENSIONALITY > 2 ) {
			 x = ExpressionObj(Param, p + 1, "", 0, &pLast).Evaluate();
			 if ( !Param.equal(pLast + 1, ",") ) 
			   Param.errorMessage(pLast + 1, 0, "Bad coordinate specifier or a missing comma");
			 p = pLast + 1;
		   }

		   if ( DIMENSIONALITY > 1 ) {
			 y = ExpressionObj(Param, p + 1, "", 0, &pLast).Evaluate();
			 if ( !Param.equal(pLast + 1, ",") ) 
			   Param.errorMessage(pLast + 1, 0, "Bad coordinate specifier or a missing comma");
			 p = pLast + 1;
		   }

		   z = ExpressionObj(Param, p + 1, "", 0, &pLast).Evaluate();
		   if ( !Param.equal(pLast + 1, "]") ) 
			   Param.errorMessage(pLast + 1, 0, "Bad coordinate specifier or a missing closing \"]\"");
 
		   if      ( DIMENSIONALITY == 1 )  { x = z; y = z    = 0.0; }
		   else if ( DIMENSIONALITY == 2 )  { x = y; y = z; z = 0.0; }

		   if (VERBOSE) { 
				fprintf(stderr, "  Interpolation object: "); 
				for (int j = p0-1; j <= pLast+1; j++) fprintf(stderr, "%s", Param[j]);
				fprintf(stderr, "\n");
			}

		   try { FO->location_to_index(x, y, z, 1); }  // Location evaluated: trip an error if more than 1 grid point out
		   catch (char *s) {
			   Param.errorMessage(p0, s, makeMessage("%s\n InterpolObj::InterpolObj => Location [%g, %g, %g] is out of bounds", s, x, y, z) );
		   }
   }
 
   token_ptr_source = &(Param.token_ptr[p0-1]);   // restore at destruction time
   token_ptr_target = Param.token_ptr[p0-1];      //  
   Param.deleteToken(p0, pLast + 1);

   oldtime = newtime = *tptr;   

   if (average)  {  result = Field->average();
                    oldres = newres = result;
					ID = StrCpy( xyz );
					Param.token_ptr[p0-1] = ID;   // *** hijack the token pointer
                    return;
   }

   nPointers = (1 << DIMENSIONALITY);
   pointers = new double *[nPointers];
   factors  = new double  [nPointers];

   int    ix = 1,   iy = 1,   iz = 1; 
   double p1 = 0.0, p2 = 0.0, p3 = 0.0, x0 = 0.0, x1 = 0.0, y0 = 0.0, y1 = 0.0, z0 = 0.0, z1 = 0.0;

   // Bracket the point (*not* the closest node!)
							   while (x > xcoord[ix] && ix < xsize - 1) ix++;  x0 = xcoord[ix-1]; x1 = xcoord[ix]; 
	if (DIMENSIONALITY > 1)  { while (y > ycoord[iy] && iy < ysize - 1) iy++;  y0 = ycoord[iy-1]; y1 = ycoord[iy]; }
	if (DIMENSIONALITY > 2)  { while (z > zcoord[iz] && iz < zsize - 1) iz++;  z0 = zcoord[iz-1]; z1 = zcoord[iz];  }

                            p1 = ( (x - x0) / (x1 - x0) );
   if ( DIMENSIONALITY > 1) p2 = ( (y - y0) / (y1 - y0) );
   if ( DIMENSIONALITY > 2) p3 = ( (z - z0) / (z1 - z0) );

   if (p1 < 1e-6) p1 = 0.0; else if (p1 > 1.0-1e-6) p1 = 1.0;
   if (p2 < 1e-6) p2 = 0.0; else if (p2 > 1.0-1e-6) p2 = 1.0;
   if (p3 < 1e-6) p3 = 0.0; else if (p3 > 1.0-1e-6) p3 = 1.0;

   double factorSum = 0.0, Xave =0.0, Yave=0.0, Zave=0.0;
   int cnt = 0;

   for (int i = 0; i < nPointers; i++)
   {
	   double factor = 1.0;
	   int    dx=0, dy=0, dz=0;
	   long   ind;
							     dx =        i & 1;  factor *= ( dx ? (1 - p1) : p1 );  ind = ix - dx;
	   if (DIMENSIONALITY > 1) { dy = (i >> 1) & 1;  factor *= ( dy ? (1 - p2) : p2 );  ind += (iy - dy) * xsize;  }
	   if (DIMENSIONALITY > 2) { dz = (i >> 2) & 1;  factor *= ( dz ? (1 - p3) : p3 );  ind += (iz - dz) * xysize; }

	   if ( (factor > 0.0) && (FO->ptype[ind] & _INSIDE_) ) 
	   {
									 Xave += factor * xcoord[ix - dx];  
			if (DIMENSIONALITY > 1)  Yave += factor * ycoord[iy - dy];
			if (DIMENSIONALITY > 2)  Zave += factor * zcoord[iz - dz];  

		   if (VERBOSE > 3) fprintf(stderr, "    Use point (%g, %g, %g) [%ld], weight=%g\n", xcoord[ix - dx], ycoord[iy - dy], zcoord[iz - dz], ind, factor);
		   factors[cnt]  = factor;
		   pointers[cnt] = field + ind;
		   factorSum += factor;
		   cnt ++;
	   }
   }

   Xave /= factorSum;  Yave /= factorSum;  Zave /= factorSum;

   if ( (nPointers = cnt) == 0 ) Param.errorMessage(pLast, 0, "Can't interpolate: point out of bounds");
   else if (VERBOSE) 
		{ fprintf(stderr, "  # Interpolated with %d points, average position=(%g, %g, %g)\n", nPointers, Xave, Yave, Zave); 
          fflush(stderr); }

   switch ( DIMENSIONALITY ) { 
		   case 1:  snprintf(xyz, 127, "%s[%g]", FO->ID, x); 
					break;
		   case 2:  snprintf(xyz, 127, "%s[%g,%g]", FO->ID, x, y);
					break;
		   default: snprintf(xyz, 127, "%s[%g,%g,%g]", FO->ID, x, y, z);
		   }

   ID = StrCpy( xyz );
   Param.token_ptr[p0-1] = ID;   // *** hijack the token pointer

   result = 0.0;
   for (int i = 0; i < nPointers; i++) result += *pointers[i] * (factors[i] /= factorSum);
   
   oldres = newres = result;
   }
   
//*******************************************************************************

InterpolArray::InterpolArray(TokenString &Param, FieldObj *Ca, BufferArray *Bufs)
{
  interpol_num = 0;
  long p0;
  int i, j;

	if (Ca) interpol_num += Param.token2_count(Ca->ID, "[");

	if (Bufs)
		for (j = 0; j < Bufs->buf_num; j++)
			interpol_num += Param.token2_count(Bufs->array[j]->ID, "[");

	array = new InterpolObj *[interpol_num];
	int count = 0;

	if (!interpol_num) return; 

	if (VERBOSE) { fprintf(stderr,"\n### Interpolating: \n"); fflush(stderr); }

	if (Ca) 
	  while ( Param.token2_count(Ca->ID, "[") )
		{
		p0 = Param.token2_index(Ca->ID, "[", 1);
		array[count] = new InterpolObj(Param, p0, Ca );
		count ++;
		}

	if (Bufs)
	  for (j = 0; j < Bufs->buf_num; j++)
		while ( Param.token2_count(Bufs->array[j]->ID, "[") )
		  {
		  p0 = Param.token2_index(Bufs->array[j]->ID, "[", 1);
		  array[count] = new InterpolObj(Param, p0, Bufs->array[j] );
		  count ++;
		  }

	if (VERBOSE) fprintf(stderr," ### Compressing interpolation array: ");

	InterpolObj *temp;

	for (i = 0; i < interpol_num - 1; i++ )
	  for  (j = i + 1; j < interpol_num; j++ )
	  if ( equal(array[i]->ID, array[j]->ID) )
		{
		if (VERBOSE) { fprintf(stderr,"."); fflush(stderr); }
			interpol_num--;
		temp = array[interpol_num]; 
			array[interpol_num] = array[j];
			array[j] = temp;
		}

    if (VERBOSE) fprintf(stderr," done\n");

}

//*******************************************************************************

