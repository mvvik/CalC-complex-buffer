/*****************************************************************************
 *
 *                        Calcium Calculator (CalC)
 *                 Copyright (C) 2001-2022 Victor Matveev
 *
 *                                vector.cpp
 *
 *                   A basic linear algebra vector class 
 *
 *****************************************************************************
 
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

 *****************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "PlatformSpecific.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>  // for compatibility with Visual C++
#include <string.h>
#include "vector.h"

// extern int xsize, ysize, zsize;

extern FILE *fopenAssure(const char *fname, const char *mode, const char *action, const char *id);
extern char *makeMessage(const char *fmt, ...);

//*******************************************************************************
//                       I M P L E M E N T A T I O N
//*******************************************************************************


VectorObj::VectorObj(const VectorObj &v)
{
  size = v.size;

  elem = new double[size];
  for (long i = 0; i < size; ++i)  elem[i] = v.elem[i];
}

//==============================================================================

VectorObj &VectorObj::operator=(const double *e)
{
	for (long i = 0; i < size; ++i) elem[i] = e[i];
	return *this;
}

VectorObj &VectorObj::operator=(const double val)
{
	for (long i = 0; i < size; ++i) elem[i] = val;
	return *this;
}

VectorObj &VectorObj::operator=(const VectorObj &v)
{
	for (long i = 0; i < size; ++i) elem[i] = v.elem[i];
	return *this;
}


//------------------------------------------------------------------------------


VectorObj *vabs(const VectorObj &v)
{
       VectorObj *vv = new VectorObj(v.size);
       for (long i = 0; i < v.size; ++i) (*vv).elem[i] = fabs(v.elem[i]);
       return vv;
}

//------------------------------------------------------------------------------


VectorObj operator-(const VectorObj &v)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i) vv.elem[i] = -v.elem[i];
       return vv;
}


VectorObj operator-(const VectorObj &v1, const VectorObj &v2)
{
       VectorObj v3(v1.size);
       for (long i = 0; i < v1.size; ++i)  v3.elem[i] = v1.elem[i]-v2.elem[i];
       return v3;
}

VectorObj operator-(const VectorObj &v, const double val)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i)  vv.elem[i] = v.elem[i] - val;
       return vv;
}

VectorObj operator-(const double val, const VectorObj &v)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i)  vv.elem[i] = val - v.elem[i];
       return vv;
}

void operator-=(VectorObj &v1, const VectorObj &v2)
{
       for (long i = 0; i < v1.size; ++i)  v1.elem[i] -= v2.elem[i];
}


void operator-=(VectorObj &v, const double val)
{
       for (long i = 0; i < v.size; ++i)  v.elem[i] -= val;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


VectorObj operator+(const VectorObj &v1, const VectorObj &v2)
{
       VectorObj v3(v1.size);
       for (long i = 0; i < v1.size; ++i)  v3.elem[i] = v1.elem[i] + v2.elem[i];
       return v3;
}

VectorObj operator+(const VectorObj &v, const double val)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i)  vv.elem[i] = v.elem[i] + val;
       return vv;
}

VectorObj operator+(const double val, const VectorObj &v)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i)  vv.elem[i] = v.elem[i] + val;
       return vv;
}

void operator+=(VectorObj &v1, const VectorObj &v2)
{
       for (long i = 0; i < v1.size; ++i)  v1.elem[i] += v2.elem[i];
       //       return v1;
}

void operator+=(VectorObj &v, const double val)
{
       for (long i = 0; i < v.size; ++i)  v.elem[i] += val;
       //return v;
}

// inline VectorObj& operator+(VectorObj &v)
// {
// return v;
// }

//*******************************************************************************

VectorObj operator%(const VectorObj &v1, const VectorObj &v2)
{
       VectorObj v3(v1.size);
       for (long i = 0; i < v1.size; ++i)  v3.elem[i] = v1.elem[i] * v2.elem[i];
       return v3;
}


double operator*(const VectorObj &v1, const VectorObj &v2)
{
       double sum = 0.0;
       for (long i = 0; i < v1.size; ++i)  sum += v1.elem[i] * v2.elem[i];
       return sum;
}

VectorObj operator*(const VectorObj &v, const double val)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i)  vv.elem[i] = v.elem[i] * val;
       return vv;
}

VectorObj operator*(const double val, const VectorObj &v)
{
       VectorObj vv(v.size);
       for (long i = 0; i < v.size; ++i)  vv.elem[i] = v.elem[i] * val;
       return vv;
}

VectorObj operator%=(VectorObj &v1, const VectorObj &v2)
{
       for (long i = 0; i < v1.size; ++i)  v1.elem[i] *= v2.elem[i];
       return v1;
}


VectorObj operator*=(VectorObj &v, const double val)
{
       for (long i = 0; i < v.size; ++i)  v.elem[i] *= val;
       return v;
}

/****************************************************************************/

void  VectorObj::print(int n) const
{
        char s1[20], s2[20];

        strcpy(s1,"%.");
        sprintf(s2,"%d",n);
        strcat(s1,s2);
        strcat(s1,"g ");

	for (long i = 0; i < size; ++i)  fprintf(stderr, s1, elem[i]);
	printf("\n");
}
/****************************************************************************/

double  VectorObj::norm() const
{
    double sum = 0.0;
	for (long i = 0; i < size; ++i) sum += elem[i]*elem[i]; 
	return sqrt(sum);
}

double  VectorObj::norm_L1() const
{
    double sum = 0.0;
	for (long i = 0; i < size; ++i) sum += elem[i]; 
	return sum;
}

//***********************************************************************************
//***********************************************************************************

double VectorObj::checkerBoardNorm1D()
{
   double  norm = 0.0;
   double  *e = elem + 1;     
   double  *e1 = e - 1, *e2 = e + 1;
  
   for (long i = 1; i < size - 1; i++)
    //if ( (ptype[i] & _V_INSIDE_) && !(ptype[i] & V_SURF_MASK) )
		norm += fabs( *(e1++) + *(e2++) - 2.0 * *(e++) );
   
   return norm;
}

//***********************************************************************************

double VectorObj::checkerBoardNorm2D(int S1)
{
   double  norm = 0.0;
   double  *e = elem + S1;     
   double  *e1 = e - 1, *e2 = e + 1, *e3 = e - S1, *e4 = e + S1;
  
   for (long i = S1; i < size - S1; i++)
    //if ( (ptype[i] & _V_INSIDE_) && !(ptype[i] & V_SURF_MASK) )
		norm += fabs( *(e1++) + *(e2++) + *(e3++) + *(e4++) - 4.0 * *(e++) );
   
   return norm;
}

//***********************************************************************************

double VectorObj::checkerBoardNorm3D(int S1, int S2)
{
   double  norm = 0.0;
   double  *e = elem + S2;     
   double  *e1 = e - 1, *e2 = e + 1, *e3 = e - S1, *e4 = e + S1, *e5 = e - S2, *e6 = e + S2;
  
   for (long i = S2; i < size - S2; i++)
    //if ( (ptype[i] & _V_INSIDE_) && !(ptype[i] & V_SURF_MASK) )
		norm += fabs( *(e1++) + *(e2++) + *(e3++) + *(e4++) + *(e5++) + *(e6++) - 6.0 * *(e++) );
   
   return norm;
}


//****************************************************************************/

double VectorObj::L_infty_norm() const
{
   double norm = 0.0, temp = 0.0;

   for (long i = 0; i < size; ++i)
     { 
     temp = fabs(elem[i]);
     if ( _isnan(temp) ) return temp;
     if (temp > norm) norm = temp;
     }

   return norm;
}

//****************************************************************************/

void VectorObj::Export(const char *filename)
{
  FILE *f = fopenAssure(filename, "wb", "Vector export", "");

  fwrite( (void *)(&size), sizeof(long),   1, f);
  fwrite( (void *)elem,    sizeof(double), size, f);

  fclose(f);

}

//****************************************************************************/

void VectorObj::import(const char *filename)
{
  FILE *f = fopenAssure(filename, "rb", "Vector import", "");

  long n;  
  fread( (void *)(&n), sizeof(long),   1, f);

  if (n != size) throw makeMessage("in Vector import: wrong number of elements (%ld vs %ld)\n",n,size);

  fread( (void *)elem, sizeof(double), size, f);

  fclose(f);

}

//****************************************************************************/
