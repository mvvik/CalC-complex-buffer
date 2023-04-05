/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2019 Victor Matveev
 *
 *                            grid.cpp
 *
 *    Spatial grid definition class GridObj
 *    Boundary Condition definition array class BCarrayObj
 *    Grid stretching routines
 *    Difference scheme utility array initialization routines (initScheme)
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
#include "syntax.h"
#include "box.h"
#include "grid.h"

//*****************************************************************************
//                    G R I D   O P E R A T I O N S
//*****************************************************************************
   
 void GridObj::set_dimensions(RegionObj &reg, int nx, int ny, int nz)
 {
 nx = int(nx + 0.5);  if (nx == 0)  nx = 1;
 ny = int(ny + 0.5);  if (ny == 0)  ny = 1;
 nz = int(nz + 0.5);  if (nz == 0)  nz = 1;

 if (VERBOSE) { 
   fprintf(stderr,  "\n### Setting up the grid: [0..%d]", nx - 1);
   if ( DIMENSIONALITY > 1) fprintf(stderr, " x [0..%d]", ny - 1);
   if ( DIMENSIONALITY > 2) fprintf(stderr, " x [0..%d]", nz - 1);
   fprintf(stderr, "\n");
 }

 xsize = nx; ysize = ny; zsize = nz;

 xysize = nx * ny;
 Size   = nx * ny * nz;

 xgrid = new double[nx+1];  xcoord = new double[nx];
 ygrid = new double[ny+1];  ycoord = new double[ny];
 zgrid = new double[nz+1];  zcoord = new double[nz];

 dxplus = new double[nx]; dxminus = new double[nx];
 dyplus = new double[ny]; dyminus = new double[ny];
 dzplus = new double[nz]; dzminus = new double[nz];
 dvx    = new double[nx]; dvy     = new double[ny];  dvz = new double[nz];

 double dx = 1.0; double dy = 1.0; double dz = 1.0;
 xmin = reg.XMin;  xmax = reg.XMax;
 ymin = reg.YMin;  ymax = reg.YMax;
 zmin = reg.ZMin;  zmax = reg.ZMax;

 if (nx > 1) dx = (xmax - xmin) / double(nx);
 if (ny > 1) dy = (ymax - ymin) / double(ny);
 if (nz > 1) dz = (zmax - zmin) / double(nz);

 int i;
 for (i = 0; i < nx; i++)  
    { xgrid[i] = dx; xcoord[i] = xmin + (double(i) + 0.5) * dx; }
 for (i = 0; i < ny; i++)  
    { ygrid[i] = dy; ycoord[i] = ymin + (double(i) + 0.5) * dy; }
 for (i = 0; i < nz; i++)  
    { zgrid[i] = dz; zcoord[i] = zmin + (double(i) + 0.5) * dz; }

 xgrid[nx] = dx; ygrid[ny] = dy; zgrid[nz] = dz;

 xmax = reg.XMax = reg.XMin + double(nx) * dx; 
 ymax = reg.YMax = reg.YMin + double(ny) * dy; 
 zmax = reg.ZMax = reg.ZMin + double(nz) * dz; 
 }

//*****************************************************************************
   
 void GridObj::set_dimensions(RegionObj &region, TokenString &params)
 {
 int  nx=0, ny=0, nz=0;
 long pos;

 int gridNum = params.token_count("grid", &pos);
 if (gridNum == 0) 
   globalError( StrCpy("Forgot to define the grid!\n    Use \"grid nx ny nz\"") );

 if ( DIMENSIONALITY == 1)
   nx = params.get_int(pos + 1);
 else  if ( DIMENSIONALITY == 2)
   params.trail_pars(pos, 'i', &nx, 'i', &ny );
 else
   params.trail_pars(pos, 'i', &nx, 'i', &ny, 'i', &nz);


 if ( DIMENSIONALITY > 2) { 
   if (nz <= 0) params.errorMessage( pos + 3, 0, "Check the 3rd grid dimension"); }
 if ( DIMENSIONALITY > 1) { 
   if (ny <= 0) params.errorMessage( pos + 2, 0, "Check the 2nd grid dimension"); } 
   if (nx <= 0) params.errorMessage( pos + 1, 0, "Check the 1st grid dimension");

 set_dimensions(region, nx, ny, nz);
 }

//*****************************************************************************
   
int GridObj::location_to_index(double *coord, double *grid, int size, double x)
{
int i = 0;

if ( x < coord[0]-0.51*grid[0] || x > coord[size-1]+0.51*grid[size] ) 
	throw makeMessage("In GridObj::location_to_index: %g outside of bounds [%g, %g]", x, xmin, xmax);

while (x > coord[i] + 0.5 * grid[i+1] && i < (size - 1) ) i++;

return i;
}

//*****************************************************************************

long GridObj::location_to_index(double x, double y, double z)  // Only called from fieldObj::location_to_index()
{
	int ix = 0, iy = 0, iz = 0;

                         ix = location_to_index(xcoord, xgrid, xsize, x);
 if (DIMENSIONALITY > 1) iy = location_to_index(ycoord, ygrid, ysize, y);
 if (DIMENSIONALITY > 2) iz = location_to_index(zcoord, zgrid, zsize, z);

return ix + xsize * iy + xysize * iz;
}

//*****************************************************************************

void GridObj::split(long ind, int &ix, int &iy, int &iz)
    {
    iz = ind / xysize;
    iy = (ind -= iz * xysize) / xsize;
    ix = ind - iy * xsize; 
    }


//*****************************************************************************

void GridObj::grid_stretch(TokenString &params)
{
  char    dir[MAX_TOKEN_LENGTH];  // longest is "theta"
  double  x1, x2;
  double  f = 1.05;
  int     n = params.token_count("stretch");

  if (n == 0) return;

  params.get_param("stretch.factor", &f);

  if (f < 1.0 || f > 2.0) 
    params.errorMessage( params.token_index("stretch.factor") + 2, 
                         makeMessage("Stretch factor value = %g, should be between 1.0 and 2.0", f));

  for (int i = 1; i <= n; i++)
    {
    params.trail_pars("stretch", i, 'S', dir, 'd', &x1, 'd', &x2, 'E');  
    try { grid_stretch(dir, f, x1, x2); }
    catch (char *str) { params.errorMessage( params.token_index("stretch", i) + 1, str); }
    }
}


/*****************************************************************************

 Non-uniform grid in [xmin, x1] and [x2, xmax];
 Stretch intervals [0, n1] and [n2, n]; 
 dx = unstretched grid distance within [n1, n2]:

     0             n1    n2                    N
  :  |    |   |  | ||||||| |  |   |     |      |   :
xmin               x1    x2                       xmax
  <---    X1   ---><- X -><---       X2         --->

 1) X = dx (n2 - n1)  or  X = dx (n2 - n1 + 1/2)

 2) X1 = x1 - xmin = dx f + dx f^2 + ... + 3/2 * dx f^n1
 
        = f dx [ (f^n1 - 1) / (f - 1) + f^(n1-1)/2 ]
 
        = f dx [ 3f^n1/2 - f^n1 / (2f) - 1 ] / (f - 1)
 
        = dx ( f^n1 (3f - 1) - 2f ) / (2f - 2)
 
 3) X2 = xmax - x2 = dx f + dx f^2 + ...  + 3/2 * dx f^(N - n2)
 
        = dx ( f^(N-n2) (3f - 1) - 2f ) / (2f - 2)

                  SOLUTION:
                  ^^^^^^^^
 I.   From (3): f^(N - n2) = [ 2 X2 (f-1) / dx + 2 f ] / (3f - 1)

 II.  From (1): f^n1 = f^(N - X/dx) / f^(N - n2) = f^(N - X/dx) (3f - 1) / [ 2 X2 (f-1) / dx + 2 f ]

 III. Substitute into (2):  X1 = dx ( f^(N - X/dx) (3f - 1) / [ 2 X2 (f-1) / dx + 2f ] (3f - 1) - 2f ) / (2f - 2)
                               = dx ( f^(N - X/dx) (3f - 1)^2 / [ 4 X2 (f-1) / dx + 4f ] - f ) / (f - 1)

 IV.  Solve III for dx, using bisection

 V.   Get n1 from (2):  n1 = log( (2 X1 (f - 1) / dx + 2f) / (3f - 1) ) / log f

 VI.  Get n2 from (1):  n2 = X / dx + n1  <== If X1=0 or X2=0: X = dx(n2 - n1 + 0.5)

*****************************************************************************/
 
void GridObj::grid_stretch(const char *dir, double factor, double x1, double x2)
{
  int       N;
  double    *grid, *coord;
  double    eps = 1.0e-14, tiny = 1.0e-15, cmin, cmax, sfn;

  if (fabs(factor - 1) < 1e-8 ) return;

       if (::equal(dir, LABEL_DIM1) ) { grid = xgrid; coord = xcoord; N = xsize - 1; cmin = xmin;  cmax = xmax; }
  else if (::equal(dir, LABEL_DIM2) ) { grid = ygrid; coord = ycoord; N = ysize - 1; cmin = ymin;  cmax = ymax; }
  else if (::equal(dir, LABEL_DIM3) ) { grid = zgrid; coord = zcoord; N = zsize - 1; cmin = zmin;  cmax = zmax; }
  else throw makeMessage("Invalid grid stretch direction %s", dir);

  if (x1 < cmin) x1 = cmin;  
  if (x2 < cmin) x2 = cmin;  
  if (x1 > cmax) x1 = cmax;
  if (x2 > cmax) x2 = cmax; 

  double X1      = x1   - cmin;
  double X2      = cmax - x2;
  double X       = x2   - x1;
  double crange  = cmax - cmin;
  double correct = 0;

  if (x2 < x1) 
    throw makeMessage("3rd argument of the \"stretch\" command (%g) should be greater than the 2nd argument (%g)",x2,x1); 

  if (VERBOSE) 
    fprintf(stderr,"  # Stretching grid with factor=%g along %s-axis leaving [%.4g,%.4g] uniform\n", factor, dir, x1, x2);

  if ( X / crange < eps ) {
	  if ( X1 / crange < eps )  { grid_stretch(factor, 0, N, cmin, cmax, grid, coord); return; } // use simpler explicit method
	  if ( X2 / crange < eps )  { grid_stretch(factor, N, N, cmin, cmax, grid, coord); return; } // use simpler explicit method
	  correct = 0.5;
	  }

  //*************************************************************************************
  
  double logf = log(factor);
  double log23f = 2*log(3*factor - 1);
  double f1 = factor - 1;
  double f2 = 2*factor; 

  double dx0 = X / double(N + 1) + tiny;         // smallest possible dx
  double dx1 = crange / double(N + 1) - tiny;    // largest possible dx
  double y0 = (X / dx0 - N - correct) * logf - log23f + log( (2*X1*f1 / dx0 + f2 )*( 2*X2*f1 / dx0 + f2 ) ); 
  // double y1 = (X / dx1 - N - correct) * logf - log23f + log( (2*X1*f1 / dx1 + f2 )*( 2*X2*f1 / dx1 + f2 ) ); 
  double dxn, yn;

  do     // Bisection begins
    {
    dxn = 0.5 * (dx0 + dx1);
	yn = (X / dxn - N - correct) * logf - log23f + log( (2*X1*f1 / dxn + f2 )*( 2*X2*f1 / dxn + f2 ) ); 
    if (yn * y0 > 0) { dx0 = dxn; y0 = yn; }
                else { dx1 = dxn; } // y1 = yn; }
    }
  while (fabs( (dx0 - dx1) / (dx0 + dx1) ) > eps );

  int n1 = int( log( (2*X1*f1/dxn + f2)/(3*factor-1) ) / logf + 0.5 );   // step V.
  int n2 = int( log( (2*X1*f1/dxn + f2)/(3*factor-1) ) / logf + X/dxn + 0.5 );   // step VI.  
  
  if (VERBOSE>3) 
     fprintf(stderr,"  # n1=%d n2=%d dx=%g eps=%.5e\n", n1, n2, dxn, yn);

  //*************************************************************************************
  
  if ((n1 == n2) && (X / crange > eps) && VERBOSE)  
     fprintf(stderr, "\n>>> Grid stretch warning: the uniform interval too small to be resolved\n");

  double dx = dxn;
  if ( n2 > n1 ) {
	if ( (X1 / crange < tiny) || (X2 / crange < tiny) ) dx = X / (n2 - n1 + 0.5);
	                                              else  dx = X / (n2 - n1);
  }
  for (int i = n1 + 1; i <= n2; i++)  grid[i] = dx; 

  double dxlocal;

  if (X1 > 0) {
     sfn = stretchFactorCorrection(X1, dx, n1, dxlocal);
	 if (n1 > 0) for (int i = n1; i >= 1;  i--)  grid[i] = dxlocal * exp( log(sfn) * (n1 + 1 - i) );  // n1 => 0  
  } 
  
  if (X2 > 0) {
	 sfn = stretchFactorCorrection(X2, dx, N-n2, dxlocal);
	 if (n2 < N) for (int i = n2 + 1; i <= N;  i++)  grid[i] = dxlocal * exp( log(sfn) * (i - n2) );  // n2 => N  
  }
 
  grid[0]   = grid[1];
  grid[N+1] = grid[N];
  
  double range = 0.5 * (grid[0] + grid[N+1]);
  for (int i = 1; i <= N; i++)  range += grid[i]; 

  double gridCheck = crange / range;
  if (VERBOSE) fprintf(stderr, "    Grid check error %.2g\n", fabs(gridCheck - 1.0) );
  
  coord[0] = cmin + 0.5 * grid[0];
  for (int i = 1; i <= N; i++)   coord[i] = coord[i-1] + grid[i];

  double stretchOrder = fabs( log(grid[n2]/grid[0]) / log(10.0) );

  if ( stretchOrder > 4.0 ) 
	  if (VERBOSE)
		  fprintf(stderr, "%s\n", makeMessage("*** Warning: largest to smallest interval ratio spans %.1g orders of magnitude", stretchOrder));

  if (VERBOSE) {
    fprintf(stderr,"    Uniform %s-interval = [%.6g, %.6g] (%d..%d)\n", 
          dir, coord[n1], coord[n2], n1, n2);
    fprintf(stderr,"    Shortest dx = %.4g; dx on the boundaries: %.4g; %.4g\n", 
          grid[n2], grid[0], grid[N]);
  }

}
//*************************************************************************************
//      0    1   2    n1    n2                    N
//   :  |    |   |  | ||||||| |  |   |     |      |   :
// dx[0]dx[1]dx[2]...                                xmax
//
//*************************************************************************************
//   X1  = dx ( f^n1 (3f - 1) - 2f ) / (2f - 2)
//  2 X1 (f - 1) / dx = f^n1 (3f - 1) - 2f
//*************************************************************************************

double GridObj::stretchFactorCorrection(double XX, double dx, int NN, double &dxlocal)
{
	  double  eps = 1.0e-14, fmax = 2.0;
	  double  y0, yn;
	  double  dd = 2 * XX / dx;
	  double  sf0 = 1.00000001;           // smallest possible stretch factor
	  double  sf1 = fmax;                 // largest possible stretch factor
	  double  sfn;
	  dxlocal = dx;

	  y0 = exp(NN * log(sf0) ) * (3*sf0 - 1) -  2*sf0 - dd * (sf0 - 1);
	  // y1 = exp(NN * log(sf1) ) * (3*sf1 - 1) -  2*sf1 - dd * (sf1 - 1);

	  do     // Bisection begins
		{
		sfn = 0.5 * (sf0 + sf1);
		yn = exp(NN * log(sfn) ) * (3*sfn - 1) -  2*sfn - dd * (sfn - 1);
		if (yn * y0 > 0) { sf0 = sfn; y0 = yn; }
		  else { sf1 = sfn; } // y1 = yn; }
		}
	  while (fabs( (sf0 - sf1) / (sf0 + sf1) ) > eps );

	  if ( sfn > (fmax - 0.01) || sfn < 1.0 ) {
		  fprintf( stderr, ">>> Grid stretch unsuccessful: grid size either too small or too large\n");
		  sfn = 1.0;
		  dxlocal = XX / (NN + 0.5);
	  }
	  if (VERBOSE) fprintf(stderr, "    Stretch-factor corrected to %g;  dx=%g\n", sfn, dxlocal);
    
	  return sfn;
}

  //********************************************************************************

void GridObj::grid_stretch(double factor, int n1, int N, double cmin, double cmax, double *grid, double *coord)
{
  if (fabs(factor - 1) < 1e-8 ) return;

  double dx = (cmax - cmin) / ( (1 + exp(double(N-1) * log(factor)) )/2 + (1 - exp(double(N) * log(factor))) / (1 - factor) );

  if (n1 == 0) {
	  grid[0] = (grid[1] = dx);
      for (int i = 2; i <= N; i++)  grid[i] = factor * grid[i - 1]; 
      grid[N+1] = grid[N];
  } else if (n1 == N) {
	  grid[N] = (grid[N+1] = dx);
      for (int i = N-1; i >= 1; i--) grid[i] = factor * grid[i + 1]; 
      grid[0] = grid[1];
  }
          
  double range = 0.5 * (grid[0] + grid[N+1]);
  for (int i = 1; i <= N; i++)  range += grid[i];
  double gridCheck = (cmax - cmin) / range;
  if (VERBOSE) fprintf(stderr, "    Grid check error = %.2g\n", fabs(gridCheck - 1.0) );
  
  coord[0] = cmin + 0.5 * grid[0];
  for (int i = 1; i <= N; i++)   coord[i] = coord[i-1] + grid[i];

  if (VERBOSE) {
	  	           fprintf(stderr,"    Exact grid stretch factor implementation \n");
	  if (n1 == 0) fprintf(stderr,"    Shortest dx = %.4g; largest dx = %.4g\n", grid[0], grid[N]);
	          else fprintf(stderr,"    Shortest dx = %.4g; largest dx = %.4g\n", grid[N], grid[0]);
  }
  
  double stretchOrder = fabs( log(grid[N]/grid[0]) / log(10.0) );
  if ( stretchOrder > 4.0 ) 
	  if (VERBOSE)
		  fprintf(stderr, "%s\n", makeMessage("*** Warning: largest to smallest interval ratio spans %.1g orders of magnitude", stretchOrder));
}

//*****************************************************************************

void GridObj::initScheme() {
int i;

       if (GEOMETRY1 > 0) initSchemePolar( GEOMETRY1 );  // Polar / Cylindrical / Spherical
  else for (i = 0; i < xsize; i++) {
      dxplus[i]  = 2.0 / xgrid[i + 1] / (xgrid[i + 1] + xgrid[i]);
      dxminus[i] = 2.0 / xgrid[i]     / (xgrid[i + 1] + xgrid[i]);
      dvx[i]     = 0.5 * (xgrid[i] + xgrid[i+1]);
      }

  if (DIMENSIONALITY == 1) return;

  if (GEOMETRY2 == 3) initSchemeTheta();
  else 
    for (i = 0; i < ysize; i++)    {
      dyplus[i]  = 2.0 / ygrid[i + 1] / (ygrid[i + 1] + ygrid[i]);
      dyminus[i] = 2.0 / ygrid[i]     / (ygrid[i + 1] + ygrid[i]);
      dvy[i]     = 0.5 * (ygrid[i] + ygrid[i+1]);
      }

  if (DIMENSIONALITY == 2) return;

  for (i = 0; i < zsize; i++)    {
      dzplus[i]  = 2.0 / zgrid[i + 1] / (zgrid[i + 1] + zgrid[i]);
      dzminus[i] = 2.0 / zgrid[i]     / (zgrid[i + 1] + zgrid[i]);
      dvz[i]     = 0.5 * (zgrid[i] + zgrid[i+1]);
  }

}

/***************************************************************************************
 
                          POLAR COORDINATES DIFFERENTIAL 
                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^                                  
 with Tortuosity:

                1  d        df
 Dr^2 f[i,j] == - --- ( D r -- ) == rminus[i] f[i-1] + rcenter[i] f[i] + rplus[i] f[i+1]
                r  dr       dr 
 where r == r[i]

  r[i+1/2] D[i+1/2] (f[i+1]-f[i]) / dr[i+1] - D[i-1/2] r[i-1/2] (f[i] - f[i-1]) / dr[i] 
= ------------------------------------------------------------------------------------- =
                           r[i]  ( dr[i] + dr[i+1] ) / 2

         2                      D[i+1/2] r[i+1/2]          D[i-1/2] r[i-1/2]         
---------------------- { f[i+1] ----------------- + f[i-1] ----------------- - f[i] (..sum..) } =
r[i] (dr[i] + dr[i+1])                dr[i+1]                     dr[i]            

----------------------------------

 Special case of Neumann boundary ( df/dx[i] + bc_lin f[i] = bc_coef ):
 Assume uniform grid at boundary for accuracy!!

  :   |      |    
vir  f[0]   f[1]  
tual

   Boundary at i=0:  f[-1] = f[0] - dr (a - 0.5 b (f[0] + f[-1]) )   where a = bc_coef + bgr b, b = bc_lin
                     f[-1] = f[0] - dr a + 0.5 b dx (f[0] + f[-1])  
                     f[-1](1 - b dr / 2) = f[0] (1 + b dr / 2) - a dr

  rplus[0] f[1] + rminus[0] (f[0]*(1 + b dr / 2) - a dx) / (1 - b dr / 2)  - (rplus[0] + rminus[0]) f[0] 
   
   rplus[0] f[1] +  (rminus[0]*b*dr / (1 - b dr / 2) - rplus[0]) f[0] - a dr rminus[0] / (1 - b dr / 2)
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 rminus[N] f[N-1] +  (rplus[N]*b*dr / (1 - b dr / 2) - rminus[N]) f[N] - a dr rplus[N] / (1 - b dr / 2)
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
NOTE: SAME AS CARTESIAN!

 Special case r[-1/2] = 0 (radial grid starts at origin):  f[-1]=f[0]

          2                        D[1/2] r[1/2]         
-------------------- (f[1] - f[0]) -------------:  rcenter = -rplus
r[0] (dr[0] + dr[1])                   dr[1]                


***************************************************************************************/

void GridObj::initSchemePolar(int dim) {

	int i;
	long double pi = 4.0 * atan((long double)1.0);
	double f;

	 switch(GEOMETRY) {     // put symmetric variable measures into this, 1st coordinate norm
							case SPHERICAL:     f = 4.0*pi; break;
							case DISC:          
							case CONICAL:       
							case CYLINDRICAL:   f = 2.0*pi; break;
							default:            f = 1.0; // 3D spherical & 3D cylindrical: d(phi) still there
	 }

	 for (i = 0; i < xsize; i++)  {
	   double r       = xcoord[i];
	   double drplus  = xgrid[i + 1];
	   double drminus = xgrid[i];
	   double rplus   = r + 0.5 * drplus;
	   double rminus  = r - 0.5 * drminus;

	   if (rminus <= 0.0) rminus = 0.0;
	   if (r <= 0.0) throw makeMessage(">> Error in initiShemePolar: r=%g not positive-definite", r);

	   dxplus[i]  = 2 * pow( rplus  / r, dim) / (drplus  * (drplus + drminus) );
	   dxminus[i] = 2 * pow( rminus / r, dim) / (drminus * (drplus + drminus) );

	   //dvx[i]  = f * (pow(rplus, dim) - pow(rminus, dim));
	   dvx[i] = f * pow(r, dim) * (drplus + drminus) / 2.0;
	   }

}


//*****************************************************************************

 void GridObj::initSchemeTheta()
 { 
 int i;

 for (i = 0; i < ysize; i++)  {
   double theta    = ycoord[i];
   double dplus    = ygrid[i + 1];
   double dminus   = ygrid[i];
   double sinplus  = sin(theta + 0.5 * dplus);
   double sinminus = sin(theta - 0.5 * dminus);
   double c        = 2 / ( sin(theta) * (dplus + dminus) );

   dyplus[i]  = c * sinplus / dplus;
   dyminus[i] = c * sinminus / dminus;

   dvy[i]  = sin(theta) * (dplus + dminus) / 2.0;
   }
 }
//*****************************************************************************

 void GridObj::kill_grid()
 { 
 delete [] xgrid;  delete [] ygrid;   delete [] zgrid; 
 delete [] xcoord; delete [] ycoord;  delete [] zcoord;

 delete [] dxplus; delete [] dxminus; 
 delete [] dyplus; delete [] dyminus; 
 delete [] dzplus; delete [] dzminus; 

 delete [] dvx;  delete [] dvy;   delete [] dvz;
 } 

//*****************************************************************************
//     B O U N D A R Y   C O N D I T I O N   D A T A   S T R U C T U R E S
//*****************************************************************************

 void  BCarrayObj::set_bc_types(TokenString &params)
   {
   bc_type_count = 0; 

   int    predefined = 3;
   double CaD = 0.2;
   double Ca0 = 0;

   params.get_param("Ca.bgr", &Ca0);
   params.get_param("Ca.D",   &CaD);

   bc_type_num = params.token_count("bc.define") + predefined;

   bc_deriv = new double[bc_type_num]; bc_lin   = new double[bc_type_num];
   bc_coef  = new double[bc_type_num]; bc_id    = new char *[bc_type_num];  bc_const = new double[bc_type_num];
   bc_pump  = new double[bc_type_num]; bc_Kn    = new double[bc_type_num];  bc_pow   = new double[bc_type_num]; 
   bc_pump2 = new double[bc_type_num]; bc_Kn2   = new double[bc_type_num];  bc_pow2  = new double[bc_type_num];

   char s[256];
 
   set_bc_type(1, 0, 0, 1, 0, 0, Ca0, CaD, "Noflux");
   set_bc_type(0, 1, 0, 0, 0, 0, Ca0, 0.0, "Dirichlet");
   set_bc_type(0, 1, 0, 0, 0, 0, Ca0, 0.0, "Bgr");

   for (int i = 0; i < bc_type_num - predefined; i++)
     {
     double lin = 0, deriv = 0, coef = 0, pump = 0, Kd = 0, Kdinv, power = 0, pump2 = 0, Kd2 = 0, power2 = 0;
	 int    nPars = params.tokens_to_eol(params.token_index("bc.define", i + 1) + 1);
	 try {
		 switch (nPars)  {
			case 2:  params.trail_pars("bc.define", i+1, 's', s, 'd', &coef);
					 set_bc_type(    0,   1, 0, 1, 1, coef, Ca0,   0, s);
					 break;
			case 3:  params.trail_pars("bc.define", i+1, 's', s, 'd', &deriv, 'd', &lin);
					 // deriv = fabs(deriv); lin = -fabs(lin);
					 set_bc_type(deriv, lin, 0, 1, 1,    0, Ca0, CaD, s);
					 break;
			case 4:  params.trail_pars("bc.define", i+1, 's', s, 'd', &deriv, 'd', &lin, 'd', &coef);  // Backward compatibility
                     set_bc_type(deriv, lin, 0, 1, 1,  coef, Ca0, 1,  s);                              // Note: D is set to 1
					 break;
			case 5:  params.trail_pars("bc.define", i+1, 's', s, 'd', &deriv, 'd', &lin, 'd', &coef, 'd', &Kdinv);  // Backward compatibility
                     if (Kdinv == 0)   set_bc_type(deriv, lin,         0, 1,       1,  coef, Ca0, 1, s);            // Note: D is set to 1
					 else if (Kdinv>0) set_bc_type(deriv, 0,   lin/Kdinv, 1, 1/Kdinv,  coef, Ca0, 1, s);
					 else throw makeMessage("Inverse surface pump affinity should be positive (1/Kd = %g < 0)", Kdinv);
					 break;
			case 7:
			case 6:  params.trail_pars("bc.define", i+1, 's', s, 'd', &deriv, 'd', &lin, 'd', &pump, 'd', &power, 'd', &Kd);
                     deriv = fabs(deriv); lin = -fabs(lin); pump = -fabs(pump);      // Signs automatically adjusted: assumed outward current
					 set_bc_type(deriv, lin, pump, power, Kd, coef, Ca0, CaD, s); 
					 break;
			case 8:  params.trail_pars("bc.define", i+1, 's', s, 'd', &deriv, 'd', &lin, 'd', &pump, 'd', &power, 'd', &Kd, 'd', &Ca0, 'd', &CaD);
                     set_bc_type(deriv, lin, pump, power, Kd, coef, Ca0, CaD, s); 
					 break;
			case 9:  params.trail_pars("bc.define", i+1, 's', s, 'd', &deriv, 'd', &lin, 'd', &pump, 'd', &power, 'd', &Kd, 'd', &pump2, 'd', &power2, 'd', &Kd2 );
					deriv = fabs(deriv); lin = -fabs(lin); 
					pump = -fabs(pump);  pump2 = -fabs(pump2);    // Signs automatically adjusted: assumed outward current
					set_bc_type(deriv, lin, pump, power, Kd, pump2, power2, Kd2, coef, Ca0, CaD, s);
					break;
			default: throw makeMessage("Wrong number of arguments in bc.define");
		 }
	 } catch(char *str) { params.errorMessage( params.token_index("bc.define", i + 1) + 1, str); }
  
     if (VERBOSE) {
		       fprintf(stderr, "\n### Defining boundary condition \"%s\": \n", s);
			   double a  = bc_deriv[predefined + i], b  = bc_lin[predefined + i],  d  = bc_coef[predefined + i]; 
			   double c  = bc_pump[predefined + i],  p  = bc_pow[predefined + i],  k  = bc_Kn[predefined + i];
			   double c2 = bc_pump2[predefined + i], p2 = bc_pow2[predefined + i], k2 = bc_Kn2[predefined + i];
			   if (a == 0.0) { 
				   fprintf(stderr, "     Dirichlet:  u - u.bgr = %g \n", d);
				   continue;
			   }
			   b = b * CaD; c = c * CaD; c2 = c2 * CaD; d = d * CaD;
			            fprintf(stderr, "     Flux(in -> out) = ");
			   if (b)   fprintf(stderr, "%g * (u - %g) ", -b, Ca0 );
			   if (c) { fprintf(stderr, "+ %g * ( u^%.1g / (u^%.1g + %g) - ",  -c, p,      p, k);
						fprintf(stderr, "%g^%.1g / (%g^%.1g + %g) ) ",        Ca0, p, Ca0, p, k); }
			   if (c2) {
				   fprintf(stderr, "+ %g * ( u^%.1g / (u^%.1g + %g) - ", -c2, p2, p2, k2);
				   fprintf(stderr, "%g^%.1g / (%g^%.1g + %g) ) ", Ca0, p2, Ca0, p2, k2);
			   }
			   if (!b && !c && !c2)  {
				   if ( !d ) fprintf(stderr, "0 (Neumann, zero flux)");
				   else fprintf(stderr, "%g (constant flux)", d);
			   }
			   fprintf(stderr, "\n");
			   
	 }
   }
 }

//*****************************************************************************
 //set_bc_type(deriv, lin, pump, power, Kd, coef, Ca0, CaD, s)

   void  BCarrayObj::set_bc_type(double a, double b, double c, double p, double k, double c2, double p2, double k2, double d, double bgr, double CaD, const char *id)
   {
   if ( fabs(a) > 1e-12 ) 
     { b = b / a / CaD; c = c / a / CaD; c2 = c2 / a / CaD;  d = d / a / CaD;  a = 1;  }
   else if ( fabs(b) > 1e-12 )
   {
	   a = 0.0; c = 0.0; c2 = 0.0; d /= b; b = 1.0;
   }  // if a = 0, then linear coefficient = 1: it's a Dirichlet b.c. 
   else throw makeMessage("One of the first two arguments to \"bc.define\" should be non-zero (a=%g, b=%g)", a, b);
   
   bc_deriv[bc_type_count] = a;   bc_lin [bc_type_count] = b;  bc_coef[bc_type_count] = d;
   bc_pump[bc_type_count]  = c;   bc_pow[bc_type_count]  = 1;  bc_Kn[bc_type_count]   = 1;
   bc_pump2[bc_type_count] = c2;  bc_pow2[bc_type_count] = 1;  bc_Kn2[bc_type_count]  = 1;

   bc_const[bc_type_count] = d + bgr * b;

   if ( fabs(c) > 0.0 || fabs(c2) > 0.0 )  {
	    if (p  < 1)         throw makeMessage("Pump non-linearity exponent should be greater than 1 (p = %g < 1)",  p );
	    if (p2 < 1)         throw makeMessage("Pump non-linearity exponent should be greater than 1 (p2 = %g < 1)", p2);
		if (k  <= 1e-6*bgr) throw makeMessage("Pump affinity should be non-negligible (Kd = %g too small)",  k );
		if (k2 <= 1e-6*bgr) throw makeMessage("Pump affinity should be non-negligible (Kd2 = %g too small)", k2);
		bc_pow  [bc_type_count] = p;  bc_Kn [bc_type_count] = pow(k,  p );
		bc_pow2[bc_type_count]  = p2; bc_Kn2[bc_type_count] = pow(k2, p2);
		bc_const[bc_type_count] = d + bgr * ( b + c * pow(bgr, p-1) / ( pow(bgr, p) + pow(k, p) ) + c2 * pow(bgr, p2 - 1) / (pow(bgr, p2) + pow(k2, p2)));
   }

   bc_id[bc_type_count] = StrCpy(id);
   bc_type_count ++;
   }

   //*****************************************************************************
 //set_bc_type(deriv, lin, pump, power, Kd, coef, Ca0, CaD, s)

   void  BCarrayObj::set_bc_type(double a, double b, double c, double p, double k, double d, double bgr, double CaD, const char* id)
   {
	   if (fabs(a) > 1e-12)
	   {
		   b = b / a / CaD; c = c / a / CaD; d = d / a / CaD;  a = 1;
	   }
	   else if (fabs(b) > 1e-12)
	   {
		   a = 0; c = 0; d /= b; b = 1;
	   }  // if a = 0, then linear coefficient = 1: it's a Dirichlet b.c. 
	   else throw makeMessage("One of the first two arguments to \"bc.define\" should be non-zero (a=%g, b=%g)", a, b);

	   bc_deriv[bc_type_count] = a;   bc_lin[bc_type_count] = b;  bc_coef[bc_type_count] = d;
	   bc_pump[bc_type_count]  = c;   bc_pow[bc_type_count] = 1;  bc_Kn[bc_type_count] = 1;
	   bc_pump2[bc_type_count] = 0.0; bc_pow2[bc_type_count] = 1;  bc_Kn2[bc_type_count] = 1;

	   bc_const[bc_type_count] = d + bgr * b;

	   if (fabs(c) > 0.0) {
		   if (p < 1)         throw makeMessage("Pump non-linearity exponent should be greater than 1 (p=%g < 1)", p);
		   if (k <= 1e-6 * bgr) throw makeMessage("Pump affinity should be non-negligible (Kd=%g too small)", k);
		   bc_pow[bc_type_count] = p;  bc_Kn[bc_type_count] = pow(k, p);
		   bc_const[bc_type_count] = d + bgr * (b + c * pow(bgr, p - 1) / (pow(bgr, p) + pow(k, p)));
	   }

	   bc_id[bc_type_count] = StrCpy(id);
	   bc_type_count++;
   }

//*****************************************************************************

   int  BCarrayObj::bcidtoint(const char *id)
   {
   for (int i=0; i<bc_type_num; i++)
     if ( strcmp(bc_id[i], id) == 0 ) return i ;
     
   throw makeMessage("Boundary condition \"%s\" is undefined",id);
   }

//*****************************************************************************

   void  BCarrayObj::set_bc_types(int n)
   {
   bc_type_num = n;
   bc_deriv = new double[n];  bc_lin = new double[n];
   bc_pow   = new double[n];  bc_Kn  = new double[n]; bc_pump  = new double[n]; 
   bc_pow2  = new double[n];  bc_Kn2 = new double[n]; bc_pump2 = new double[n];
   bc_id    = new char *[n];
   bc_const = new double[n]; 
   }

  //*****************************************************************************

   void  BCarrayObj::kill_bc()
   {
   for (int i = 0; i < bc_type_num; i++) delete bc_id[i];

   delete [] bc_deriv; delete [] bc_lin;  delete [] bc_coef;  delete [] bc_id;

   delete[] bc_pump;  delete[] bc_Kn;  delete[] bc_pow;
   delete[] bc_pump2; delete[] bc_Kn2; delete[] bc_pow2;

   bc_type_count = 0;  // reset the static counter in set_bc_type
   }



//*****************************************************************************
//                        S T A T U S      W I D G E T S
//*****************************************************************************
