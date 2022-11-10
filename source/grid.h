/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2019 Victor Matveev
 *
 *                            grid.h
 *
 *    Spatial grid definition class GridObj
 *    Boundary Condition definition array class BCarrayObj
 *    Grid stretching routines
 *    Difference scheme utility array initialization routines (initScheme)
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

#ifndef CALC_GRID_H_included
#define CALC_GRID_H_included

//*****************************************************************************
//                    G R I D   O P E R A T I O N S
//*****************************************************************************

class RegionObj;

class GridObj {

 public:

 int     xsize, ysize, zsize;
 long    Size, xysize;

 double  *xgrid, *ygrid, *zgrid;
 double  *xcoord, *ycoord, *zcoord;
 double  *dxplus, *dxminus;
 double  *dyplus, *dyminus;
 double  *dzplus, *dzminus;
 double  *dvx, *dvy, *dvz;
 double  xmin, xmax, ymin, ymax, zmin, zmax;

 GridObj(RegionObj &reg, TokenString &params) { 
      set_dimensions(reg, params);  // read grid dimensions from parameter file
      grid_stretch(params);         // make non-uniform grid (read parameters from file)
      // grid_clamp(params);        // clamp grid nodes to specified locations
      initScheme();                 // set up finite difference support arrays
 }

 ~GridObj() { kill_grid(); }

 void set_dimensions(RegionObj &region, int nx, int ny, int nz);
 void set_dimensions(RegionObj &region, TokenString &params);

 void exportGrid(RegionObj *region);
 void grid_stretch(const char *dir, double factor, double x1, double x2);
 void grid_stretch(double factor, int n1, int N, double cmin, double cmax, double *grid, double *coord);
 void grid_stretch(TokenString &params);
 // void grid_clamp(TokenString &params);
 double stretchFactorCorrection(double XX, double dx, int NN, double &dxlocal);
 
 void initScheme();
 void initSchemePolar(int power);
 void initSchemeTheta();

 void kill_grid();
 long location_to_index(double x, double y, double z);
 int  location_to_index(double *, double *, int, double);

 double nearest_location(double *coord, double *grid, int size, double c)
   { return coord[ location_to_index(coord, grid, size, c) ]; }

 void split(long, int &, int &, int &);

};

//*****************************************************************************
//     B O U N D A R Y   C O N D I T I O N   D A T A   S T R U C T U R E S
//*****************************************************************************

class BCarrayObj {

 public:
  
   double  *bc_deriv, *bc_lin, *bc_coef, *bc_pump, *bc_pow, *bc_Kn, *bc_const;
   char    **bc_id;
   int     bc_type_num;
   int     bc_type_count;

   BCarrayObj(TokenString &params) { bc_type_count = 0; set_bc_types(params); }
   ~BCarrayObj() { kill_bc(); }

   void  set_bc_types(TokenString &params);
   void  set_bc_types(int n);
   void  set_bc_type(double a, double b, double c, double p, double k, double d, double Ca0, double CaD, const char *id = "");
   int   bcidtoint(const char *);
   void  kill_bc();

};

//*****************************************************************************

#endif
