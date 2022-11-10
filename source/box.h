/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2019 Victor Matveev
 *
 *                            box.h
 *
 *               Geometry definition classes:
 * 
 *               Generalized box class BoxObj
 *               Box as an obstacle: class ObstBoxObj
 *               Set of boxes and obstacles: RegionObj
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

#ifndef CALC_BOX_H_included
#define CALC_BOX_H_included

#define  SURF_XMIN   1
#define  SURF_XMAX   2
#define  SURF_YMIN   4
#define  SURF_YMAX   8
#define  SURF_ZMIN   16
#define  SURF_ZMAX   32
#define  SURF_MASK   63
#define  _INSIDE_    64
#define  INSIDE_MASK 127

//***************************************************************
//*                   C L A S S   B O X
//***************************************************************

class GridObj;

//***************************************************************

class VolumeObjClass
{
protected:

  double p1, p2, p3, p4, p5, p6;  // Any volume is defined by at most 6 parameteres
  
public:

  bool   *isInside;
  long   tokenPosition;

  static GridObj *Grid;
  static void bindToGrid(GridObj &grid) { Grid = &grid; }
  void   computeFormulas(TokenString &pars);

  char   volumeType;
  int    isObstacle;
  double xmin, xmax, ymin, ymax, zmin, zmax;


  VolumeObjClass() { isInside = 0; isObstacle = 0; tokenPosition = 0;  p1 = p2 = p3 = p4 = p5 = p6 = xmin = xmax = ymin = ymax = zmin = zmax = 0.0; };
   VolumeObjClass(TokenString &pars, long pos, int isObst);
   VolumeObjClass(const VolumeObjClass &v) { this->set(v); }
  ~VolumeObjClass() {  };

  int point_type(int, int, int);

  
  void setPars(int isObst, char vType, double pp1, double pp2, double pp3=0.0, double pp4=0.0, double pp5=0.0, double pp6=0.0)
    {volumeType = vType;  isObstacle = isObst;  p1 = pp1;  p2 = pp2;  p3 = pp3;  p4 = pp4;  p5 = pp5;  p6 = pp6; }

  void setSize(double x0, double x1, double y0, double y1, double z0, double z1) { xmin = x0;  xmax = x1;  ymin = y0;  ymax = y1;  zmin = z0;  zmax = z1; }

  void set (const VolumeObjClass &B) { setPars(B.isObstacle, B.volumeType, B.p1, B.p2, B.p3, B.p4, B.p5, B.p6);
                                       if ( !(B.isObstacle) ) setSize(B.xmin, B.xmax, B.ymin, B.ymax, B.zmin, B.zmax); 
									   isInside = B.isInside;
									   tokenPosition = B.tokenPosition;
  }

  VolumeObjClass &operator=(const VolumeObjClass &B) { set(B); return *this; }

  int inside(int, int, int);
};

//***************************************************************
//*                   C L A S S   R E G I O N
//***************************************************************

class RegionObj
{
protected: 

 double *xgrid, *ygrid, *zgrid;
 double *xcoord, *ycoord, *zcoord;

 int    xsize, ysize, zsize;
 long   xysize;

 int    encl_num, obst_num, total;
 
public:

 double XMax, XMin, YMax, YMin, ZMax, ZMin;

 static GridObj *Grid;
 static void bindToGrid(GridObj &grid) { Grid = &grid; VolumeObjClass::bindToGrid(grid); }

 VolumeObjClass    *enclosures;
 VolumeObjClass    *obstacles;
 
 int    get_box_num()      { return     total; }
 int    get_surface_num()  { return 2 * DIMENSIONALITY * total; }

 void   computeFormulas(TokenString &pars);

 // void initialize(int vNum, int oNum);
 // RegionObj(int encls, int obstcls) { initialize(encls, obstcls); };
 
 RegionObj(TokenString &);
 ~RegionObj();

 int  point_type(int, int, int, int fieldObstNum, class VolumeObjClass *OBO);
 // int  point_type(long,          int fieldObstNum, class VolumeObjClass *OBO);
};

void print_type(int x);

#endif
