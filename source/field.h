/*****************************************************************************
 *
 *                      Calcium Calculator (CalC)
 *                Copyright (C) 2001-2021 Victor Matveev
 *
 *                                field.h
 *
 *  Definitions of the main FieldObj concentration field class (derived from
 *  the VectorObj class), the daughter BufferObj buffer concentration field 
 *  class, and the compound buffer array class BufferArray.
 *
 *  FieldObj methods include all low-level ADI difference schemes routines.
 *  Highest-level numeric engine routines are defined in "calc.cpp"       
 *       
 *  FieldObj and BufferArray variables are member of the main SimulationObj
 *  object, defined in simulation.h / simulation.cpp
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

#ifndef CALC_FIELD_H_included
#define CALC_FIELD_H_included

#define SOURCE_EPS  1.0e-4
#define MIN_CURRENT 1.0e-12

const int  VARY_MASK = 127;
const int  REGION_MASK = 255;  // DIRICHLET + VARY_MASK
const int  _OUTSIDE_ = 256;

const int  BC_ID_MASK = 127;  // 5 bits for bc type id - at most 32 different surface types
const int  XSHIFT = 9;
const int  YSHIFT = 16;
const int  ZSHIFT = 23;

typedef int  (*FieldMethod) (class FieldObj &,  VectorObj &, class BufferArray &, double dt, double);
typedef int  (*BufferMethod)(class BufferObj &, VectorObj &, class FieldObj &,    double dt, double); 

double gaussian(double, double, double, int);
double square(double, double, double, int);

extern int GEOMETRY;
extern int DIMENSIONALITY;

double calcium_gain( FieldObj *Ca, BufferArray *Bufs);


//****************************************************************************
//*                         C L A S S   F I E L D
//****************************************************************************

class FieldObj: public VectorObj
{ 

protected:

  //static  double (*integrate)(double, double, double, int);

  static  double *dxplus, *dxminus, *dvx;
  static  double *dyplus, *dyminus, *dvy;
  static  double *dzplus, *dzminus, *dvz;

  static  double *right, *diag, *sup, *sub;

  void    tridiag(int n);

  int     *source_x0,    *source_nx;
  int     *source_y0,    *source_ny;
  int     *source_z0,    *source_nz;
  double  **source_xsum, **source_ysum, **source_zsum;

public:

  static  double* bc_deriv, * bc_lin,  * bc_coef, * bc_const;
  static  double* bc_pump,  * bc_pow,  * bc_Kn;
  static  double* bc_pump2, * bc_pow2, * bc_Kn2;

  static  char   **bc_id;
  static  int    bc_type_num;

  static  int   xsize, ysize, zsize;
  static  long  Size, xysize;

  static  double *xgrid, *ygrid, *zgrid;
  static  double *xcoord, *ycoord, *zcoord;

  static class RegionObj     *Region;
  static class GridObj       *Grid;
  static class BCarrayObj    *BCarray;

  static RegionObj *get_region() { return Region; }
  static void setStaticData(RegionObj &r, GridObj &GO, BCarrayObj &BCA);

  static void    init_tridiag(int);
  static void    kill_tridiag();

  struct TermStruct *Currents;
  class ExpressionObj *Current;

  static  char SAME_CURRENT;

  void getCurrents(TokenString &, int simID, VarList *VL, double *);
  void evaluateCurrents();
  void killCurrents();
  void printCurrents(FILE *out = (FILE *)stderr);

  double     D;
  double    *Diff;
  VectorObj *kuptake;

  bool    tortDefined; // true if tortuosity is defined 
  int     source_num;

  long           *ptype;
  int            *bccond;
  int            fieldObstNum;
  VolumeObjClass *fieldObstArray;

  double  ICa;

  double  bgr; // initial average;
  char    *ID;
  double  Time;

  void allocate(const char *id = "Field");

  FieldObj(int=0, double=1.0, double tot=0.0, const char *id = "Ca");
  FieldObj(TokenString &, const char *id = "Ca");
  FieldObj(const FieldObj &);

  ~FieldObj();

  FieldObj& operator=(const FieldObj &f);
  FieldObj& operator=(const VectorObj &f);
  FieldObj& operator=(double (*)(double, double, double));
  FieldObj& operator=(const double);

  int point_type(int x, int y, int z) { return Region->point_type(x, y, z, fieldObstNum, fieldObstArray); }
  // int point_type(long ind)            { return Region->point_type(ind,     fieldObstNum, fieldObstArray); }

  void   init_boundary();
  signed long location_to_index(double, double, double, bool);

  void set_source(int, double, double, double, double, double, double);
  void set_source_xyz(int id, double c0, int ix, double sigma, 
                      double *coord, double *grid, double *dvArray, int N,
                      int *source_xyz0, int *source_nxyz, double **source_xyzsum);

  void adjust_sources(FieldObj f);
  void add_sources(VectorObj &, double);
  void add_sources(double dt) { add_sources(*this, dt); }

  void set_bc(int boxid, int, int, int, int, int, int);
  void set_bc(int boxid, int);
  void set_bc(int boxid, const char *, const char *, const char *,
                         const char *, const char *, const char *);
  void set_bc(int boxid, const char *);
  int  get_bc(int boxid, int surface) { return bccond[boxid * DIMENSIONALITY * 2 + surface]; }  

  double get_value(double, double, double);

  long *get_ptype() { return ptype; }

  double getD() { return D; }
  void   setD(double d) { D = d; }

  // double getGhost(long i);

  void nablaX(long i, int ix, long bcvar, double &, double &, double &, double &);
  void nablaY(long i, int iy, long bcvar, double &, double &, double &, double &, double grid = 1.0);
  void nablaZ(long i, int iz, long bcvar, double &, double &, double &, double &, double grid = 1.0);

  void Run1D(double, double, double *, double *);

  void Run2Dx(double, double, double, double *, double *);
  void Run2Dy(double, double, double, double *, double *);

  void Run3Dx(double, double, double, double, double *, double *);
  void Run3Dy(double, double, double, double, double *, double *);
  void Run3Dz(double, double, double, double, double *, double *);

  void cleanBoundaries(); 
  void print_ptype(long i);
  void print_ptype();

  double interpolate(double x, double y, double z);

  //  void init(double T=0.0) { Time = T; *this = total; };

  double average1D (double *V = 0);
  double average2D (double *V = 0);
  double average3D (double *V = 0);

  double average(double *V = 0)  { 
	 switch(DIMENSIONALITY) { case 1:  return average1D (V);
						      case 2:  return average2D (V);
						      case 3:  return average3D (V);
							  default: return 0.0;
	 }
  }

  double gain()  { 
    double V;
    double c = average(&V) - bgr;
    return c * V;
  }    

  void   print3d(int n=3);
  friend void print3d(long *);
};


//*****************************************************************************
//                         B U F F E R   C L A S S
//*****************************************************************************

class BufferObj : public FieldObj
{
 protected:
  
 public:

  VectorObj *total;
  ExpressionObj *kplus;
  ExpressionObj *kminus;
  char coopType;

  ~BufferObj() { delete total; };

  BufferObj(TokenString &, const char *, const char coopFlag = 0);

  BufferObj(const BufferObj &b) : FieldObj(b)
    { kminus = b.kminus; kplus = b.kplus; coopType = b.coopType; 
      total = new VectorObj(*(b.total)); }

  BufferObj& operator=(const BufferObj &b)
    { 
    *(FieldObj *)this = b; 
    *total = *(b.total);
    kplus = b.kplus; kminus = b.kminus; coopType = b.coopType;
    return *this; 
    }

  BufferObj& operator=(const FieldObj &f)
    { *(FieldObj *)this = f; return *this; }

  BufferObj& operator=(const VectorObj &v)
    { *(VectorObj *)this = v; return *this; }

  BufferObj& operator=(double (*f)(double, double, double))
    { *(FieldObj *)this = f; return *this; }

  BufferObj& operator=(const double val)
    { *(VectorObj *)this = val; return *this; }

  void setMembraneBound(double Total, int C_not_N,  int boxid);
  void setMembraneLayer(double Total, double depth, int boxid);

};

//****************************************************************************
//*                   C L A S S   F I E L D   A R R A Y
//****************************************************************************

class BufferArray
{
private:

public:

  char *ID;
  int buf_num;
  int nonCoopNum;
  BufferObj **array;

  BufferArray() { array = 0; ID = 0; buf_num = 0; nonCoopNum = 0; }
  BufferArray(TokenString &params);
  BufferArray(const BufferArray &bufs)  { allocateCopy(bufs); }

  ~BufferArray()  { destroy(); }
 
  void operator=(const BufferArray &bufs) {  // beware: group ID ignored
    if ( (buf_num = bufs.buf_num) > 0 )
      for (int i = 0; i < buf_num; i++) *array[i] = *bufs.array[i];
  }

  void getCurrents(TokenString &TS, int simID, VarList *VL, double *tptr)  {  
    if (buf_num) for (int i = 0; i < buf_num; i++) array[i]->getCurrents(TS, simID, VL, tptr);
  };

  void addTime(double dt)  {  
    if (buf_num) for (int i = 0; i < buf_num; i++) array[i]->Time += dt;
  };

  void setTime(double T)  {  
    if (buf_num) for (int i = 0; i < buf_num; i++) array[i]->Time = T;
  };

  void killCurrents()  {  
    if (buf_num) for (int i = 0; i < buf_num; i++) array[i]->killCurrents();
  };


  void allocateCopy(const BufferArray&);

  void destroy()  {  if (buf_num) {
                       delete [] ID;
                       for (int i = 0; i < buf_num; i++) delete array[i];
                       delete [] array;
		     }
  }

  double gain() {
    double bound = 0.0;
    for (int i = 0; i < buf_num; i++) 
		bound += (array[i]->coopType - 1) * array[i]->gain();
  return bound;
  }
};

//****************************************************************************


#endif
