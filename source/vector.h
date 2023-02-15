/*****************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *                 Copyright (C) 2001-2019 Victor Matveev
 *
 *                                vector.h
 *
 *                   A basic linear algebra vector class 
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

#ifndef CALC_VectorObj_H_included
#define CALC_VectorObj_H_included

extern int VERBOSE;

class VectorObj {

 protected:

 public:

  double *elem;
  long   size;

  VectorObj(long n = 0)  { elem = new double[size = n]; };

  ~VectorObj() { delete [] elem; }

  VectorObj(const VectorObj &v);

  void set_size(long n) { size = n; }
  long get_size()  { return size; };

  double  operator[](long index) const  { return elem[index]; };
  double& operator[](long index) { return elem[index]; };

  double *operator()() const  { return elem; };

  VectorObj& operator=(const double *e);
  VectorObj& operator=(const double val);
  VectorObj& operator=(const VectorObj &v);

  void import(const char *);
  void Export(const char *);

  friend VectorObj *vabs(const VectorObj &v1);

  friend VectorObj operator-(const VectorObj &v);
  friend VectorObj operator-(const VectorObj &v1, const VectorObj &v2);
  friend VectorObj operator-(const VectorObj &v1, const double val);
  friend VectorObj operator-(const double val, const VectorObj &v1);

  friend void operator-=(VectorObj &v1, const VectorObj &v2);
  friend void operator-=(VectorObj &v1, const double val);

  friend VectorObj operator+(const VectorObj &v1, const VectorObj &v2);
  friend VectorObj operator+(const VectorObj &v1, const double val);
  friend VectorObj operator+(const double val, const VectorObj &v1);

  friend void operator+=(VectorObj &v1, const VectorObj &v2);
  friend void operator+=(VectorObj &v1, const double val);

  friend double operator*(const VectorObj &v1, const VectorObj &v2);
  friend VectorObj operator%(const VectorObj &v1, const VectorObj &v2);

  friend VectorObj operator*(const VectorObj &v1, const double val);
  friend VectorObj operator*(const double val, const VectorObj &v1);
  friend VectorObj operator%=(VectorObj &v1, const VectorObj &v2);
  friend VectorObj operator*=(VectorObj &v1, const double val);

  double norm() const;
  double normMod(int n) const;
  //double deriv_norm() const;
  double L_infty_norm() const;
  double norm_L1() const;

  double checkerBoardNorm1D();
  double checkerBoardNorm2D(int S1);
  double checkerBoardNorm3D(int S1, int S2);

  double checkerBoardNorm(int nDIMS, int S1 = 0, int S2 = 0) { 
     switch(nDIMS) { case 1: return checkerBoardNorm1D();
                     case 2: return checkerBoardNorm2D(S1);
					 case 3: return checkerBoardNorm3D(S1, S2);
					 default: return 0.0;
	 }
  }

  void   print(int n=4) const;
};

#endif
