/************************************************************************
*
*                   Calcium Calculator (CalC)
*            Copyright (C) 2001-2019 Victor Matveev
*
*                            box.cpp
*
*               Geometry definition classes:
* 
*               Generalized box class BoxObj
*               Box as an obstacle: class ObstBoxObj
*               Set of boxes and obstacles: RegionObj
*
**************************************************************************
 
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
#include <math.h>
#include <stdarg.h>
#include "syntax.h"
#include "grid.h"
#include "box.h"

class GridObj *VolumeObjClass::Grid = 0;
class GridObj *RegionObj::Grid = 0;

extern void setFieldByFunction(TokenString &TS, long p, bool    *Diff, const char *errStr);


//***************************************************************
//                I M P L E M E N T A T I O N S
//***************************************************************

int VolumeObjClass::inside(int ix, int iy, int iz) {
	int b;
	double dist;

							 if (ix < 0 || ix >= Grid->xsize) return isObstacle << 6; // = 64 = _INSIDE_ if isObstacle==1
	if ( DIMENSIONALITY > 1) if (iy < 0 || iy >= Grid->ysize) return isObstacle << 6; // = 64 = _INSIDE_ if isObstacle==1
	if ( DIMENSIONALITY > 2) if (iz < 0 || iz >= Grid->zsize) return isObstacle << 6; // = 64 = _INSIDE_ if isObstacle==1

	switch (volumeType) 
	{
	case 'b':  // Box    
		                         b =     ((Grid->xcoord[ix] >= p1) & (Grid->xcoord[ix] <= p2));
		if (DIMENSIONALITY > 1)  b = b & ((Grid->ycoord[iy] >= p3) & (Grid->ycoord[iy] <= p4));
	    if (DIMENSIONALITY > 2)  b = b & ((Grid->zcoord[iz] >= p5) & (Grid->zcoord[iz] <= p6));
	    return (isObstacle ^ b) << 6;

	case 's': // Sphere
		dist = sqrt( (Grid->xcoord[ix] - p1) * (Grid->xcoord[ix] - p1) +
				 	 (Grid->ycoord[iy] - p2) * (Grid->ycoord[iy] - p2) +
					 (Grid->zcoord[iz] - p3) * (Grid->zcoord[iz] - p3) );
		return (isObstacle ^ (dist <= p4)) << 6;

	case 'd': // Disk
		dist = sqrt( (Grid->xcoord[ix] - p1) * (Grid->xcoord[ix] - p1) +
					 (Grid->ycoord[iy] - p2) * (Grid->ycoord[iy] - p2) );
		return (isObstacle ^ (dist <= p3)) << 6;

	case 'c': // Cylinder
		dist = sqrt( (Grid->xcoord[ix] - p1) * (Grid->xcoord[ix] - p1) +
		  			 (Grid->ycoord[iy] - p2) * (Grid->ycoord[iy] - p2) );
		return (isObstacle ^ ((dist <= p5) & (Grid->zcoord[iz] >= p3) & (Grid->zcoord[iz] <= p4))) << 6;

	case 'f': // Formula
		int intermediate = isInside[ix + iy*Grid->xsize + iz*Grid->xysize];
		return (isObstacle ^ intermediate) << 6;
	}

	return 0;
}


//***************************************************************

int VolumeObjClass::point_type(int ix, int iy, int iz)
{
	int  surface;

	if ( (surface = inside(ix, iy, iz)) > 0 )
	{
		if (! inside(ix - 1, iy, iz) ) surface |= SURF_XMIN;
		if (! inside(ix + 1, iy, iz) ) surface |= SURF_XMAX;
		if (DIMENSIONALITY == 1) return surface;
		if (! inside(ix, iy - 1, iz) ) surface |= SURF_YMIN;
		if (! inside(ix, iy + 1, iz) ) surface |= SURF_YMAX;
		if (DIMENSIONALITY == 2) return surface;
		if (! inside(ix, iy, iz - 1) ) surface |= SURF_ZMIN;
		if (! inside(ix, iy, iz + 1) ) surface |= SURF_ZMAX;
	}

	return surface;
}
//***************************************************************

VolumeObjClass::VolumeObjClass(TokenString &pars, long pos, int isObst)
{
	int    nPars = pars.tokens_to_eol(pos + 1);
	if (VERBOSE) isObst ? fprintf(stderr, "    Adding obstacle:  ") : fprintf(stderr, "    Adding enclosure: ");
	isInside = 0;

	switch (DIMENSIONALITY)
		{
		case 1:   // 1D
			if (nPars < 2) pars.errorMessage(pos+1, 0, "At least two arguments needed in 1D geometry volume definition");
			pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'E');
			this->setSize(p1, p2, 0.0, 0.0, 0.0, 0.0);
			this->setPars(isObst, 'b', p1, p2);
			if (VERBOSE) fprintf(stderr,"interval [%g,%g] \n", p1, p2);
			if (nPars == 2) break; 
			this->volumeType = 'f';
			this->tokenPosition = pos+3;
			break;

		case 2:  // 2D
			switch (nPars)
			{
			case 4: // it's a box
				pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'd', &p4, 'E');
				this->setPars(isObst, 'b', p1, p2, p3, p4);
				this->setSize(p1, p2, p3, p4, 0.0, 0.0);
				if (VERBOSE) fprintf(stderr,"box [%g, %g] x [%g, %g]", p1, p2, p3, p4);
				return;

			case 3: // it's a disk
				pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'E');
				if ( p3 < 0 ) pars.errorMessage(pos+1, makeMessage("Check disk dimensions: radius=%g < 0", p2));
				this->setPars(isObst, 'd', p1, p2, p3);
				this->setSize(p1-p3, p1+p3, p2-p3, p2+p3, 0.0, 0.0);
				if (VERBOSE) fprintf(stderr,"disk: center=[%g, %g], radius=%g \n", p1, p2, p3);
				break;

			default: // Process an analytic expression
				if (nPars < 4) pars.errorMessage(pos+1, 0, "At least 3 arguments needed in 2D volume definition ");
				if (!isObst) { 
					pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'd', &p4, 'E'); 
					this->setSize(p1, p2, p3, p4, 0.0, 0.0);
					pos += 4; 
					if (VERBOSE) fprintf(stderr,"Analytic shape within box [%g, %g] x [%g, %g]", p1, p2, p3, p4);
				}
				if ( !pars.equal(pos+1, "=") ) pars.errorMessage(pos+1, 0, "Expected equal sign followed by formula");
				this->setPars(isObst, 'f', p1, p2, p3, p4);
				this->tokenPosition = pos+1;
			}
			break;

		case 3:  // 3D
			switch (nPars)
			{
			case 6: // it's a box
				pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'd', &p4, 'd', &p5, 'd', &p6, 'E');
				this->setPars(isObst, 'b', p1, p2, p3, p4, p5, p6);
				this->setSize(p1, p2, p3, p4, p5, p6);
				if (VERBOSE) fprintf(stderr, "box [%g, %g] x [%g, %g] x [%g, %g]\n", p1, p2, p3, p4, p5, p6);
				break;

			case 4: // it's a sphere
				volumeType = 's';
				pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'd', &p4, 'E');
				if ( p4 < 0 ) pars.errorMessage(pos+1, makeMessage("Check sphere dimensions: radius=%g < 0", p2));
				this->setPars(isObst, 's', p1, p2, p3, p4);
				this->setSize(p1-p4, p1+p4, p2-p4, p2+p4, p3-p4, p3+p4);
				if (VERBOSE) fprintf(stderr, "sphere: center=[%g, %g, %g], radius=%g \n", p1, p2, p3, p4);
				break;

			case 5: // it's a cylinder
				pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'd', &p4, 'd', &p5, 'E');
				if ( p3 > p4 || p5 < 0 ) pars.errorMessage(pos+1, makeMessage("Check cylinder dimensions: zmax=%g < zmin=%g or radius=%g < 0", p3, p4, p5));
				this->setPars(isObst, 'c', p1, p2, p3, p4, p5);
				this->setSize(p1-p5, p1+p5, p2-p5, p2+p5, p3, p4);
				if (VERBOSE) fprintf(stderr,"cylinder: base center=(%g, %g, %g), z in [%g, %g], radius=%g \n", p1, p2, p3, p3, p4, p5);
				break;

			default: // Process an analytic expression
				if (nPars < 6) pars.errorMessage(pos+1, 0, "At least 4 arguments needed in 3D volume definition ");
				if (!isObst) { 
					pars.trail_pars(pos, 'd', &p1, 'd', &p2, 'd', &p3, 'd', &p4, 'd', &p5, 'd', &p6, 'E');
					this->setSize(p1, p2, p3, p4, p5, p6);
					pos += 6; 
					if (VERBOSE) fprintf(stderr,"Analytic shape within box [%g, %g] x [%g, %g] x [%g, %g]\n", p1, p2, p3, p4, p5, p6);
				}
				if ( !pars.equal(pos+1, "=") ) pars.errorMessage(pos+1, 0, "Expected equal sign followed by formula");
				this->setPars(isObst, 'f', p1, p2, p3, p4, p5, p6);
				this->tokenPosition = pos+1;
			}
	} 
}

//***************************************************************

void   VolumeObjClass::computeFormulas(TokenString &pars) 
{
	  if (volumeType != 'f') return;

	  isInside = new bool[Grid->Size];
	  if (VERBOSE) isObstacle ? fprintf(stderr, "  # Obstacle formula: ") : fprintf(stderr, "  # Volume formula: ");

	  switch (DIMENSIONALITY)
		{
		case 1:  
			setFieldByFunction(pars, tokenPosition, isInside, "Wrong conditional expression for 1D volume");
			break;

		case 2:  
			setFieldByFunction(pars, tokenPosition, isInside, "Wrong conditional expression for 2D volume");
			break;

		case 3:  
			setFieldByFunction(pars, tokenPosition, isInside, "Wrong conditional expression for 3D volume");
			
	} 
}
//***************************************************************
//*   C L A S S   R E G I O N  I M P L E M E N T A T I O N S
//***************************************************************

RegionObj::RegionObj(TokenString &pars)
{

	encl_num = pars.token_count("volume");
	obst_num = pars.token_count("obstacle");
	total    = encl_num + obst_num;

	if (VERBOSE)       fprintf(stderr,"\n#### Building geometry with %d enclosures and %d obstacles: \n", encl_num, obst_num); 
	if (!encl_num) globalError(StrCpy("No diffusion volume definitions in script file"));

	// initialize( encl_num, obst_num );

	enclosures =            new VolumeObjClass[encl_num];
	obstacles  = obst_num ? new VolumeObjClass[obst_num] : 0;

	              for (int i = 0; i < encl_num; i++)  enclosures[i] = VolumeObjClass(pars, pars.token_index("volume",   i + 1), 0);
	if (obst_num) for (int i = 0; i < obst_num; i++)   obstacles[i] = VolumeObjClass(pars, pars.token_index("obstacle", i + 1), 1);

	XMin = enclosures[0].xmin;   XMax = enclosures[0].xmax;
	YMin = enclosures[0].ymin;   YMax = enclosures[0].ymax;
	ZMin = enclosures[0].zmin;   ZMax = enclosures[0].zmax;

	for (int i = 0; i < encl_num; i++) 
	{
		if (enclosures[i].xmin >= enclosures[i].xmax) 
			pars.errorMessage(pars.token_index("volume", i + 1) + 1, makeMessage("Check dimensions of volume %d: xmin >= xmax\n", i+1));
		if (enclosures[i].xmin < XMin) XMin = enclosures[i].xmin;  
		if (enclosures[i].xmax > XMax) XMax = enclosures[i].xmax; 
		if (DIMENSIONALITY == 1) continue;
		if (enclosures[i].ymin >= enclosures[i].ymax) 
			pars.errorMessage(pars.token_index("volume", i + 1) + 3, makeMessage("Check dimensions of volume %d: ymin >= ymax\n", i+1));
		if (enclosures[i].ymin < YMin) YMin = enclosures[i].ymin;  
		if (enclosures[i].ymax > YMax) YMax = enclosures[i].ymax; 
		if (DIMENSIONALITY == 2) continue;
		if (enclosures[i].zmin >= enclosures[i].zmax) 
			pars.errorMessage(pars.token_index("volume", i + 1) + 5, makeMessage("Check dimensions of volume %d: zmin >= zmax\n", i+1));
		if (enclosures[i].zmin < ZMin) ZMin = enclosures[i].zmin;  
		if (enclosures[i].zmax > ZMax) ZMax = enclosures[i].zmax; 
	}
}

//***************************************************************

/*
void RegionObj::initialize(int vNum, int oNum) {  

	encl_num   = vNum; 
	obst_num   = oNum;

	total = encl_num + obst_num ;

	enclosures = new VolumeObjClass[encl_num];
	for (int i = 0; i < encl_num) enclosures[i].isObstacle = 0;

	obstacles  = new VolumeObjClass[obst_num];
	for (int i = 0; i < obst_num)  obstacles[i].isObstacle = 1;
}
*/

//***************************************************************

RegionObj::~RegionObj()
{
	for (int i = 0; i < encl_num; i++) 
		if (enclosures[i].volumeType == 'f' && enclosures[i].isInside)  delete [] enclosures[i].isInside;

	for (int i = 0; i < obst_num; i++) 
		if (obstacles[i].volumeType == 'f' && obstacles[i].isInside)  delete [] obstacles[i].isInside;

	delete [] enclosures;
	delete [] obstacles;
}

//***************************************************************
// point_type bit fields are:
// bit # is:    ...8     |   7    |  6   |  5   |  4   |  3   |  2   |  1  
// meaning:   box number | inside | zmax | zmin | ymax | ymin | xmax | xmin
//***************************************************************


int RegionObj::point_type(int ix, int iy, int iz, int fieldObstNum, class VolumeObjClass *OBO)
{
	int    ptype, boxid = 0;
	int    encl_bnd = SURF_MASK;
	int    obst_bnd = 0;  
	int    i, j, inside = 0;

	for (i = 0; i < encl_num; i++)
	{  
		ptype = enclosures[i].point_type(ix, iy, iz); 
		if (ptype & _INSIDE_)    
		{                                      //  Point is inside if it is
			inside = _INSIDE_;                     //  inside any of the enclosures.
			encl_bnd &= ptype;                     //  Inside point is on a boundary only
			if (ptype  & SURF_MASK )  boxid = i;   //  if it belongs to the same boundary
		}                                      //  of all the enclosures containing it
	}    

	if (encl_bnd == SURF_MASK)  encl_bnd = 0;  

	for (j = 0; j < obst_num; j++)
	{  
		ptype = obstacles[j].point_type(ix, iy, iz);
		if (ptype == 0) inside = 0;              // Point is outside if it belongs
		else                                   // to any of the obstacles.
			if (ptype  & SURF_MASK )  
			{ 
				obst_bnd |= (ptype & SURF_MASK);    
				boxid = j + encl_num; 
			}                                 
	}                                       


	if (fieldObstNum)
		for (j = 0; j < fieldObstNum; j++)
		{  
			ptype = OBO[j].point_type(ix, iy, iz);
			if (ptype == 0) inside = 0;              // Point is outside if it belongs
			else                                   // to any of the obstacles.
				if (ptype  & SURF_MASK )  
				{                                  
					obst_bnd |= (ptype & SURF_MASK);                 
					boxid = j + encl_num +  obst_num; 
				}                                 
		}                                        

	if ( inside == 0 ) return 0;

	int bnd = encl_bnd | obst_bnd;

	if (bnd == 0) return inside;
	/*
	if ( ((bnd & SURF_XMIN) && (bnd & SURF_XMAX)) ||
	((bnd & SURF_YMIN) && (bnd & SURF_YMAX)) ||
	((bnd & SURF_ZMIN) && (bnd & SURF_ZMAX)) )
	throw makeMessage("Ambiguity: two surfaces touch. Please enlarge the obstacle dimensions, extending it beyond the domain boundary, to avoid surface ambiguity");
	else */  
	return  _INSIDE_ + (boxid << 7) + bnd;

}
//***************************************************************

void   RegionObj::computeFormulas(TokenString &pars) 
{
	for (int i = 0; i < encl_num; i++) enclosures[i].computeFormulas(pars);
	for (int i = 0; i < obst_num; i++)  obstacles[i].computeFormulas(pars);

}
//***************************************************************
/*
int RegionObj::point_type(long ind, int fieldObstNum, class ObstBoxObj *OBO)
{
	int iz = ind / Grid->xysize;
	ind -=  iz * Grid->xysize;
	int iy = ind / Grid->xsize;
	int ix = ind - iy * Grid->xsize;
	return point_type(ix, iy, iz, fieldObstNum, OBO);
}
*/
//***************************************************************

void print_type(int x)
{
	if (x == 0) fprintf(stderr,"[OUT] ");
	else
	{
		if (x & _INSIDE_) fprintf(stderr,"[IN");
		if (!(x & SURF_MASK)) fprintf(stderr,"] ");
		else
		{
			if (x & SURF_XMIN) fprintf(stderr," XMIN");
			if (x & SURF_XMAX) fprintf(stderr," XMAX");
			if (x & SURF_YMIN) fprintf(stderr," YMIN");
			if (x & SURF_YMAX) fprintf(stderr," YMAX");
			if (x & SURF_ZMIN) fprintf(stderr," ZMIN");
			if (x & SURF_ZMAX) fprintf(stderr," ZMAX");
			fprintf(stderr," %d] ", x>>9);
		}
	}
}
#ifdef CALC_DEBUG
#endif
//***************************************************************
//                       T H E   E N D
//***************************************************************

