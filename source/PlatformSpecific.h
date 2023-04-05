/************************************************************************
 *
 *                   Calcium Calculator (CalC)
 *            Copyright (C) 2001-2022 Victor Matveev
 *                    LBM/NIDDK/NIH and DMS/NJIT
 *
 *                      PlatformSpecific.h
 *
 *  This header contains *only* platform-specific pre-processor directives
 * 
 *  OpenGl/GLUT/FreeGlut is only specified for APPLE and WINDOWS platforms
 *  If you know how to link FreeGLut to other platforms, make changes below
 *
 ************************************************************************/

//  If GLUT library is NOT installed, uncomment the following line:
// (or use an alternative makefile target: "make noGraphs")

// #define _NO_GLUT_

#ifdef _WIN32
	#ifndef _NO_GLUT_
		#include <GL/freeglut.h>
	#endif
#else

    #define _isnan(X)  isnan(X)
    #define _finite(X) isfinite(X)

    #ifdef __APPLE__
      #include <OpenGL/glu.h>
      #include <GLUT/GLUT.h>
    #endif

#endif


 /************************************************************************
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
