#
#   Calcium Calculator * CalC version 7.9.5 * makefile
#   Victor Matveev (C) Copyright 2001-2019 (GPLv3) 
#   
#   Department of Mathematical Sciences
#   New Jersey Institute of Technology
#
#   E-mail: matveev@njit.edu
#   Http://web.njit.edu/~matveev
#
#   Developed at LBM, NIDDK, NIH and DMS/NJIT
#   Supported in part by NSF DMS0417416, DMS0817703, DMS517085
#
#   Please use "gmake" to compile if "make" fails
#
# ************************************************************************
# 
#    This file is part of Calcium Calculator (CalC).
#
#    CalC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CalC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CalC.  If not, see <https://www.gnu.org/licenses/>
#
# ************************************************************************

D = source
VPATH = . : $D 

flags = -O ${CXXFLAGS}
llibs = -lm ${LDFLAGS}
CXX = g++

objects  = syntax.o vector.o table.o peak.o box.o grid.o field.o \
           interpol.o fplot.o gate.o loop.o simulation.o markov.o

includes = box.h vector.h field.h fplot.h syntax.h \
           gate.h table.h peak.h interpol.h loop.h grid.h simulation.h markov.h

sources  = calc.cpp box.cpp vector.cpp field.cpp fplot.cpp syntax.cpp \
           calc.cpp table.cpp peak.cpp interpol.cpp loop.cpp grid.cpp simulation.cpp markov.cpp

calc:   moveobjects ${objects} ${includes} ${sources} 
	   $(CXX) $D/calc.cpp ${objects} ${llibs} ${flags} -o calc
	   -mv *o source/
	   echo '*** Compilation completed successfully ***'

moveobjects: 
	   -mv source/*o ./

syntax.o : syntax.cpp syntax.h 
	   $(CXX) $D/syntax.cpp ${flags} -c

vector.o : vector.cpp vector.h
	   $(CXX) $D/vector.cpp ${flags} -c          

table.o  : table.cpp table.h syntax.h
	   $(CXX) $D/table.cpp  ${flags} -c

peak.o  :  peak.cpp peak.h syntax.h
	   $(CXX) $D/peak.cpp  ${flags} -c

markov.o : markov.cpp markov.h syntax.h
	   $(CXX) $D/markov.cpp  ${flags} -c

box.o :    box.cpp box.h grid.h syntax.h 
	   $(CXX) $D/box.cpp ${flags} -c

grid.o :   grid.cpp grid.h box.h syntax.h 
	   $(CXX) $D/grid.cpp ${flags} -c

field.o :  field.cpp field.h box.h grid.h vector.h syntax.h  
	   $(CXX) $D/field.cpp ${flags} -c  

interpol.o : interpol.cpp syntax.h field.h box.h grid.h vector.h
	   $(CXX) $D/interpol.cpp  ${flags} -c

gate.o  :  gate.cpp gate.h syntax.h field.h box.h grid.h vector.h table.h peak.h interpol.h simulation.h markov.h
	   $(CXX) $D/gate.cpp  ${flags} -c

fplot.o :  fplot.cpp fplot.h box.h field.h vector.h syntax.h gate.h table.h peak.h interpol.h simulation.h markov.h
	   $(CXX) $D/fplot.cpp ${flags} -c

simulation.o : simulation.cpp simulation.h box.h field.h vector.h syntax.h gate.h table.h peak.h interpol.h fplot.h markov.h
	   $(CXX) $D/simulation.cpp ${flags} -c

loop.o  :  loop.cpp loop.h fplot.h box.h field.h vector.h syntax.h gate.h table.h peak.h interpol.h markov.h
	   $(CXX) $D/loop.cpp  ${flags} -c



