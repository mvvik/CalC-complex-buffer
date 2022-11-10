******************************************************************************
 
                        Calcium Calculator (CalC)
                  Copyright (C) 2001-2022 Victor Matveev
                         LBM/NIDDK/NIH and DMS/NJIT

                             matveev@njit.edu
                       http://web.njit.edu/~matveev

                             October 10, 2022

 ************************************************************************
 
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

 ************************************************************************

==============================================================================
                             1.  COMPILATION
==============================================================================

NOTE: platform-specific definitions are contained in the short source header
      file called "PlatformSpecific.h".

1) Windows: 
   
   Use MSDN Visual Studio to compile the code. The Visual Studio "Project" file
   is contained within the "source" folder. 

   Alternatively, you can install the Linux emulator "cygwin" (http://www.cygwin.com")
   on your Windows machine, and follow the UNIX installation instructions below. 

2) Mac OS X users: you will need to install Xcode and XQuartz. There should be some 
   version of C++ compiler such as g++, and some kind of make utility to complile the 
   code. Then, run "make". That is all.
   
3) Non-Mac UNIX / LINUX / CygWin

   Replace Makefile with Makefile.no_glut (remove the extension). Then, run "make"
   as usual. This will compile a version of the program without OpenGl/GLUT/FreeGlut 
   graphics. Alternatively, install freeglut and change the linker option in the 
   Makefile to properly link your installation of GLUT/FreeGlut library. See also 
   "PlatformSpecific.h" header file in the source directory.

==============================================================================
                      2.  EXECUTION (ALL PLATFORMS)
==============================================================================

Windows OS only: make sure that the included "freeglut.dll" is present in the folder 
with your calc executable and script file.

If a model script does not use command-line parameters, simply click on the 
executable and enter the script file name when prompted (the script file should
reside in the same directory).  

If a script utilizes command-line parameters, then open a shell (for Windows, 
go to Start Menu -> Run -> "cmd"), go to the directory where the executable 
"cwinxx.exe" is residing, and type 

   cwinxxx.exe filename [ parlist ]
 
where "fileName" is the name of the script file describing the simulation, and
optional parList is a space-separated list of command-line parameters (see
http://web.njit.edu/~matveev/calc/manual.html#pars).

In order to monitor program output and error messages, include the statement 
"verbose = 4" (or higher) in your script: this will prevent CalC from 
auto-terminating upon completing the simulation.

Note that script debugging is easier when using UNIX, since output stream 
flushing is problematic in Windows, and diagnostic error messages are sometimes 
not displayed.

==============================================================================
                             3.  DATA OUTPUT:
==============================================================================

o All platforms: 

Simulation results can be saved to a file in real time using "mute" plot 
statements, ASCII files are produced that are readable by any graphics package 
such as MATLAB (Mathworks, Inc). See demo scripts and refer to the manual
at http://web.njit.edu/~matveev/calc/manual.html#mute for details

The "binary" plot type allows to save an entire concentration field at several 
time points during the simulation, and can be read and displayed using MATLAB
via MATLAB scripts provided in the "examples" directory and on the on-line
demo script page.

o Windows and Mac OS X:

Specify "plot.method gl" within your script to make time-dependent variable 
plots, and 1D and 2D concentration plots (consult the manual). On Windows, make
sure that "freeglut.dll" is in your home directory. On Mac OS X, GLUT is 
preinstalled (but depricated). On other platforms, you have to have GLUT/freeglut 
installed on your computer, and change the linker directive in the Mekefile.

o UNIX platforms (including cygwin and Mac OS X):

You can use the "xmgrace" graphics application, which allows viewing the 
simulation results in real time. This involves piping the program output into 
xmgrace, by executing one of the following commans:

  calc scriptFileName | xmgrace -pipe

where "calc" is the name of your CalC executable, "scriptfileName" is the name 
of the simulation script. In this case the script file should contain an 
instruction "plot.method xmgr".

A "homebrew" installation of xmgrace is easy: 
   https://formulae.brew.sh/formula/grace

o MATLAB integration:

Like any system command, CalC can be launched from MATLAB simply by executing
the command "system('.\cwinxx.exe script.par')". You can then collect the data
by reading the "mute plot" files. See CalC demo script web page for further details.

==============================================================================
                             3. DOCUMENTATION:
==============================================================================

o  A hyper-text CalC script syntax manual can be found at URL
   http://www.calciumcalculator.org

o  Also, examine the included commented example script files (file extension 
   ".par"): they reside in the "examples" directory.

o  Example script files are also available at
   http://web.njit.edu/~matveev/calc/scripts.html

******************************************************************************
==============================================================================
