******************************************************************************
 
                        Calcium Calculator (CalC)
                  Copyright (C) 2001-2019 Victor Matveev
                         LBM/NIDDK/NIH and DMS/NJIT

                             matveev@njit.edu
                       http://web.njit.edu/~matveev

                               May 20, 2019

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

1) Windows: 

   If the stand-alone executables do not run, download and unzip the file 
   "CalC_win32_xxx.zip" (where "xxx" denotes the current version number, e.g.
   "CalC_win32_6.9.5.zip"). It contains the C++ source code (subdirectory "source"), 
   and the example script files (subdirectory "examples"). You have to use MSDN
   Visual Studio to compile the code. Compilation will produce an executable 
   named "Calcium Calculator.exe". 

   Alternatively, you can install the Linux emulator "cygwin" (http://www.cygwin.com")
   on your Windows machine, and follow the UNIX installation instructions below. This
   would enable you to use the "xmgrace" graphics program to view the simulation
   results in real time, provided you install the xmgrace, make, g++ and X11 packages
   as part of the cygwin installation.

2) UNIX platforms (including cygwin and OS X)

  If the stand-alone executables do not run (see below), place the downloaded source 
  archive file into a new separate directory, then unpack and compile it using:

    gunzip calc_unix_xxx.tgz
    tar xf calc_unix_xxx.tar
    make
  
  Here "xxx" denotes the current version number (e.g. "calc_unix_6.9.5.tgz")
  If "make" does not do anything, try "gmake". Compilation will produce an 
  executable named "calc" or "calc.exe"

3) Mac OS X users: you will need to install Xcode, as well as the g++ (gcc) compiler 
   and the make utility in order to be able to compile the Unix version of the 
   code.

==============================================================================
                             2.  EXECUTION:
==============================================================================

1) Windows: 

If a model script does not use command-line parameters, simply click on the 
executable and enter the script file name when prompted (the script file should
reside in the same directory). In order to monitor program output and error 
messages, include the statement "verbose = 4" (or higher) in your script: this 
will prevent CalC from auto-terminating upon completing the simulation.

If a script utilizes command-line parameters, then open a DOS window 
(Start Menu -> Run -> "cmd"), go to the directory where the executable 
"cwinxx.exe" is residing, and type 

   cwinxx.exe filename [ parlist ]
 
where "fileName" is the name of the script file describing the simulation, and
optional parList is a space-separated list of command-line parameters (see
http://web.njit.edu/~matveev/calc/manual.html#pars).

Note that script debugging is easier with a unix- or cygwin-compiled code, 
since output stream flushing is problematic in Windows, and diagnostic error 
messages are sometimes not displayed.

2) Mac OS X (and other UNIX platforms)

First, make sure the executable has permissions to run by typing

  chmod +x cmacxxx

where "xxx" is the CalC version number (e.g. "chmod +x cmac671"). Then, execute

  cmacxxx fileName

where "fileName" is the name of the script file describing the simulation. For
other UNIX platforms, note that the CalC executable may have a different name
(e.g. "ccyg671" for cygwin, or "calc" if you compile the source code yourself)


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

o UNIX platforms (including cygwin and Mac OS X):

It is recommended that you use the "xmgr" or "xmgrace" graphics applications,
which allow viewing the simulation results in real time. This involves piping
the program output into xmgr or xmgrace, by executing one of the following:

  calc scriptFileName | xmgr -pipe 
  calc scriptFileName | xmgrace -pipe

where "scriptfileName" is the name of the simulation script. In this case the 
script file should contain an instruction "plot.method xmgr".

  Xmgrace can be downloaded from
  ftp://plasma-gate.weizmann.ac.il/pub/grace/

  Xmgr can be downloaded from 
  ftp://plasma-gate.weizmann.ac.il/pub/xmgr4/

o Windows platforms: 

Use "mute" and/or "binary" plots to produce data files, which can then be plotted
using MATLAB. Alternatively, install "cygwin" (http://www.cygwin.com"), a Linux
emulator for Windows, which would enable you to install the unix version of
CalC, and use its xmgrace graphics capabilities (see above).

o MATLAB integration:

Like any system command, CalC can be launched from MATLAB simply by executing
the command "system('.\cwinxx.exe script.par')". See CalC demo script web page for 
further details.
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
