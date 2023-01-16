******************************************************************************
## README: Calcium Calculator (CalC)
- Developed by [Victor Matveev](http://www.victormatveev.org), Department of Mathematical Sciences, NJIT
- Initial development (2001-2003): Laboratory of Biological Modeling, NIDDK, NIH
- Email any bug reports (make sure to include your script) to: matveev@njit.edu
- Mirror maintained at http://www.calciumcalculator.org
******************************************************************************
CalC ("Calcium Calculator") is a modeling tool for simulating intracellular calcium diffusion and buffering. CalC solves continuous reaction-diffusion PDEs describing the entry of calcium into a volume through point-like channels, and its diffusion, buffering and binding to calcium receptors. CalC uses a variation of the Alternating Direction Implicit (ADI) finite difference method, which is quite CPU-time efficient, and accurate to 2nd order in time and space. Time-step is varied adaptively during the simulation. Other main features are:
- CalC is platform-independent (Windows, OS X, Linux, cygwin, etc.)
- CalC is operated by a simple script language (with optional flow-control functionality).
- CalC is easily combined with MATLAB without any special modifications (see below).
- CalC allows simulations in any geometry: cartesian 3D, 2D or 1D, polar, spherical, cylindrical, conical, etc.
- CalC allows an arbitrary number of calcium buffers, with a single or two calcium binding sites per molecule
- CalC scripts can integate ordinary differential equations as well, e.g. to model calcium dependent exocytosis.
- CalC results can be viewed in real time using **xmgrace** or **freeglut** libraries (see below)


If you use CalC in your published work, please cite [2002 Biophys J article](https://pubmed.ncbi.nlm.nih.gov/12202362/) article, and please send me a reference for inclusion in the [CalC publication list](https://web.njit.edu/~matveev/calc/calc_pub.html) upon publication. CalC is provided on an as-is basis, but I will respond to any bug reports or technical questions.
******************************************************************************
#### License
CalC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

CalC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should see a copy of the GNU General Public License in this repository.  If not, [follow this URL](<https://www.gnu.org/licenses/>)
******************************************************************************
### Executables

Executables for the latest versions of **Windows** and **macOS** are contained in the **executables** folder of this repository.  Note that the name of the executable file varies between different systems (you can rename it as you like, obviously). In this document, executable is referred to by the name **calc**.  If the executable is not working on your OS, follow **compilation** instructions below. Otherwise, proceeed to the **execution** section of this README file.

******************************************************************************
### Compilation

NOTE: you can find all platform-specific definitions in the short header file called **PlatformSpecific.h**, contained within the source folder.

##### 1) Windows:   

Use MSDN Visual Studio to compile the code. The Visual Studio "Project" file **CalC.vcxproj** contains all the necessary dependencies (it resides within the source folder). Alternatively, you can install the Linux emulator [cygwin](http://www.cygwin.com") on your Windows machine, and follow the UNIX installation instructions below. 

##### 2) macOS

Install XQuartz and Xcode. The latter will have a C++ compiler such as g++, and some kind of make utility to complile the code. Then, run **make**. That is all.
   
##### 3) Non-Mac UNIX / LINUX / CygWin

Unless you have [freeglut](https://freeglut.sourceforge.net/) installed, replace **Makefile** with **Makefile.no_glut** (and remove the extension **.no_glut**). Then, run **make**, as usual. This will compile a version of the program without run-time OpenGl/GLUT/FreeGlut graphics.

Alternatively, install [freeglut](https://freeglut.sourceforge.net/) and change the linker option in the **Makefile** to properly link your installation of GLUT/FreeGlut library. See also **PlatformSpecific.h** header file in the source directory.

******************************************************************************
### Execution (all platforms)

Windows OS only, CalC versions **x.10.1** or higher: make sure that the provided **freeglut.dll** is present in the folder containing your calc executable.

If a model script does not use command-line parameters, simply click on the executable and enter the script file name when prompted (the script file should reside in the same directory). Alternatively, name your CalC script **DefaultScripts.txt**, and it will be executed automatically after you lauch the executable.

If your script utilizes command-line parameters, then open a shell (for Windows, launch Start Menu -> Run -> cmd), go to the directory where the executable **calc** is residing, and type 

    calc filename  parList
 
where **calc** is the name of the executable (replace with correct executable name -- see **executables** folder or compilation instructions above), **fileName** is the name of the script file describing the simulation, and **parList** is an optional space-separated list of command-line parameters (see [manual](http://web.njit.edu/~matveev/calc/manual.html#pars)).

In order to monitor program output and error messages, include the statement **verbose = 4** (or higher verbosity level) in your script: this will prevent CalC from auto-terminating upon completing the simulation.

******************************************************************************
### Data Output

##### File plots (all platforms): 
Simulation results can be saved to files in real time using [mute plot](http://web.njit.edu/~matveev/calc/manual.html#method_mute) statements, ASCII or binary files are produced that are readable by any graphics-capable language such as MATLAB (Mathworks, Inc). See [demo scripts](https://web.njit.edu/~matveev/calc/scripts.html) and refer to the [manual](http://web.njit.edu/~matveev/calc/manual.html#method_mute) for details.

The [binary](http://web.njit.edu/~matveev/calc/manual.html#binary) plot type allows to save an entire concentration field at several time points during the simulation, and can be read and displayed using MATLAB via scripts provided in the **examples** directory and on the [demo script page](https://web.njit.edu/~matveev/calc/scripts.html)

#####  Real-time OpenGL plots:

Include command [plot.method gl](https://web.njit.edu/~matveev/calc/manual.html#method_gl) within your script to make real-time variable plots (or 1D and 2D concentration plots) in your OS window. On Windows OS, make sure that the **freeglut.dll** dynamic library is in the same folder as your executable (it is provided in the repository). On macOS, GLUT is preinstalled (but depricated). On other platforms, you have to have GLUT/freeglut installed on your computer, and change the linker directive in the Makefile appropriately.

NOTE: Graphics buffer flushing appears flaky with freeglut, and sometimes graphs are not updated until the entire script runs to completion. Further, including gl plots may significantly slow down the script execution. Therefore, gl plots are useful mostly for initial script debugging. 

#####  Xmgrace plots (UNIX platforms):

You can use the [xmgrace](https://plasma-gate.weizmann.ac.il/Grace/) graphics application instead of **freeglut**, which also allows viewing the simulation results in real time. This involves piping the program output into 
xmgrace, by executing the following command:

    calc scriptFileName | xmgrace -pipe

where **calc** is the name of your CalC executable, **scriptfileName** is the name of the simulation script. In this case the script file should contain an instruction **plot.method xmgr**.

A [**homebrew**](https://formulae.brew.sh/formula/grace) installation of xmgrace is quite easy to perform. 

******************************************************************************
###  MATLAB integration:

Like any system program, **CalC** can be launched from **MATLAB** (Mathworks, Inc) simply by executing the command

    system('.\calc ScriptFileName')
    
You can then collect the data by reading the **mute plot** or **binary** files (see above about data output). See CalC demo script web page for further details.

******************************************************************************
### Documentation

A hyper-text CalC script syntax manual can be found at [this URL](https://web.njit.edu/~matveev/calc/manual.html)
Also, examine the included commented example script files (file extension ".par"): they reside in the **examples** folder. Example script files are also available at http://web.njit.edu/~matveev/calc/scripts.html
******************************************************************************
