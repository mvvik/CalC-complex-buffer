/*****************************************************************************
 *
 *                       Calcium Calculator (CalC)
 *              Copyright (C) 2001-2019 Victor Matveev
 *   
 *                               fplot.h
 *
 *  Defines all plot objects, including mute (file dump) plots and xmgr/xmgrace pipes
 *
 *                           Class Hierarchy
 *                           ^^^^^^^^^^^^^^^
 *            /- StateDump                            
 *           |                                      
 *           |-- MutePointPlot     /-- FieldPlot2D ------ MutePlot2D 
 *           |                    |                 
 *  PlotObj -|-- FieldPlotObj ----|--- FieldDump 
 *           |                    |
 *           |                    |--- FieldDumpT ("binary", or "4D" plot)
 *           |                    |                    
 *           |                     \-- FieldPlot1D \----- MutePlot1D 
 *            \                                    |   
 *             \_____ XmgrPointPlot                |    
 *             /                                   | 
 *            /                                     \_____ 
 *  XmgrPlot --------------------------------------------- Xmgr1Dplot
 *
 *
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

#ifndef CALC_FPLOT_H_included
#define CALC_FPLOT_H_included
#define _CRT_SECURE_NO_DEPRECATE

#define METHOD_XMGR 1
#define METHOD_MUTE 2

#define XMGR_DIVISIONS 6

#define xmgr_dxmargin 0.12
#define xmgr_dymargin 0.12

#define xmgr_tmargin 0.1
#define xmgr_bmargin 0.1

#define xmgr_lmargin 0.1
#define xmgr_rmargin 0.06

#define dx_legend 0.2
#define dy_legend 0.1

#define title_color 4 


double format_double(double number0, double number1, int divisions, int sigdigits, char *fstr);
double format_digits(double number, int &digits);

class SimulationObj; // forward class declaration 

//*******************************************************************************
//                 C L A S S   P L O T O B J E C T
//*******************************************************************************

class PlotObj
{
protected:
public:

  static int    UPDATE_STEPS;
  static int    UPDATE_STEPS_1D;
  static int    UPDATE_STEPS_BINARY;
  static double UPDATE_ACCURACY;

  int   plot_height, win_height, plot_width,  win_width;
  char  win_title[255];
  long  winid;

  char     x_label[80];
  double   x_value;

  char  log_plot, mute;

  double   fmin, fmax;
  double   *fptr;

  PlotObj(double *ptr, char islog, const char *txt = "", const char *label = "")
    {
    fptr = ptr; log_plot = islog;
    if (txt) strcpy(win_title,txt); else win_title[0] = 0;
    if (label) strcpy(x_label,label); else x_label[0] = 0;
    x_value = fmax = fmin = 0.0;
    };

  virtual ~PlotObj() { };

  virtual void draw() = 0;
  virtual void redraw() = 0;
  void    draw_grid(double, double, double, double, char);

  virtual double get_value(long ind = 0)
   {
   double v = fptr[ind];
   if (!log_plot) return v;
      else return ( (v <= 0.0) ? 0.0 : log(v)/log(10.0) );
   };

};

//*******************************************************************************
//                   C L A S S    S T A T E    D U M P
//*******************************************************************************

class StateDump : public PlotObj
{
protected:

  bool   complete;
  double exportTime;
  char   fileName[255];
  class  SimulationObj *Sim;

public:

 FILE   *file;

 StateDump( SimulationObj *sim, double t,
            const char *fname = 0, const char *txt = 0) : PlotObj(0, 0, txt) {
    strcpy(fileName, fname); 
    Sim = sim;
    exportTime = t;
    complete = false;
    };

 ~StateDump() { }

 void    draw();
 void    redraw() { draw(); }
};

//*******************************************************************************
//               C L A S S   M U T E   P O I N T   P L O T
//*******************************************************************************

class MutePointPlot : public PlotObj
{
protected:

  double *Time;
  double fold;
  double tscale;
  char   fileName[255];
  double x_temp, f_temp;
  long   counter, bufferSize;
  double *Tbuffer, *Ybuffer;

public:

 FILE   *file;

 MutePointPlot(double *ptr, double *tptr, char islog,  double t, 
               const char *fname = "", const char *txt = "");

 ~MutePointPlot() { fclose(file); delete [] Tbuffer; delete [] Ybuffer; };

 void    draw();
 void    redraw()  { x_value = 0; draw(); }
 void    pushValue(double t, double y);
};


//*******************************************************************************
//                 C L A S S   F I E L D P L O T O B J E C T
//*******************************************************************************

class FieldPlotObj : public PlotObj
{
protected:
public:

 FieldObj *field; 
 long     field_index;

 FieldPlotObj(FieldObj *f, char islog, const char *txt = 0, const char *label = 0):
   PlotObj((*f)(), islog, txt, label)  { field = f; field_index = 0; };

 double get_value(long ind = 0);

};
//*******************************************************************************

class FieldDump : public FieldPlotObj
{
protected:

  bool   complete;
  double exportTime;
  char   fileName[255];

public:

 FILE   *file;

 FieldDump(FieldObj *f, double t, const char *fname = 0, const char *txt = 0) : FieldPlotObj(f, 0, txt)
    {
    strcpy(fileName, fname);
	sprintf(win_title, "%s (dump into %s at time %g)", f->ID, fileName, t);
    exportTime = t;
    complete = false;
    }

 // ~FieldDump() { }

 void    draw()  {
   if (field->Time < exportTime) { complete = false; return; }
   if ( complete ) return;
   complete = true;
   field->Export(fileName);
   };

 void    redraw() { draw(); }
};
//*******************************************************************************

class FieldDumpT : public FieldPlotObj
{
protected:

  double timeBetweenSaves;
  double oldTime, newTime;
  char   fileName[255];

public:

 FILE   *file;

 FieldDumpT(FieldObj *f, double total, const char *fname = 0, const char *txt = 0) : FieldPlotObj(f, 0, txt)
    {
		strcpy(fileName, fname);
		sprintf(win_title, "%s (binary plot into file %s)", field->ID, fileName);
		timeBetweenSaves = total / double(UPDATE_STEPS_BINARY);
		oldTime  = field->Time;
		newTime  = oldTime + timeBetweenSaves;
        FILE *ff = fopenAssure(fileName, "wb", "continuous field dump", "header initialization");
        header(ff);
		fclose(ff);
		dump();
    }
   
   //**********************************************

	 void  header(FILE *ff) {
								fwrite( (void *)(&GEOMETRY),     sizeof(int),   1, ff);
								fwrite( (void *)(&field->xsize), sizeof(int),   1, ff);
		if (DIMENSIONALITY > 1)	fwrite( (void *)(&field->ysize), sizeof(int),   1, ff);
		if (DIMENSIONALITY > 2) fwrite( (void *)(&field->zsize), sizeof(int),   1, ff);

								fwrite( (void *)(field->xcoord), sizeof(double), field->xsize, ff);
		if (DIMENSIONALITY > 1)	fwrite( (void *)(field->ycoord), sizeof(double), field->ysize, ff);
		if (DIMENSIONALITY > 2) fwrite( (void *)(field->zcoord), sizeof(double), field->zsize, ff);
	 }

   //**********************************************

	 void  chuckHeader(FILE *ff, void *buffer) {
								fread( (void *)buffer, sizeof(int),    DIMENSIONALITY + 1, ff);
								fread( (void *)buffer, sizeof(double), field->xsize,       ff);
		if (DIMENSIONALITY > 1)	fread( (void *)buffer, sizeof(double), field->ysize,       ff);
		if (DIMENSIONALITY > 2) fread( (void *)buffer, sizeof(double), field->zsize,       ff);
	 }
	   
   //**********************************************

	 void  dump() {
	   FILE *ff = fopenAssure(fileName, "ab", "continuous field dump", "binary write");

	   fwrite( (void *)&field->Time, sizeof(double), 1,           ff);
	   fwrite( (void *) field->elem, sizeof(double), field->size, ff);

	   fclose(ff);
	 }

 //**********************************************

	 void  draw();
  
 //**********************************************

  ~FieldDumpT() { dump(); }
 
  void    redraw() { draw(); }
};

//*******************************************************************************
//                      C L A S S   1 D - P L O T
//*******************************************************************************

class FieldPlot1D : public FieldPlotObj
{
protected:

  long    incr, num;
  double  *grid, *coord;

public:

 FieldPlot1D(FieldObj *f, char islog, const char *dir, double coord1, double coord2, 
             const char *label = "time=");

 virtual void    get_range();
};


//*******************************************************************************

class MutePlot1D : public FieldPlot1D
{
public:
 
 FILE   *file;
 char   fileName[255];
 double total, exportTime;
 double tscale;
 bool   complete;

 long   counter, bufferSize;
 double *Tbuffer,  *Ybuffer;

 MutePlot1D(FieldObj *f, char islog, const char *dir, double coord1, double coord2, double T,
            const char *fname = 0, const char *txt = 0) : FieldPlot1D(f, islog, dir, coord1, coord2)  
 {
   exportTime = total = T;
   x_value = -1e-12;
   complete = false;
   tscale = double(UPDATE_STEPS_1D) / total;
   strcpy(fileName, fname);
   file = 0;

   bufferSize = UPDATE_STEPS_1D;
   Tbuffer = new double[UPDATE_STEPS_1D + 1];
   Ybuffer = new double[num * (UPDATE_STEPS_1D + 1)];
   counter = 0;
 };

 ~MutePlot1D() { delete [] Tbuffer;  delete [] Ybuffer;  if (file) fclose(file);  };

 void    draw();
 void    pushRow(double t, double *F, int columns = 3);
 void    redraw() { /* x_value = 0.0; */ draw(); }; 

};


//*******************************************************************************
//                C L A S S   2 D   F I E L D   P L O T
//*******************************************************************************

class FieldPlot2D : public FieldPlotObj
{
protected:

  long    incr1, incr2, num1, num2;
  double  *grid1, *grid2;
  double  *coord1, *coord2;

public:

 FieldPlot2D(FieldObj *f, char islog, const char *dir, double coord, const char *label = "time=");

 virtual void  get_range();
};

//*******************************************************************************

class MutePlot2D : public FieldPlot2D
{
public:

	bool   complete;
	double exportTime;
	char   fileName[255];

	MutePlot2D(FieldObj *f, char islog, const char *dir, double coord, double T,
		const char *fname = 0, const char *txt = 0) :
	FieldPlot2D(f, islog, dir, coord)
	{
		exportTime = T;
		strcpy(fileName,fname);
	};

	~MutePlot2D() {  };

	void    draw() 
	{  
		if (field->Time < exportTime) { complete = false; return; }
		if ( complete ) return;
		complete = true;
		long ind = field_index;
		FILE *file = fopenAssure(fileName, "w", "plot", win_title);
		for (int j = 0; j < num2; j++) {
			for (int i = 0; i < num1; i++) {
				ind = field_index + j * incr2 + i * incr1;
				fprintf(file,"%.6g %.6g %.9g \n", coord1[i], coord2[j], get_value(ind));
			}
			fprintf(file,"\n");
		}
		fclose(file);
	}

	void    redraw() { draw(); }; 
};

//*******************************************************************************
//           B A S E   C L A S S   X M G R   P I P E   P L O T
//*******************************************************************************
// xmgr colors: 0-white 1-black 2-red 3-green 4-blue 5-yellow 6-brown 7-gray 
// 8-violet 9-cyan 10-magenta 11-orange 12-indigo 13-maroon 14-turquise 15 -green4

class XmgrPlot
{
protected:

  static int XMGR_STEPS;
  static int graph_num;
  static int graph_count;

  static int    rows;
  static int    cols;
  static int    divs;
  static double xscale;
  static double yscale;
 
  static FILE *file;

  int    steps;

  int    set_num;
  int    set_count;
  int    graph_id;

  int    view_left, view_right, view_top, view_bottom;

  double first, last;

  bool   logScale;

public:

 static void init(TokenString &, int gnum, FILE *f = stdout, int = 0, int = 0);
 static void processStrings(TokenString &);

 XmgrPlot(int sets, char is_log, const char *xlabel = "time");
 ~XmgrPlot() {};

 int setColor(int);
};

//*******************************************************************************
//               C L A S S   X M G R   P I P E   P L O T
//*******************************************************************************
// xmgr colors: 0-white 1-black 2-red 3-green 4-blue 5-yellow 6-brown 7-gray 
// 8-violet 9-cyan 10-magenta 11-orange 12-indigo 13-maroon 14-turquise 15 -green4

class XmgrPointPlot : public XmgrPlot, public PlotObj
{
protected:

  double **tptrs, **fptrs;
  double *x_old;
  double *f_old;
  double *x_temp;
  double *f_temp;

  double **tArray, **fArray;
  int    *pointCount, *maxPoints;

  double tScale;

public:

 XmgrPointPlot(int sets, char islog, double T = 0, const char *xlabel = "time");

 void set_set(double *p, double *t, const char *txt = 0);

 ~XmgrPointPlot() { for (int i=0; i<set_num; i++) { delete [] tArray[i]; delete [] fArray[i]; }
                    delete [] tArray; delete [] fArray; delete [] pointCount; delete [] maxPoints;
                    delete [] tptrs;  delete [] fptrs;  delete [] x_old;      delete [] f_old; 
					delete [] x_temp; delete [] f_temp; 
                   } 

 void    get_range();
 void    draw();
 void    redraw()  { for (int i = 0; i < set_num; i++) x_old[i] = x_temp[i] = 0.0; draw(); }
 void    truncate(int set);
 void    addMemory(int);
};


//*******************************************************************************
//            C L A S S   1 D   X M G R   P I P E   P L O T  
//*******************************************************************************

class Xmgr1Dplot : public XmgrPlot, public FieldPlot1D
{
protected:

 double setTime[30];

public:

 Xmgr1Dplot(FieldObj *f, char islog, const char *dir, double coord1, double coord2, double T,  int m = 5);

 void    draw();
 void    redraw() { /*x_value = 0.0;*/ draw(); }; 
 int     setColor(int);
};

//*******************************************************************************
//            C L A S S   1 D   X M G R   P I P E   P L O T  
//*******************************************************************************

class Xmgr2Dplot : public XmgrPlot, public FieldPlot2D
{
protected:

 double setTime[30];
 double Theta;
 double FscaleFactor;

 double cosine, sine;
 int  nBinsTotal, nBinsX, xBin;
 double *Ybuffer;
 double Xmax, Ymax, Xscale, Yscale, binScale, binMult;

public:

   Xmgr2Dplot(FieldObj *f, char islog, const char *dir, double coord, double T, double angle=1.1, double factor = 0.4, int nBins = 400);
  ~Xmgr2Dplot() { delete [] Ybuffer; }

 void    draw();
 void    redraw() { /*x_value = 0.0;*/ draw(); }; 
 int     setColor(int i = 0) { return 2; }

};

//*******************************************************************************
//                      C L A S S   P L O T   A R R A Y
//*******************************************************************************

class PlotArray
{
 protected:

 public:

  PlotObj **array;
  int     plot_num;
  int     count;
  char    method;

  PlotArray(class SimulationObj  &);

  int get_plot_num(TokenString &);

  PlotArray(int n)   { count = 0;
                       plot_num = n;
                       if (n) array = new PlotObj *[n]; }
  ~PlotArray();

  void set_plot(PlotObj *plot)  {
     if (count >= plot_num) 
        throw makeMessage("in PlotArray: cannot set plot #%d: total number = %d", count, plot_num); 
     else array[count++] = plot;
     }

  void draw_plot(int n)   { if (n < count) array[n]->draw();   }
  void redraw_plot(int n) { if (n < count) array[n]->redraw(); }
  void draw_all()         { if (count) for (int i=0; i<count; i++) array[i]->draw();   }
  void redraw_all()       { if (count) for (int i=0; i<count; i++) array[i]->redraw(); }

};

//*****************************************************************************
//                      S T A T U S      W I D G E T
//*****************************************************************************


class RunStatusString {

  int    statusLength, widgetLength, widgetCount;
  char   *widgetStr, *wholeStr;
  char   *timeStr, *extraStr;
  double *timePtr, *extraPtr;
  double total, t0;
  double realTime;

 public:

  RunStatusString(int length, const char *s1, double *p1, double T, const char *s2=0, double *p2=0);
  ~RunStatusString();
  
  void reset(double T = 0);
  void update(char c = 0);
};



//*******************************************************************************

#endif
