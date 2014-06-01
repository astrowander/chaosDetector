#include <GL/glut.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <limits>	
#include <cmath>
//#define double long double

#define font GLUT_BITMAP_HELVETICA_18
using namespace std;

void plot(int argc, char**argv, const vector< vector<double> > &_x, const vector<vector<double> >  &_y, char* _xsign, char* _ysign, bool lines, bool _xlog, bool _ylog, double _xmin, double _xmax, double _ymin, double _ymax, bool automaticBorder);
