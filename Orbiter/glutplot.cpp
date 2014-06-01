#include "glutplot.h"
//using namespace std;

double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300, dx, dy;
bool dots_or_lines, y_logarithmic, x_logarithmic;
vector<vector<double> > xaxis, yaxis;
char *xsign, *ysign;
double xTab, yTab;
const int ncolors=1;
const double colors[][3] = {{0.7,0.1,0.1}, {0.1,0.7,0.7}, {0.1,0.1,0.7}, {0.7, 0.1, 0.7}, {0.7, 0.7, 0.1}};
//const double colors[][3] = {{0.7,0.1,0.1}, {0.1,0.7,0.7}, {0.7, 0.7, 0.1} };
void glWrite(double x, double y, string ss)
{
	glRasterPos2f(x,y);
	for (int i=0; i<ss.length(); i++)
	{
		glutBitmapCharacter(font, ss[i]);
	}
}

void glWriteNumber(double x, double y, double mydouble)
{
	ostringstream ss(ostringstream::out);
	//ss.precision(2);
	ss << mydouble;
	glWrite(x,y,ss.str());
}

double scr(double z)
{
	if (y_logarithmic) return ymin + (ymax-ymin)*log10(z/ymin)/log10(ymax/ymin);
	else return z;
}

double x_scr(double z)
{
	if (x_logarithmic) return xmin + (xmax-xmin)*log10(z/xmin)/log10(xmax/xmin);	
	else return z;
}
void calculateX()
{
	double temp = log10(xmax-xmin), exponent, fractpart =  modf(temp, &exponent);
	exponent=floor(temp);
	if ((fabs(fractpart)>0.65 && temp<0) || (fabs(fractpart)<0.4 && temp>0)) exponent-=1;
	
	dx=pow (10, exponent );
	if (floor((xmax-xmin)/dx)<=15)  dx/=10;
	double x = dx*ceil( xmin / dx), length;

	for(int i=1; x <= xmax; i++) 
	{
		glColor3f(0.7,0.7,0.7);
		glLineWidth(1); 
		length=1;
		double intpart, mod = modf(x/(dx*10), &intpart);
	//	cout<<"**xx**"<<endl<< x << " "<< dx <<" "<< x/dx << " " <<fabs( 1- mod) <<endl<< "**xx**" << endl; 
		if ( fabs(mod) < 1e-6 || fabs(1-mod) < 1e-6 || fabs(1+mod) < 1e-6) {
			
			glColor3f(0,0.5,0);
			if (x!=0) {
				glWriteNumber(x - xTab/4, ymin - 3*yTab/4 , dx*int(x/dx+0.5* (x/fabs(x)) ));
				glWriteNumber(x - xTab/4, ymax + yTab/2 , dx*int(x/dx+0.5* (x/fabs(x)) ));
			}
			else {
				glWriteNumber(x - xTab/4, ymin - 3*yTab/4 , 0);
				glWriteNumber(x - xTab/4, ymax + yTab/2 , 0);
			 }
			glColor3f(0.5,0.5,0.5); 
			glLineWidth(2); 
			length=1.5;
		}
		
		glBegin(GL_LINES);
			glVertex2f(x, ymax+2*yTab/8);
			glVertex2f(x, ymax+(2-length)*yTab/8);
		glEnd();
		glBegin(GL_LINES);
			glVertex2f(x, ymin-2*yTab/8);
			glVertex2f(x, ymin-(2-length)*yTab/8);
		glEnd();
		
		x+=dx;
	}		
}
void calculateLogX()
{
	dx = pow(10, floor(log10(xmin)));
	double x = dx*ceil(xmin / dx)-dx, length;
		
	for(int i=1; x <= xmax; i++) 
	{
		glColor3f(0.7,0.7,0.7);
		glLineWidth(1); 
		length=1;
		double intpart, mod = modf(log10(x), &intpart);
		//cout<<"***"<<endl<< x << " "<< dx <<" "<< x/dx << " " <<fabs(mod) <<endl<< "***" << endl; 
		if ( fabs(mod)<1e-6) {
			
			glColor3f(0,0.5,0);
			glWriteNumber(x_scr(x) - xTab/4, ymin - 3*yTab/4 , x);
			glWriteNumber(x_scr(x) - xTab/4, ymax + yTab/2 , x);
						 
			glColor3f(0.5,0.5,0.5); 
			glLineWidth(2); 
			length=1.5;
			if (x_logarithmic && x!=dx) dx*=10;
		}
		
		glBegin(GL_LINES);
			glVertex2f(x_scr(x), ymax+2*yTab/8);
			glVertex2f(x_scr(x), ymax+(2-length)*yTab/8);
		glEnd();
		glBegin(GL_LINES);
			glVertex2f(x_scr(x), ymin-2*yTab/8);
			glVertex2f(x_scr(x), ymin-(2-length)*yTab/8);
		glEnd();
		
		x+=dx;
	}
}

void calculateY()
{
	double temp = log10(ymax-ymin), exponent, fractpart =  modf(temp, &exponent);
	exponent=floor(temp);
	if ((fabs(fractpart)>0.65 && temp<0) || (fabs(fractpart)<0.4 && temp>0)) exponent-=1;
	
	dy=pow (10, exponent );
	//ytick= floor((ymax-ymin)/dy);
	if (floor((ymax-ymin)/dy)<15)  dy/=10;
	double y = dy*ceil( ymin / dy), length;
	cout<<y <<" "<<ymax<< endl;
	for(int i=1; y <= ymax; i++) 
	{
		glColor3f(0.7,0.7,0.7);
		glLineWidth(1); 
		length=1;
		double intpart, mod = modf(y/(dy*10), &intpart);
		//cout<<"***"<<endl<< y << " "<< dy <<" "<< y/(dy*10) << " " << (mod) <<endl<< "***" << endl; 
		if ( fabs(mod) < 1e-6 || fabs(mod-1) < 1e-6 || fabs(mod+1) < 1e-6) {
			
			glColor3f(0,0.5,0);
			if (y!=0) {
					glWriteNumber(xmin-xTab*0.8, y - yTab/20, dy*int(y/dy+0.5* (y/fabs(y)) )); 
					glWriteNumber(xmax+xTab*0.4, y - yTab/20, dy*int(y/dy+0.5* (y/fabs(y)) )); 
			}
			else {
					glWriteNumber(xmin-xTab*0.8, y - yTab/20, 0); 
					glWriteNumber(xmax+xTab*0.4, y - yTab/20, 0); 
			}			 
			glColor3f(0.5,0.5,0.5); 
			glLineWidth(2); 
			length=1.5;
		}
		
		glBegin(GL_LINES);
			glVertex2f(xmin-2*xTab/8, y);
			glVertex2f(xmin-(2-length)*xTab/8, y);
		glEnd();
		glBegin(GL_LINES);
			glVertex2f(xmax+2*xTab/8, y);
			glVertex2f(xmax+(2-length)*xTab/8, y);
		glEnd();
		
		y+=dy;
	}		
}
void calculateLogY()
{
	dy = pow(10, floor(log10(ymin)));
	double y= dy*ceil(ymin / dy)-dy, length;
	//cout<<"***"<<endl<< y << " "<< dy <<" "<< y/dy << " " << ymin <<endl<< "***" << endl; 	
	for(int i=1; y <= ymax; i++) 
	{
		glColor3f(0.7,0.7,0.7);
		glLineWidth(1); 
		length=1;
		double intpart, mod = modf(log10(y), &intpart);
		//cout<<"***"<<endl<< y << " "<< dy <<" "<< y/dy << " " <<fabs(mod) <<endl<< "***" << endl; 
		if ( fabs(mod) < 1e-6 || fabs(mod-1) < 1e-6 || fabs(mod+1) < 1e-6) {
			
			glColor3f(0,0.5,0);
			glWriteNumber(xmin-xTab*0.8, scr(y) - yTab/20, y); 
			glWriteNumber(xmax+xTab*0.4, scr(y) - yTab/20, y); 
						 
			glColor3f(0.5,0.5,0.5); 
			glLineWidth(2); 
			length=1.5;
			if (y_logarithmic && y!=dy) dy*=10;
		}
		
		glBegin(GL_LINES);
			glVertex2f(xmin-2*xTab/8, scr(y));
			glVertex2f(xmin-(2-length)*xTab/8, scr(y));
		glEnd();
		glBegin(GL_LINES);
			glVertex2f(xmax+2*xTab/8, scr(y));
			glVertex2f(xmax+(2-length)*xTab/8, scr(y));
		glEnd();

		
		y+=dy;
	}
}

bool outOfBorders(double x, double y)
{
	if (x<xmin || x>xmax || y < ymin || y>ymax) 
	return true; 
	else return false;
}
void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);	
	
	glLineWidth(2.0);
	glPointSize(1);
	
	int width = glutGet(GLUT_SCREEN_WIDTH) - 300;
	for (int k=0; k<xaxis.size(); k++) {
		glColor3f(colors[k % ncolors][0],colors[k % ncolors][1],colors[k % ncolors][2]);
		double magn = pow(10, log10(xaxis[k].size()-2)/width), old_i = 0;

		for (int i = 1; i<xaxis[k].size()-2; (x_logarithmic && i>int( 1+1/(magn-1) ) ) ? i=int(i*magn) : i++)
		{
			
			if (!outOfBorders(xaxis[k][i],yaxis[k][i])) {
			
			if (dots_or_lines)
			{
				glBegin(GL_LINES);
					glVertex2f(x_scr(xaxis[k][old_i]), scr(yaxis[k][old_i]));
					glVertex2f(x_scr(xaxis[k][i]), scr(yaxis[k][i]));
				glEnd();
			}
			else
			{
				glBegin(GL_POINTS);
					glVertex2f(x_scr(xaxis[k][i]), scr(yaxis[k][i]));
					//cout<<xaxis[k][i] <<" " << yaxis[k][i]<< endl;
				glEnd();
			}
			}
			old_i=i;	
		}
	}

	glLineWidth(1);
	glColor3f(0.7,0.7,0.7);
	glBegin(GL_LINES);
			glVertex2f(xmax+2*xTab/8, ymax+2*yTab/8);
			glVertex2f(xmax+2*xTab/8, ymin-2*yTab/8);
	glEnd();
	glBegin(GL_LINES);
			glVertex2f(xmax+2*xTab/8, ymin-2*yTab/8);
			glVertex2f(xmin-2*xTab/8, ymin-2*yTab/8);
	glEnd();
	glBegin(GL_LINES);
			glVertex2f(xmin-2*xTab/8, ymin-2*yTab/8);
			glVertex2f(xmin-2*xTab/8, ymax+2*yTab/8);
	glEnd();
	glBegin(GL_LINES);
			glVertex2f(xmin-2*xTab/8, ymax+2*yTab/8);
			glVertex2f(xmax+2*xTab/8, ymax+2*yTab/8);
	glEnd();	
	
	(x_logarithmic) ? calculateLogX() : calculateX();
	(y_logarithmic) ? calculateLogY() : calculateY();
	
	glColor3f(0.5,0.5,0.5);
	glWrite(xmax+xTab/3,ymin-3*yTab/4,xsign);
	glWrite(xmin - xTab*0.8,ymax+3*yTab/4, ysign);
	
    glFlush();
}

void plot(int argc, char**argv, const vector< vector<double> > &_x, const vector<vector<double> > &_y, char* _xsign, char* _ysign, bool lines, bool _xlog, bool _ylog, double _xmin, double _xmax, double _ymin, double _ymax, bool automaticBorder)
{
glutInit(&argc,argv);
glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA);
glutInitWindowPosition(50, 25);
glutInitWindowSize(1300, 700);
glutCreateWindow("Green window");
glClearColor(1.0,1.0,1.0,1.0);

dots_or_lines = lines;
xaxis = _x;
yaxis = _y;
xsign = _xsign;
ysign = _ysign;
x_logarithmic = _xlog;
y_logarithmic = _ylog;

if (automaticBorder)
{
	for (int i=0; i<_x.size(); i++) {
		int size = _x[i].size() - 2;
		if (_x[i][size]   < xmin) xmin = _x[i][size];
		if (_x[i][size+1] > xmax) xmax = _x[i][size+1];
		if (_y[i][size]   < ymin) ymin = _y[i][size];
		if (_y[i][size+1] > ymax) ymax = _y[i][size+1];
	}
}
else
{
	xmin=_xmin;
	xmax=_xmax;
	ymin=_ymin;
	ymax=_ymax;
}
cout<< xmin << " " << xmax << " " << ymin<<" " << ymax << endl;
xTab = 150*(xmax-xmin)/glutGet(GLUT_SCREEN_WIDTH);
yTab = 150*(ymax-ymin)/glutGet(GLUT_SCREEN_HEIGHT);
glOrtho(xmin-xTab,xmax+xTab,ymin-yTab,ymax+yTab,-1.0,1.0);
glutDisplayFunc( draw );
glutMainLoop();
}

