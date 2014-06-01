#include <cmath>
#include <valarray>
#include <iostream>
//#define double long double
using namespace std;
typedef valarray<double> (*pt2Func) (double,const valarray<double> &);

//valarray<double> euler(double &t, const valarray<double> &x, double dt, pt2Func eq);
//void rkutta4(double &t, valarray<double> &x, double dt, pt2Func eq);
//valarray<double> merson5(double &t, const valarray<double> &x, double &dt, double ode_eps, pt2Func eq);
//void dopri5(double &t, valarray<double> &x, double &dt, double ode_eps, pt2Func eq);
//void butcher6(double &t, valarray<double> &x, double &dt, pt2Func eq);
void fehlberg8(double &t, valarray<double> &x, double &dt, double hmin, double tolerance, pt2Func eq);
void fastFehlberg8(double &t, valarray<double> &x, double &dt, double tolerance, pt2Func eq, valarray<double> *k);
void modFastFehlberg8(double &t, valarray<double> &x, double &h, double hmin ,double tolerance, pt2Func eq, valarray<double> *k);
void verner9(double &t, valarray<double> &x, double &dt, double hmin, double tolerance, pt2Func eq);
//void everhart15(double &t, valarray<double> &x, double &dt, double tol, int ni, int ns, int nf, pt2FortranSub force, bool start); 
//valarray<double> dopri8(double &t, const valarray<double> &x, double &dt, pt2Func eq);
