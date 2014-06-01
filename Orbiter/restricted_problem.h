#include "ode.h"
#include <vector>
#include <valarray>
#include <ctime>
#include <sstream>
//#define double long double
using namespace std;

void integrate(double C0, double x0, double dt, double tmax, vector<double>& xaxis, vector<double>& yaxis, int whichIndicator, double tol, double &indicator, double &CPU_time, bool evolution);
void findmax(vector<double> &sol);
