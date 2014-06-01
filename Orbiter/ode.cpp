#include "ode.h"
#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )
/*valarray<double> euler(double &t, const valarray<double>& x, double dt, pt2Func eq)
{
		t+=dt;
		return eq(t, x)*dt;
}

void rkutta4(double &t, valarray<double> &x, double dt, pt2Func eq)
{
	
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim);
	k1=eq(t,x)*dt;
	t+=dt/2.0;
	k2=eq(t,x+k1*0.5)*dt;
	k3=eq(t,x+k2*0.5)*dt;
	t+=dt/2.0;
	k4=eq(t,x+k3)*dt;
	x += (k1 + k2*2.0 + k3*2.0 + k4)/6.0;
}

valarray<double> merson5(double &t, const valarray<double> &x, double &dt, double ode_eps, pt2Func eq)
{
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), delta;
	dt*=2.0;
	do
	{
		dt/=2.0;
		k1=eq(t,x)*dt/3.0;
		k2=eq(t+dt/3.0,x+k1)* dt/3.0;
		k3=eq(t+dt/3.0,x+k1/2.0+k2/2.0)* dt/3.0;
		k4=eq(t+dt/2.0,x+k1*0.375+k3*1.125)* dt/3.0;
		k5=eq(t+dt/2.0, x+k1*1.5 - k2*4.5 +k4*6.0) *dt/3.0;
		delta = k1 - k2*4.5 + k4*4.0 -k5*0.5;
	}
	while (pow(delta,2.0).sum()>5*ode_eps);
	t+=dt;
	if (pow(delta,2.0).sum() < 5*ode_eps/32) dt*=2;
	return (k1+k4*4.0+k5)/2.0;	
}

void dopri5(double &t, valarray<double> &x, double &h, double ode_eps, pt2Func eq)
{
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), k7(dim);
	const double two_thirds = 2.0 / 3.0;
  const double one_seventwoninths = 1.0 / 729.0;
   const double one_twoninesevenzero = 1.0 / 2970.0;
  const double one_twofivetwozero = 1.0 / 2520.0;
   const double one_fiveninefourzero = 1.0 / 5940.0;

      double h5 = 0.2 * h;

   k1 = eq(t, x);
   k2 = eq(t+h5, x + h5 * k1);
   k3 = eq(t+0.3*h, x + h * ( 0.075 * k1 + 0.225 * k2) );
   k4 = eq(t+0.6*h, x + h * ( 0.3 * k1 - 0.9 * k2 + 1.2 * k3) );
   k5 = eq(t+two_thirds * h,  x + one_seventwoninths * h * ( 226.0 * k1 
                                - 675.0 * k2 + 880.0 * k3 + 55.0 * k4) );
   k6 = eq(t+h, x + one_twoninesevenzero * h * ( - 1991.0 * k1 + 7425.0 * k2
                               - 2660.0 * k3 - 10010.0 * k4 + 10206.0 * k5) );
   x +=  one_fiveninefourzero * h * ( 341.0 * k1 + 3800.0 * k3
                                   - 7975.0 * k4 + 9477.0 * k5 + 297.0 * k6 );
   t+=h;
   //return one_twofivetwozero * (77.0 * k1 - 400.0 * k3 + 1925.0 * k4
     //                                          - 1701.0 * k5 + 99.0 * k6);

}*/

/*void butcher6(double &t, valarray<double> &x, double &dt, pt2Func eq)
{
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), k7(dim);
	k1=eq(t, x)*dt;
	k2=eq(t +dt/2.0, 	x +k1*(1.0/2.0))*dt;
	k3=eq(t +2.0*dt/3.0, x +k1*(2.0/9.0) 	 +k2*(4.0/9.0))*dt;
	k4=eq(t +dt/3.0, 	x +k1*(7.0/36.0)	 +k2*(2.0/9.0)   -k3*(1.0/12.0))*dt;
	k5=eq(t +5.0*dt/6.0, x -k1*(35.0/144.0) -k2*(55.0/36.0) +k3*(35.0/48.0)  +k4*(15.0/8.0))*dt;
	k6=eq(t +dt/6.0, 	x -k1*(1.0/360.0)  -k2*(11.0/36.0) -k3*(1.0/8.0)    +k4*(1.0/2.0)    +k5*(1.0/10.0))*dt;
	k7=eq(t +dt, 	x -k1*(41.0/260.0) +k2*(22.0/13.0) +k3*(43.0/156.0) -k4*(118.0/39.0) +k5*(32.0/195.0) +k6*(80.0/39.0))*dt;

	t+=dt;
	x += k1*(13.0/200.0)+k3*(11.0/40.0)+k4*(11.0/40.0)+k5*(4.0/25.0)+k6*(4.0/25.0)+k7*(13.0/200.0); 		
*/

void fehlberg8(double &t, valarray<double> &x, double &h, double hmin ,double tolerance, pt2Func eq)
{
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), 
		k7(dim), k8(dim), k9(dim), k10(dim), k11(dim), k12(dim), k13(dim), err(dim);
	 const double c_1_11 = 41.0 / 840.0;
    const double c6 = 34.0 / 105.0;
    const double c_7_8= 9.0 / 35.0;
    const double c_9_10 = 9.0 / 280.0;

    const double a2 = 2.0 / 27.0;
    const double a3 = 1.0 / 9.0;
    const double a4 = 1.0 / 6.0;
    const double a5 = 5.0 / 12.0;
    const double a6 = 1.0 / 2.0;
    const double a7 = 5.0 / 6.0;
    const double a8 = 1.0 / 6.0;
    const double a9 = 2.0 / 3.0;
    const double a10 = 1.0 / 3.0;

    const double b31 = 1.0 / 36.0;
    const double b32 = 3.0 / 36.0;
    const double b41 = 1.0 / 24.0;
    const double b43 = 3.0 / 24.0;
    const double b51 = 20.0 / 48.0;
    const double b53 = -75.0 / 48.0;
    const double b54 = 75.0 / 48.0;
    const double b61 = 1.0 / 20.0;
    const double b64 = 5.0 / 20.0;
    const double b65 = 4.0 / 20.0;
    const double b71 = -25.0 / 108.0;
    const double b74 =  125.0 / 108.0;
    const double b75 = -260.0 / 108.0;
    const double b76 =  250.0 / 108.0;
    const double b81 = 31.0/300.0;
    const double b85 = 61.0/225.0;
    const double b86 = -2.0/9.0;
    const double b87 = 13.0/900.0;
    const double b91 = 2.0;
    const double b94 = -53.0/6.0;
    const double b95 = 704.0 / 45.0;
    const double b96 = -107.0 / 9.0;
    const double b97 = 67.0 / 90.0;
    const double b98 = 3.0;
    const double b10_1 = -91.0 / 108.0;
    const double b10_4 = 23.0 / 108.0;
    const double b10_5 = -976.0 / 135.0;
    const double b10_6 = 311.0 / 54.0;
    const double b10_7 = -19.0 / 60.0;
    const double b10_8 = 17.0 / 6.0;
    const double b10_9 = -1.0 / 12.0;
    const double b11_1 = 2383.0 / 4100.0;
    const double b11_4 = -341.0 / 164.0;
    const double b11_5 = 4496.0 / 1025.0;
    const double b11_6 = -301.0 / 82.0;
    const double b11_7 = 2133.0 / 4100.0;
    const double b11_8 = 45.0 / 82.0;
    const double b11_9 = 45.0 / 164.0;
    const double b11_10 = 18.0 / 41.0;
    const double b12_1 = 3.0 / 205.0;
    const double b12_6 = - 6.0 / 41.0;
    const double b12_7 = - 3.0 / 205.0;
    const double b12_8 = - 3.0 / 41.0;
    const double b12_9 = 3.0 / 41.0;
    const double b12_10 = 6.0 / 41.0;
    const double b13_1 = -1777.0 / 4100.0;
    const double b13_4 = -341.0 / 164.0;
    const double b13_5 = 4496.0 / 1025.0;
    const double b13_6 = -289.0 / 82.0;
    const double b13_7 = 2193.0 / 4100.0;
    const double b13_8 = 51.0 / 82.0;
    const double b13_9 = 33.0 / 164.0;
    const double b13_10 = 12.0 / 41.0;
   
    const double err_factor  = -41.0 / 840.0;
    const double err_exponent = 1.0 / 7.0;	
    
    double module_err,yy, scale ;	

	h*=2.0;
	do
	{
	h/=2.0;
	double h2_7 = a2 * h;
   k1 = eq(t, x);
   k2 = eq(t+h2_7, x + h2_7 * k1);
   k3 = eq(t+a3*h, x + h * ( b31*k1 + b32*k2) );
   k4 = eq(t+a4*h, x + h * ( b41*k1 + b43*k3) );
   k5 = eq(t+a5*h, x + h * ( b51*k1 + b53*k3 + b54*k4) );
   k6 = eq(t+a6*h, x + h * ( b61*k1 + b64*k4 + b65*k5) );
   k7 = eq(t+a7*h, x + h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6) );
   k8 = eq(t+a8*h, x + h * ( b81*k1 + b85*k5 + b86*k6 + b87*k7) );
   k9 = eq(t+a9*h, x + h * ( b91*k1 + b94*k4 + b95*k5 + b96*k6
                                                          + b97*k7 + b98*k8) );
   k10 = eq(t+a10*h, x + h * ( b10_1*k1 + b10_4*k4 + b10_5*k5 + b10_6*k6
                                          + b10_7*k7 + b10_8*k8 + b10_9*k9 ) );
   k11 = eq(t+h, x + h * ( b11_1*k1 + b11_4*k4 + b11_5*k5 + b11_6*k6
                           + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10 ) );
   k12 = eq(t, x + h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8
                                                 + b12_9*k9 + b12_10 * k10 ) );
   k13 = eq(t+h, x + h * ( b13_1*k1 + b13_4*k4 + b13_5*k5 + b13_6*k6
                + b13_7*k7 + b13_8*k8 + b13_9*k9 + b13_10*k10 + k12 ) );

   err =  err_factor * (k1 + k11 - k12 - k13);
   //module_err = sqrt(err[0]*err[0]+err[1]*err[1]+err[2]*err[2]+err[3]*err[3]);   
module_err = sqrt(err[2]*err[2]+err[3]*err[3])/sqrt(x[2]*x[2]+x[3]*x[3]);   
   }
	while (module_err>5*tolerance);
	t+=h;
	x+=  h * (c_1_11 * (k1 + k11)  + c6 * k6 + c_7_8 * (k7 + k8) 
                                           + c_9_10 * (k9 + k10) );
	if (module_err < 5*tolerance/32) h*=2;       
}	

void fastFehlberg8(double &t, valarray<double> &x, double &h, double tolerance, pt2Func eq, valarray<double> *k)
{
	int dim=x.size();
	
	 const double c_1_11 = 41.0 / 840.0;
    const double c6 = 34.0 / 105.0;
    const double c_7_8= 9.0 / 35.0;
    const double c_9_10 = 9.0 / 280.0;

    const double a2 = 2.0 / 27.0;
    const double a3 = 1.0 / 9.0;
    const double a4 = 1.0 / 6.0;
    const double a5 = 5.0 / 12.0;
    const double a6 = 1.0 / 2.0;
    const double a7 = 5.0 / 6.0;
    const double a8 = 1.0 / 6.0;
    const double a9 = 2.0 / 3.0;
    const double a10 = 1.0 / 3.0;

    const double b31 = 1.0 / 36.0;
    const double b32 = 3.0 / 36.0;
    const double b41 = 1.0 / 24.0;
    const double b43 = 3.0 / 24.0;
    const double b51 = 20.0 / 48.0;
    const double b53 = -75.0 / 48.0;
    const double b54 = 75.0 / 48.0;
    const double b61 = 1.0 / 20.0;
    const double b64 = 5.0 / 20.0;
    const double b65 = 4.0 / 20.0;
    const double b71 = -25.0 / 108.0;
    const double b74 =  125.0 / 108.0;
    const double b75 = -260.0 / 108.0;
    const double b76 =  250.0 / 108.0;
    const double b81 = 31.0/300.0;
    const double b85 = 61.0/225.0;
    const double b86 = -2.0/9.0;
    const double b87 = 13.0/900.0;
    const double b91 = 2.0;
    const double b94 = -53.0/6.0;
    const double b95 = 704.0 / 45.0;
    const double b96 = -107.0 / 9.0;
    const double b97 = 67.0 / 90.0;
    const double b98 = 3.0;
    const double b10_1 = -91.0 / 108.0;
    const double b10_4 = 23.0 / 108.0;
    const double b10_5 = -976.0 / 135.0;
    const double b10_6 = 311.0 / 54.0;
    const double b10_7 = -19.0 / 60.0;
    const double b10_8 = 17.0 / 6.0;
    const double b10_9 = -1.0 / 12.0;
    const double b11_1 = 2383.0 / 4100.0;
    const double b11_4 = -341.0 / 164.0;
    const double b11_5 = 4496.0 / 1025.0;
    const double b11_6 = -301.0 / 82.0;
    const double b11_7 = 2133.0 / 4100.0;
    const double b11_8 = 45.0 / 82.0;
    const double b11_9 = 45.0 / 164.0;
    const double b11_10 = 18.0 / 41.0;
    const double b12_1 = 3.0 / 205.0;
    const double b12_6 = - 6.0 / 41.0;
    const double b12_7 = - 3.0 / 205.0;
    const double b12_8 = - 3.0 / 41.0;
    const double b12_9 = 3.0 / 41.0;
    const double b12_10 = 6.0 / 41.0;
    const double b13_1 = -1777.0 / 4100.0;
    const double b13_4 = -341.0 / 164.0;
    const double b13_5 = 4496.0 / 1025.0;
    const double b13_6 = -289.0 / 82.0;
    const double b13_7 = 2193.0 / 4100.0;
    const double b13_8 = 51.0 / 82.0;
    const double b13_9 = 33.0 / 164.0;
    const double b13_10 = 12.0 / 41.0;
   
    const double err_factor  = -41.0 / 840.0;
    const double err_exponent = 1.0 / 7.0;	
    
    double module_err,yy, scale ;	

	h*=2.0;
	do
	{
	h/=2.0;
	double h2_7 = a2 * h;
   k[1] = eq(t, x);
   k[2] = eq(t+h2_7, x + h2_7 *k[1]);
   k[3] = eq(t+a3*h, x + h * ( b31*k[1] + b32*k[2]) );
   k[4] = eq(t+a4*h, x + h * ( b41*k[1] + b43*k[3]) );
   k[5] = eq(t+a5*h, x + h * ( b51*k[1] + b53*k[3] + b54*k[4]) );
   k[6] = eq(t+a6*h, x + h * ( b61*k[1] + b64*k[4] + b65*k[5]) );
   k[7] = eq(t+a7*h, x + h * ( b71*k[1] + b74*k[4] + b75*k[5] + b76*k[6]) );
   k[8] = eq(t+a8*h, x + h * ( b81*k[1] + b85*k[5] + b86*k[6] + b87*k[7]) );
   k[9] = eq(t+a9*h, x + h * ( b91*k[1] + b94*k[4] + b95*k[5] + b96*k[6]
                                                         + b97*k[7] + b98*k[8]) );
   k[10] = eq(t+a10*h, x + h * ( b10_1*k[1] + b10_4*k[4] + b10_5*k[5] + b10_6*k[6]
                                          + b10_7*k[7] + b10_8*k[8] + b10_9*k[9] ) );
   k[11] = eq(t+h, x + h * ( b11_1*k[1] + b11_4*k[4] + b11_5*k[5] + b11_6*k[6]
                           + b11_7*k[7] + b11_8*k[8] + b11_9*k[9] + b11_10 *k[10] ) );
   k[12] = eq(t, x + h * ( b12_1*k[1] + b12_6*k[6] + b12_7*k[7] + b12_8*k[8]
                                                 + b12_9*k[9] + b12_10 *k[10] ) );
   k[13] = eq(t+h, x + h * ( b13_1*k[1] + b13_4*k[4] + b13_5*k[5] + b13_6*k[6]
                + b13_7*k[7] + b13_8*k[8] + b13_9*k[9] + b13_10*k[10] + k[12] ) );
   k[0] =  err_factor * (k[1] + k[11] - k[12] - k[13]);
   //module_err = sqrt([0]k[14][0]+k[14][1]k[14][1]+k[14][2]k[14][2]+k[14][3]k[14][3]);   
//module_err = sqrt(k[0][2]*k[0][2]+k[0][3]*k[0][3])/sqrt(x[2]*x[2]+x[3]*x[3]);   
   module_err = sqrt((k[0]*k[0]).sum())/sqrt((x*x).sum());	
   }
	while (module_err>5*tolerance);
	t+=h;
	x+=  h * (c_1_11 * (k[1] + k[11])  + c6 * k[6] + c_7_8 * (k[7] + k[8]) 
                                           + c_9_10 * (k[9] + k[10]) );
	if (module_err < 5*tolerance/32) h*=2;       
}

void modFastFehlberg8(double &t, valarray<double> &x, double &h, double hmin ,double tolerance, pt2Func eq, valarray<double> *k)
{
	int dim=x.size();
	
	 const double c_1_11 = 41.0 / 840.0;
    const double c6 = 34.0 / 105.0;
    const double c_7_8= 9.0 / 35.0;
    const double c_9_10 = 9.0 / 280.0;

    const double a2 = 2.0 / 27.0;
    const double a3 = 1.0 / 9.0;
    const double a4 = 1.0 / 6.0;
    const double a5 = 5.0 / 12.0;
    const double a6 = 1.0 / 2.0;
    const double a7 = 5.0 / 6.0;
    const double a8 = 1.0 / 6.0;
    const double a9 = 2.0 / 3.0;
    const double a10 = 1.0 / 3.0;

    const double b31 = 1.0 / 36.0;
    const double b32 = 3.0 / 36.0;
    const double b41 = 1.0 / 24.0;
    const double b43 = 3.0 / 24.0;
    const double b51 = 20.0 / 48.0;
    const double b53 = -75.0 / 48.0;
    const double b54 = 75.0 / 48.0;
    const double b61 = 1.0 / 20.0;
    const double b64 = 5.0 / 20.0;
    const double b65 = 4.0 / 20.0;
    const double b71 = -25.0 / 108.0;
    const double b74 =  125.0 / 108.0;
    const double b75 = -260.0 / 108.0;
    const double b76 =  250.0 / 108.0;
    const double b81 = 31.0/300.0;
    const double b85 = 61.0/225.0;
    const double b86 = -2.0/9.0;
    const double b87 = 13.0/900.0;
    const double b91 = 2.0;
    const double b94 = -53.0/6.0;
    const double b95 = 704.0 / 45.0;
    const double b96 = -107.0 / 9.0;
    const double b97 = 67.0 / 90.0;
    const double b98 = 3.0;
    const double b10_1 = -91.0 / 108.0;
    const double b10_4 = 23.0 / 108.0;
    const double b10_5 = -976.0 / 135.0;
    const double b10_6 = 311.0 / 54.0;
    const double b10_7 = -19.0 / 60.0;
    const double b10_8 = 17.0 / 6.0;
    const double b10_9 = -1.0 / 12.0;
    const double b11_1 = 2383.0 / 4100.0;
    const double b11_4 = -341.0 / 164.0;
    const double b11_5 = 4496.0 / 1025.0;
    const double b11_6 = -301.0 / 82.0;
    const double b11_7 = 2133.0 / 4100.0;
    const double b11_8 = 45.0 / 82.0;
    const double b11_9 = 45.0 / 164.0;
    const double b11_10 = 18.0 / 41.0;
    const double b12_1 = 3.0 / 205.0;
    const double b12_6 = - 6.0 / 41.0;
    const double b12_7 = - 3.0 / 205.0;
    const double b12_8 = - 3.0 / 41.0;
    const double b12_9 = 3.0 / 41.0;
    const double b12_10 = 6.0 / 41.0;
    const double b13_1 = -1777.0 / 4100.0;
    const double b13_4 = -341.0 / 164.0;
    const double b13_5 = 4496.0 / 1025.0;
    const double b13_6 = -289.0 / 82.0;
    const double b13_7 = 2193.0 / 4100.0;
    const double b13_8 = 51.0 / 82.0;
    const double b13_9 = 33.0 / 164.0;
    const double b13_10 = 12.0 / 41.0;
   
    const double err_factor  = -41.0 / 840.0;
    const double err_exponent = 1.0 / 7.0;	
    
    double module_err,yy, scale ;	

	
	do
	{
	
	double h2_7 = a2 * h;
   k[1] = eq(t, x);
   k[2] = eq(t+h2_7, x + h2_7 *k[1]);
   k[3] = eq(t+a3*h, x + h * ( b31*k[1] + b32*k[2]) );
   k[4] = eq(t+a4*h, x + h * ( b41*k[1] + b43*k[3]) );
   k[5] = eq(t+a5*h, x + h * ( b51*k[1] + b53*k[3] + b54*k[4]) );
   k[6] = eq(t+a6*h, x + h * ( b61*k[1] + b64*k[4] + b65*k[5]) );
   k[7] = eq(t+a7*h, x + h * ( b71*k[1] + b74*k[4] + b75*k[5] + b76*k[6]) );
   k[8] = eq(t+a8*h, x + h * ( b81*k[1] + b85*k[5] + b86*k[6] + b87*k[7]) );
   k[9] = eq(t+a9*h, x + h * ( b91*k[1] + b94*k[4] + b95*k[5] + b96*k[6]
                                                         + b97*k[7] + b98*k[8]) );
   k[10] = eq(t+a10*h, x + h * ( b10_1*k[1] + b10_4*k[4] + b10_5*k[5] + b10_6*k[6]
                                          + b10_7*k[7] + b10_8*k[8] + b10_9*k[9] ) );
   k[11] = eq(t+h, x + h * ( b11_1*k[1] + b11_4*k[4] + b11_5*k[5] + b11_6*k[6]
                           + b11_7*k[7] + b11_8*k[8] + b11_9*k[9] + b11_10 *k[10] ) );
   k[12] = eq(t, x + h * ( b12_1*k[1] + b12_6*k[6] + b12_7*k[7] + b12_8*k[8]
                                                 + b12_9*k[9] + b12_10 *k[10] ) );
   k[13] = eq(t+h, x + h * ( b13_1*k[1] + b13_4*k[4] + b13_5*k[5] + b13_6*k[6]
                + b13_7*k[7] + b13_8*k[8] + b13_9*k[9] + b13_10*k[10] + k[12] ) );
   k[0] =  err_factor * (k[1] + k[11] - k[12] - k[13]);
   //module_err = sqrt([0]k[14][0]+k[14][1]k[14][1]+k[14][2]k[14][2]+k[14][3]k[14][3]);   
module_err = sqrt((k[0]*k[0]).sum())/sqrt((x*x).sum());   
if (module_err>5*tolerance) h *=0.5*pow(tolerance/module_err,0.125);
   }
	while (module_err>5*tolerance);
	t+=h;
	x+=  h * (c_1_11 * (k[1] + k[11])  + c6 * k[6] + c_7_8 * (k[7] + k[8]) 
                                           + c_9_10 * (k[9] + k[10]) );
	h *=0.5*pow(tolerance/module_err,0.125);
	//if (module_err < 5*tolerance/32) h*=2;       
}	
	
void verner9(double &t, valarray<double> &x, double &h, double hmin, double tolerance, pt2Func eq)
{
	#define SQRT6 2.449489742783178098197284074705891

                  // c2 = c3 = ... = c7 = 0, c15 = c16 = 0 //
    int dim=x.size();
    const double c1 = 103.0 / 1680.0;
    const double c8 = -27.0 / 140.0;
    const double c9 = 76.0 / 105.0;
    const double c10 = -201.0 / 280.0;
    const double c11 = 1024.0 / 1365.0;
    const double c12 = 3.0 / 7280.0;
    const double c13 = 12.0 / 35.0;
    const double c14 = 9.0 / 280.0;

                          // a1 = 0, a14 = a16 = 1 //

    const double a2 = 1.0 / 12.0;
    const double a3 = 1.0 / 9.0;
    const double a4 = 1.0 / 6.0;
    const double a5 = 2.0 * (1.0 + SQRT6) / 15.0;
    const double a6 = (6.0 + SQRT6) / 15.0;
    const double a7 = (6.0 - SQRT6) / 15.0;
    const double a8 = 2.0 / 3.0;
    const double a9 = 1.0 / 2.0;
    const double a10 = 1.0 / 3.0;
    const double a11 = 1.0 / 4.0;
    const double a12 = 4.0 / 3.0;
    const double a13 = 5.0 / 6.0;
    const double a15 = 1.0 / 6.0;

             // b21 = 1/12, remaining missing bij's j < i are 0 //

    const double b31 = 1.0 / 27.0;
    const double b32 = 2.0 / 27.0;
    const double b41 = 1.0 / 24.0;
    const double b43 = 3.0 / 24.0;
    const double b51 = (4.0 + 94.0 * SQRT6) / 375.0;
    const double b53 = -(282.0 + 252.0 * SQRT6) / 375.0;
    const double b54 = (328.0 + 208.0 * SQRT6) / 375.0;
    const double b61 = (9.0 - SQRT6) / 150.0;
    const double b64 = (312.0 + 32.0 * SQRT6) / 1425.0;
    const double b65 = (69.0 + 29.0 * SQRT6) / 570.0;
    const double b71 = (927.0 - 347.0 * SQRT6) / 1250.0;
    const double b74 = (-16248.0 + 7328.0 * SQRT6) / 9375.0;
    const double b75 = (-489.0 + 179.0 * SQRT6) / 3750.0;
    const double b76 = (14268.0 - 5798.0 * SQRT6) / 9375.0;
    const double b81 = 4.0 / 54.0;
    const double b86 = (16.0 - SQRT6) / 54.0;
    const double b87 = (16.0 + SQRT6) / 54.0;
    const double b91 = 38.0 / 512.0;
    const double b96 = (118.0 - 23.0 * SQRT6) / 512.0;
    const double b97 = (118.0 + 23.0 * SQRT6) / 512.0;
    const double b98 = - 18.0 / 512.0;
    const double b10_1 = 11.0 / 144.0;
    const double b10_6 = (266.0 - SQRT6) / 864.0;
    const double b10_7 = (266.0 + SQRT6) / 864.0;
    const double b10_8 = - 1.0 / 16.0;
    const double b10_9 = - 8.0 / 27.0;
    const double b11_1 = (5034.0 - 271.0 * SQRT6) / 61440.0;
    const double b11_7 = (7859.0 - 1626.0 * SQRT6) / 10240.0;
    const double b11_8 = (-2232.0 + 813.0 * SQRT6) / 20480.0;
    const double b11_9 = (-594.0  + 271.0 * SQRT6) / 960.0;
    const double b11_10 = (657.0 - 813.0 * SQRT6) / 5120.0;
    const double b12_1 = (5996.0 - 3794.0 * SQRT6) / 405.0;
    const double b12_6 = (-4342.0 - 338.0 * SQRT6) / 9.0;
    const double b12_7 = (154922.0 - 40458.0 * SQRT6) / 135.0;
    const double b12_8 = (-4176.0 + 3794.0 * SQRT6) / 45.0;
    const double b12_9 = (-340864.0 + 242816.0 * SQRT6) / 405.0;
    const double b12_10 = (26304.0 - 15176.0 * SQRT6) / 45.0;
    const double b12_11 = -26624.0 / 81.0;
    const double b13_1 = (3793.0 + 2168.0 * SQRT6) / 103680.0;
    const double b13_6 = (4042.0 + 2263.0 * SQRT6) / 13824.0;
    const double b13_7 = (-231278.0 + 40717.0 * SQRT6) / 69120.0;
    const double b13_8 = (7947.0 - 2168.0 * SQRT6) / 11520.0;
    const double b13_9 = (1048.0 - 542.0 * SQRT6) / 405.0;
    const double b13_10 = (-1383.0 + 542.0 * SQRT6) / 720.0;
    const double b13_11 = 2624.0 / 1053.0;
    const double b13_12 = 3.0 / 1664.0;
    const double b14_1 = -137.0 / 1296.0;
    const double b14_6 = (5642.0 - 337.0 * SQRT6) / 864.0;
    const double b14_7 = (5642.0 + 337.0 * SQRT6) / 864.0;
    const double b14_8 = -299.0 / 48.0;
    const double b14_9 = 184.0 / 81.0;
    const double b14_10 = -44.0 / 9.0;
    const double b14_11 = -5120.0 / 1053.0;
    const double b14_12 = -11.0 / 468.0;
    const double b14_13 = 16.0 / 9.0;
    const double b15_1 = (33617.0 - 2168.0 * SQRT6) / 518400.0;
    const double b15_6 = (-3846.0 + 31.0 * SQRT6) / 13824.0;
    const double b15_7 = (155338.0 - 52807.0 * SQRT6) / 345600.0;
    const double b15_8 = (-12537.0 + 2168.0 * SQRT6) / 57600.0;
    const double b15_9 = (92.0 + 542.0 * SQRT6) / 2025.0;
    const double b15_10 = (-1797.0 - 542.0 * SQRT6) / 3600.0;
    const double b15_11 = 320.0 / 567.0;
    const double b15_12 = -1.0 / 1920.0;
    const double b15_13 = 4.0 / 105.0;
    const double b16_1 = (-36487.0 - 30352.0 * SQRT6) / 279600.0;
    const double b16_6 = (-29666.0 - 4499.0 * SQRT6) / 7456.0;
    const double b16_7 = (2779182.0 - 615973.0 * SQRT6) / 186400.0;
    const double b16_8 = (-94329.0 + 91056.0 * SQRT6) / 93200.0;
    const double b16_9 = (-232192.0 + 121408.0 * SQRT6) / 17475.0;
    const double b16_10 = (101226.0 - 22764.0 * SQRT6) / 5825.0;
    const double b16_11 = - 169984.0 / 9087.0;
    const double b16_12 = - 87.0 / 30290.0;
    const double b16_13 =  492.0 / 1165.0;
    const double b16_15 =  1260.0 / 233.0;

                           // e2 = 0, ..., e7 =0 //

    const double e1 = -1911.0 / 109200.0;
    const double e8 = 34398.0 / 109200.0;
    const double e9 = -61152.0 / 109200.0;
    const double e10 = 114660.0 / 109200.0;
    const double e11 = -114688.0 / 109200.0;
    const double e12 = -63.0 / 109200.0;
    const double e13 = -13104.0 / 109200.0;
    const double e14 = -3510.0 / 109200.0;
    const double e15 = 39312.0 / 109200.0;
    const double e16 = 6058.0 / 109200.0;

   valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), 
		k7(dim), k8(dim), k9(dim), k10(dim), k11(dim), k12(dim), k13(dim), k14(dim), k15(dim), k16(dim), err(dim);

double module_err,yy, scale ;	

	h*=2.0;
	do
	{
	h/=2.0;

   double h12 = a2 * h;
   double h6 = a3 * h;

   k1 = eq(t, x);
   k2 = eq(t+h12, x + h12 * k1);
   k3 = eq(t+a3*h, x + h * ( b31*k1 + b32*k2) );
   k4 = eq(t+a4*h, x + h * ( b41*k1 + b43*k3) );
   k5 = eq(t+a5*h, x + h * ( b51*k1 + b53*k3 + b54*k4) );
   k6 = eq(t+a6*h, x + h * ( b61*k1 + b64*k4 + b65*k5) );
   k7 = eq(t+a7*h, x + h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6) );
   k8 = eq(t+a8*h, x + h * ( b81*k1 + b86*k6 + b87*k7) );
   k9 = eq(t+a9*h, x + h * ( b91*k1 + b96*k6 + b97*k7 + b98*k8) );
   k10 = eq(t+a10*h, x + h * ( b10_1*k1 + b10_6*k6 + b10_7*k7 + b10_8*k8
                                                                + b10_9*k9 ) );
   k11 = eq(t+a11*h, x + h * ( b11_1*k1 + b11_7*k7 + b11_8*k8 + b11_9*k9
                                                            + b11_10 * k10 ) );
   k12 = eq(t+a12*h, x + h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8
                                  + b12_9*k9 + b12_10 * k10 + b12_11 * k11 ) );
   k13 = eq(t+a13*h, x + h * ( b13_1*k1 + b13_6*k6 + b13_7*k7 + b13_8*k8
                         + b13_9*k9 + b13_10*k10 + b13_11*k11 + b13_12*k12 ) );
   k14 = eq(t+h, x + h * ( b14_1*k1 + b14_6*k6 + b14_7*k7 + b14_8*k8
            + b14_9*k9 + b14_10*k10 + b14_11*k11 + b14_12*k12 + b14_13*k13 ) );
   k15 = eq(t+a15*h, x + h * ( b15_1*k1 + b15_6*k6 + b15_7*k7 + b15_8*k8
            + b15_9*k9 + b15_10*k10 + b15_11*k11 + b15_12*k12 + b15_13*k13 ) );
   k16 = eq(t+h, x + h * ( b16_1*k1 + b16_6*k6 + b16_7*k7 + b16_8*k8
            + b16_9*k9 + b16_10*k10 + b16_11*k11 + b16_12*k12 + b16_13*k13
                                                               + b16_15*k15) );
  err = e1*k1 + e8*k8 + e9*k9 + e10*k10 + e11*k11 + e12*k12 + e13*k13
                                                 + e14*k14 + e15*k15 + e16*k16;	
  module_err = sqrt((err*err).sum()/(x*x).sum());   
   }
	while (module_err>5*tolerance);
	t+=h;
	 x = x +  h * ( c1 * k1 + c8 * k8 + c9 * k9 + c10 * k10 + c11 * k11
                                         + c12 * k12 + c13 * k13 + c14 * k14 );
	if (module_err < 5*tolerance/32) h*=2;         
}

/*void everhart15(double &t, valarray<double> &x, double &dt, double tol, int ni, int ns, int nf, pt2FortranSub force, bool start)
{
	int nv = x.size();
	RADA_15_(&t, &x[0], &nv, &dt, &tol, &ni, &ns, &nf, force, &start);
} */

/*valarray<double> dopri8(double &t, valarray<double> x, double &dt, pt2Func eq)
{
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), 
		k7(dim), k8(dim), k9(dim), k10(dim), k11(dim), k12(dim), k13(dim);
	k1=eq(t, x)*dt;

	k2=eq(t + dt/18.0,x +k1*(1.0/18.0))*dt;

	k3=eq(t+dt/12.0,  x +k1*(1.0/48.0) +k2*(1.0/16.0))*dt;

	k4=eq(t+dt/8.0,   x +k1*(1.0/32.0) +k3*(3.0/32.0))*dt;

	k5=eq(t+ 5.0*dt/16.0, x + k1*(5.0/16.0) + k3*(-75.0/64.0) + k4*(75.0/64.0))*dt;

	k6=eq(t+ 3.0*dt/8.0, x +k1*(3.0/80.0) +k4*(3.0/16.0) +k5*(3.0/20.0))*dt;

	k7 = eq(t + 59.0*dt/400.0, x+ k1*(29443841.0/614563906.0) + k4*(77736538.0/692538347.0) 
		+ k5* (-28693883.0/1125000000.0) + k6*(23124283.0/1800000000.0))*dt;

	k8 = eq(t + 93.0*dt/200.0, x+ k1*(16016141.0/946692911.0) + k4*(61564180.0/158732637.0) 
		+ k5*(22789713.0/633445777.0) + k6*(545815736.0/2771057229.0) + k7*(-180193667.0/1043307555.0))*dt;

	k9 = eq(t +5490023248.0*dt/9719169821.0, x +k1*(39632708.0/573591083.0) +k4*(-433636366.0/683701615.0) +k5*(-421739975.0/2616292301.0) 
		+ k6*(100302831.0/723423059.0) + k7*(790204164.0/839813087.0) + k8*(800635310.0/3783071287.0))*dt;

	k10 = eq(t+13.0*dt/20.0, x + k1*(246121993.0/1340847787.0) + k4*(-37695042795.0/15268766246.0) + k5*(-309121744.0/1061227803.0)
		+k6*(-1299083.0/490766935.0) + k7*(6005943493.0/2108947869.0) + k8*(393006217.0/1396673457.0) 
		+ k9*(123872331.0/1001029789.0))*dt;

	k11 = eq(t+1201146811.0*dt/1299019798.0, x +  k1*(-1028468189.0/846180014.0) + k4*(8478235783.0/508512852.0) 
		+ k5*(1311729495.0/1432422823.0) +k6*(-10304129995.0/1701304382.0) + k7*(-48777925059.0/304793956.0) 
		+ k8*(15336726248.0/1032824649.0)+ k9*(-45442868181.0/3398467696.0) + k10*(3065993473.0/597172653.0))*dt;

	k12 = eq(t+dt, x + k1*(185892177.0/718116043.0)+ k4*(-3185094517.0/667107341.0) + k5*(-477755414.0/1098053517.0) 
			 +k6*(-703635378.0/230739211.0)	+ k7*(5731566787.0/1027545527.0) + k8*(5232866602.0/850066563.0) 
			+ k9*(-4093664535.0/808688257.0) +k10*(3962137247.0/1805957418.0)+ k11*(65686358.0/487910083.0) ) *dt;

	k13 = eq(t+dt, x+ k1*(403863854.0/491063109.0) + k4*(-5068492393.0/434740067.0) + k5*(-411421997.0/543043805.0) 
			+ k6*(652783627.0/914296604.0)	+ k7*(11173962825.0/925320556.0) + k8*(-13158990841.0/6184727034.0) 
			+ k9*(3936647629.0/1978049680.0) +k10*(-160528059.0/685178525.0) + k11*(248638103.0/1413531060.0)) *dt;

	//x = x + k1*(14005451.0/335480064.0) + k6*(-59238493.0/1068277825.0) + k7*(181606767.0/758867731.0) + k8*(561292985.0/797845732.0)
	//	+k9*(-1041891430.0/1371343529.0) + k10*(760417239.0/1151165299.0) + k11*(118820643.0/751138087.0) 
	//	+ k12*(-528747749.0/2220607170.0)+ k13*(1.0/4.0);
	x=x+ k1*(13451932.0/455176623.0) + k6*(-808719846.0/976000145.0) + k7*(1757004468.0/5645159321.0) + k8*(656045339.0/265891186.0)
		+ k9*(-3867574721.0/1518517206.0) + k10*(465885868.0/322736535.0) + k11*(53011238.0/667516719.0) + k12*(2.0/45.0);
	t+=dt;
	return x;	  
}*/

/*valarray<double> dopri8(double &t, const valarray<double> &x, double &dt, pt2Func eq)
{
	int dim=x.size();
	valarray<double> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), 
		k7(dim), k8(dim), k9(dim), k10(dim), k11(dim), k12(dim), k13(dim);
	k1=eq(t + dt, x)*dt;

	k2=eq(t + dt/18.0,x +k1*(1.0/18.0))*dt;

	k3=eq(t+dt/12.0,  x +k1*(1.0/48.0) +k2*(1.0/16.0))*dt;

	k4=eq(t+dt/8.0,   x +k1*(1.0/32.0) +k3*(3.0/32.0))*dt;

	k5=eq(t+ 5.0*dt/16.0, x + k1*(5.0/16.0) + k3*(-75.0/64.0) + k4*(75.0/64.0))*dt;

	k6=eq(t+ 3.0*dt/8.0, x +k1*(3.0/80.0) +k4*(3.0/16.0) +k5*(3.0/20.0))*dt;

	k7 = eq(t + 59.0*dt/400.0, x+ k1*(29443841.0/614563906.0) + k4*(77736538.0/692538347.0) 
		+ k5* (-28693883.0/1125000000.0) + k6*(23124283.0/1800000000.0))*dt;

	k8 = eq(t + 93.0*dt/200.0, x+ k1*(16016141.0/946692911.0) + k4*(61564180.0/158732637.0) 
		+ k5*(22789713.0/633445777.0) + k6*(545815736.0/2771057229.0) + k7*(-180193667.0/1043307555.0))*dt;

	k9 = eq(t +5490023248.0*dt/9719169821.0, x +k1*(39632708.0/573591083.0) +k4*(-433636366.0/683701615.0) +k5*(-421739975.0/2616292301.0) 
		+ k6*(100302831.0/723423059.0) + k7*(790204164.0/839813087.0) + k8*(800635310.0/3783071287.0))*dt;

	k10 = eq(t+13.0*dt/20.0, x + k1*(246121993.0/1340847787.0) + k4*(-37695042795.0/15268766246.0) + k5*(-309121744.0/1061227803.0)
		+k6*(-1299083.0/490766935.0) + k7*(6005943493.0/2108947869.0) + k8*(393006217.0/1396673457.0) 
		+ k9*(123872331.0/1001029789.0))*dt;

	k11 = eq(t+1201146811.0*dt/1299019798.0, x +  k1*(-1028468189.0/846180014.0) + k4*(8478235783.0/508512852.0) 
		+ k5*(1311729495.0/1432422823.0) +k6*(-10304129995.0/1701304382.0) + k7*(-48777925059.0/304793956.0) 
		+ k8*(15336726248.0/1032824649.0)+ k9*(-45442868181.0/3398467696.0) + k10*(3065993473.0/597172653.0))*dt;

	k12 = eq(t+dt, x + k1*(185892177.0/718116043.0)+ k4*(-3185094517.0/667107341.0) + k5*(-477755414.0/1098053517.0) 
			 +k6*(-703635378.0/230739211.0)	+ k7*(5731566787.0/1027545527.0) + k8*(5232866602.0/850066563.0) 
			+ k9*(-4093664535.0/808688257.0) +k10*(3962137247.0/1805957418.0)+ k11*(65686358.0/487910083.0) ) *dt;

	k13 = eq(t+dt, x+ k1*(403863854.0/491063109.0) + k4*(-5068492393.0/434740067.0) + k5*(-411421997.0/543043805.0) 
			+ k6*(652783627.0/914296604.0)	+ k7*(11173962825.0/925320556.0) + k8*(-13158990841.0/6184727034.0) 
			+ k9*(3936647629.0/1978049680.0) +k10*(-160528059.0/685178525.0) + k11*(248638103.0/1413531060.0)) *dt;

	//x = x + k1*(14005451.0/335480064.0) + k6*(-59238493.0/1068277825.0) + k7*(181606767.0/758867731.0) + k8*(561292985.0/797845732.0)
	//	+k9*(-1041891430.0/1371343529.0) + k10*(760417239.0/1151165299.0) + k11*(118820643.0/751138087.0) 
	//	+ k12*(-528747749.0/2220607170.0)+ k13*(1.0/4.0);
	t+=dt;
	return k1*(13451932.0/455176623.0) + k6*(-808719846.0/976000145.0) + k7*(1757004468.0/5645159321.0) + k8*(656045339.0/265891186.0)
		+ k9*(-3867574721.0/1518517206.0) + k10*(465885868.0/322736535.0) + k11*(53011238.0/667516719.0) + k12*(2.0/45.0);		  
}*/
