
#include "restricted_problem.h"

static int dim=14;

static double mu = 0.5;
static double a_sat=0.001;	
static double n=1.0;
valarray<double>rb0({0.0, 0.0, -mu, 0.0});
valarray<double>rb1({0.0, 0.0, 1-mu, 0.0});
static double globalFLI=0.0, globalOFLI=0.0;
bool orthogonal;
double norm2(const valarray<double> &x)
{
	return x[2]*x[2]+x[3]*x[3];
}

double norm(const valarray<double> &x)
{
	return sqrt(x[2]*x[2]+x[3]*x[3]);
}

double norm3(const valarray<double> &x)
{
	double result = x[2]*x[2]+x[3]*x[3];
	return result*sqrt(result);
}

double norm5(const valarray<double> &x)
{
	double result = x[2]*x[2]+x[3]*x[3];
	return result*result*sqrt(result);
}



double module(const valarray<double> &r)
{
	double sum=0;
	for (int i=0; i<r.size(); i++) sum+=r[i]*r[i];
	return sqrt(sum);
}

double square_module(const valarray<double> &r)
{
	double sum=0;
	for (int i=0; i<r.size(); i++) sum+=r[i]*r[i];
	return sum;
}

double scalarProduct(const valarray<double> &v1 , const valarray<double> &v2)
{
	if (v1.size() != v2.size()) return 0;
	double sum=0;
	for (int i=0; i<v1.size(); i++) sum+=v1[i]*v2[i];
	return sum; 
}


valarray<double> onlyForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];
	return f;
}

void allIndicator(double &FLI, double &MEGNO, double &average_MEGNO, double &SALI, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	FLI+=log(dev);
	globalFLI=FLI;
	//FLI=log(dev);
	dev = sqrt(x[8]*x[8]+x[9]*x[9]+x[10]*x[10]+x[11]*x[11]);
	for(int i=8; i<12; i++) x[i]/=dev;			
		
	double d1 = sqrt((x[8]-x[4])*(x[8]-x[4]) + (x[9]-x[5])*(x[9]-x[5]) + (x[10]-x[6])*(x[10]-x[6]) + (x[11]-x[7])*(x[11]-x[7])   );
	double d2 = sqrt((x[8]+x[4])*(x[8]+x[4]) + (x[9]+x[5])*(x[9]+x[5]) + (x[10]+x[6])*(x[10]+x[6]) + (x[11]+x[7])*(x[11]+x[7])   );
	SALI = (d1<=d2) ? d1 : d2;
	//SALI=d1;
	MEGNO=2*log(FLI - x[12]/t);
	average_MEGNO = 2*(x[12]-x[13])/t;	
}

valarray<double> fliForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];

	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);

	double D14 = 2*n;
	double D23 = -2*n;

	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6]; 
	return f;
}

void liIndicatorNoNorm(double &LI, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	LI= (t == 0) ? 0 : log(dev)/t;
}

void fliIndicator(double &FLI, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	FLI+=log(dev);		
}

void fliIndicatorNoNorm(double &FLI, double t, const valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	FLI=log(dev);
}

void ofliIndicator(double &OFLI, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	
	valarray<double> delta( {x[4], x[5], x[6], x[7]} );
	valarray<double> ff( fliForce(t,x));
	valarray<double> f({ff[2],ff[3],ff[0],ff[1]});
	double f_module =f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3];
	valarray<double> d0 = delta - f * (scalarProduct(delta,f)/f_module);
	OFLI+=log(1 + ((dev-1)/dev)*sqrt((d0*d0).sum()));
	for(int i=4; i<8; i++) x[i]/=dev;	
}	

/*void predOfliIndicator(double &OFLI, double t, valarray<double> &x, double &pred_dev)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	
	valarray<double> delta( {x[4], x[5], x[6], x[7]} );
	valarray<double> ff( fliForce(t,x));
	valarray<double> f({ff[2],ff[3],ff[0],ff[1]});
	double f_module =f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3];
	valarray<double> d0 = delta - f * (scalarProduct(delta,f)/f_module);
	OFLI+=log( pred_dev * sqrt((d0*d0).sum()) );
	for(int i=4; i<8; i++) x[i]/=dev;	
	pred_dev*=dev;
}*/
valarray <double> project(const valarray<double> &r, const valarray<double> f)
{
	double f_module =f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3];
	return r - f * (scalarProduct(r,f)/f_module);
}

void ofliIndicatorNoNorm(double &OFLI, double t, valarray<double> &x)
{
	valarray<double> delta( {x[4], x[5], x[6], x[7]} );
	valarray<double> ff( fliForce(t,x));
	valarray<double> f({ff[2],ff[3],ff[0],ff[1]});
	OFLI=log10(sqrt((project(delta,f)*project(delta,f)).sum()));
}	

valarray<double> megnoForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);
	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];	
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6];
	double delta = sqrt(r[4]*r[4]+r[5]*r[5]+r[6]*r[6]+r[7]*r[7]);
	f[8]=  globalFLI; //teta
	f[9]= (t==0) ? 0 : r[8]/t;
	return f;
}

void megnoIndicator(double &MEGNO, double &average_MEGNO, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	globalFLI+=log(dev);
	MEGNO=2*log(globalFLI - x[8]/t);
	average_MEGNO = 2*(x[8]-x[9])/t;	
}
/*valarray<double> omegnoForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);
	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
	f[1]= - 2 * r[0]*n + r[3] - mu*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];	
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6];
	
	f[8] = globalFLI; //teta
	f[9]= (t==0) ? 0 : r[8]/t;
	return f;
}*/
void omegnoIndicatorNoNormalize(double &MEGNO, double &average_MEGNO, double t, valarray<double> &x)
{
	valarray<double> delta( {x[4], x[5], x[6], x[7]} );
	valarray<double> ff( onlyForce(t,x));
	valarray<double> f({ff[2],ff[3],ff[0],ff[1]});
	globalFLI=log(sqrt((project(delta,f)*project(delta,f)).sum()));
	//globalFLI=log(dev);
	MEGNO=2*log(globalFLI - x[8]/t);
	average_MEGNO = 2*(x[8]-x[9])/t;
}

void megnoIndicatorNoNormalize(double &MEGNO, double &average_MEGNO, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	globalFLI=log(dev);
	MEGNO=2*log(globalFLI - x[8]/t);
	average_MEGNO = 2*(x[8]-x[9])/t;	
}

valarray<double> modMegnoForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);
	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];	
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6];
	valarray<double> delta({r[4], r[5], r[6], r[7]});
	valarray<double> dotdelta({f[4], f[5], f[6], f[7]});
	f[8] = r[9]; //teta
	f[9] = scalarProduct(dotdelta, delta) / (delta*delta).sum();
	f[10] = (t==0) ? 0: r[8]/t;
	return f;
}

void modMegnoIndicator(double &MEGNO1, double &MEGNO2, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	MEGNO1=2*(x[9] - x[8]/t);
	MEGNO2 = 2*(x[8]-x[10])/t;	
}

valarray<double> classicMegnoForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);
	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];	
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6];
	valarray<double> delta({r[4], r[5], r[6], r[7]});
	valarray<double> dotdelta({f[4], f[5], f[6], f[7]});
	//f[8] = r[9]; //teta
	f[8] = t*scalarProduct(dotdelta, delta) / (delta*delta).sum();
	f[9] = (t==0) ? 0: 2*r[8]/t;
	return f;
}

void classicMegnoIndicator(double &MEGNO1, double &MEGNO2, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	MEGNO1=2*x[8]/t;
	MEGNO2 = x[9]/t;	
}

valarray<double> megno_omegnoForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);
	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];	
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6];
	f[8]=  globalFLI; //teta
	f[9]= (t==0) ? 0 : r[8]/t;
	f[10]=  globalOFLI; //teta
	f[11]= (t==0) ? 0 : r[10]/t;
	return f;
}

void megno_omegnoIndicator(double &MEGNO_OMEGNO1, double &MEGNO_OMEGNO2, double t, valarray<double> &x)
{
	double MEGNO1=0.0, MEGNO2=0.0, OMEGNO1=0.0, OMEGNO2=0.0;

	valarray<double> delta( {x[4], x[5], x[6], x[7]} );
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	globalFLI=log(dev);
	MEGNO1=2*log(globalFLI - x[8]/t);
	MEGNO2 = 2*(x[8]-x[9])/t;

	valarray<double> ff( onlyForce(t,x));
	valarray<double> f({ff[2],ff[3],ff[0],ff[1]});
	globalOFLI=log(sqrt((project(delta,f)*project(delta,f)).sum()));
	
	OMEGNO1 = fabs(2*log(globalOFLI - x[10]/t));
	OMEGNO2 = fabs(2*(x[10]-x[11])/t);

	MEGNO_OMEGNO2 = MEGNO2 - MEGNO2*(exp(4.61*(MEGNO2-OMEGNO2))-1)/(exp(4.61*MEGNO2)-1);
	MEGNO_OMEGNO1 = MEGNO1 - MEGNO1*(exp(4.61*(MEGNO1-OMEGNO1))-1)/(exp(4.61*MEGNO1)-1);
}


valarray<double> saliForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);

	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6]; 
	f[8] = r[10];
	f[9] = r[11];
	f[10] = D11*r[8]+D12*r[9]+D14*r[11];
	f[11] = D21*r[8]+D22*r[9]+D23*r[10]; 
	
	return f;
}

void saliIndicator(double &SALI, double t, valarray<double> &r)
{
	double dev1 = sqrt(r[4]*r[4]+r[5]*r[5]+r[6]*r[6]+r[7]*r[7]);
	double dev2 = sqrt(r[8]*r[8]+r[9]*r[9]+r[10]*r[10]+r[11]*r[11]);
	
	double d1 = sqrt((r[8]/dev2-r[4]/dev1)*(r[8]/dev2-r[4]/dev1) + (r[9]/dev2-r[5]/dev1)*(r[9]/dev2-r[5]/dev1) + (r[10]/dev2-r[6]/dev1)*(r[10]/dev2-r[6]/dev1) + (r[11]/dev2-r[7]/dev1)*(r[11]/dev2-r[7]/dev1)   );
	double d2 = sqrt((r[8]/dev2+r[4]/dev1)*(r[8]/dev2+r[4]/dev1) + (r[9]/dev2+r[5]/dev1)*(r[9]/dev2+r[5]/dev1) + (r[10]/dev2+r[6]/dev1)*(r[10]/dev2+r[6]/dev1) + (r[11]/dev2+r[7]/dev1)*(r[11]/dev2+r[7]/dev1)   );
	//f[12] = (d1<=d2) ? d1 : d2;
	SALI = (d1<=d2) ? d1 : d2;
}

/*void gali2Indicator(double &GALI2, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	dev = sqrt(x[8]*x[8]+x[9]*x[9]+x[10]*x[10]+x[11]*x[11]);
	for(int i=8; i<12; i++) x[i]/=dev;	
	double d1 = sqrt((x[8]-x[4])*(x[8]-x[4]) + (x[9]-x[5])*(x[9]-x[5]) + (x[10]-x[6])*(x[10]-x[6]) + (x[11]-x[7])*(x[11]-x[7])   );
	double d2 = sqrt((x[8]+x[4])*(x[8]+x[4]) + (x[9]+x[5])*(x[9]+x[5]) + (x[10]+x[6])*(x[10]+x[6]) + (x[11]+x[7])*(x[11]+x[7])   );

	GALI2 = d1*d2/2;
}*/

valarray<double> modSaliForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);

	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6]; 
	f[8] = r[10];
	f[9] = r[11];
	f[10] = D11*r[8]+D12*r[9]+D14*r[11];
	f[11] = D21*r[8]+D22*r[9]+D23*r[10]; 

	double dev1 = sqrt(r[4]*r[4]+r[5]*r[5]+r[6]*r[6]+r[7]*r[7]);
	double dev2 = sqrt(r[8]*r[8]+r[9]*r[9]+r[10]*r[10]+r[11]*r[11]);
	
	double d1 = sqrt((r[8]/dev2-r[4]/dev1)*(r[8]/dev2-r[4]/dev1) + (r[9]/dev2-r[5]/dev1)*(r[9]/dev2-r[5]/dev1) + (r[10]/dev2-r[6]/dev1)*(r[10]/dev2-r[6]/dev1) + (r[11]/dev2-r[7]/dev1)*(r[11]/dev2-r[7]/dev1)   );
	double d2 = sqrt((r[8]/dev2+r[4]/dev1)*(r[8]/dev2+r[4]/dev1) + (r[9]/dev2+r[5]/dev1)*(r[9]/dev2+r[5]/dev1) + (r[10]/dev2+r[6]/dev1)*(r[10]/dev2+r[6]/dev1) + (r[11]/dev2+r[7]/dev1)*(r[11]/dev2+r[7]/dev1)   );
	f[12] = (d1<=d2) ? d1 : d2;

	return f;
}

void modSaliIndicator(double &SALI, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	dev = sqrt(x[8]*x[8]+x[9]*x[9]+x[10]*x[10]+x[11]*x[11]);
	for(int i=8; i<12; i++) x[i]/=dev;			
		
	SALI = (t!=0) ? x[12]/t : 0;
}

void modSaliIndicatorNoNorm(double &SALI, double t, valarray<double> &x)
{
	SALI = (t!=0) ? x[12]/t : 0;
}
valarray<double> gali2Force(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);

	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6]; 
	f[8] = r[10];
	f[9] = r[11];
	f[10] = D11*r[8]+D12*r[9]+D14*r[11];
	f[11] = D21*r[8]+D22*r[9]+D23*r[10]; 
	double dev1 = sqrt(r[4]*r[4]+r[5]*r[5]+r[6]*r[6]+r[7]*r[7]);
	double dev2 = sqrt(r[8]*r[8]+r[9]*r[9]+r[10]*r[10]+r[11]*r[11]);
	
	double d1 = sqrt((r[8]/dev2-r[4]/dev1)*(r[8]/dev2-r[4]/dev1) + (r[9]/dev2-r[5]/dev1)*(r[9]/dev2-r[5]/dev1) + (r[10]/dev2-r[6]/dev1)*(r[10]/dev2-r[6]/dev1) + (r[11]/dev2-r[7]/dev1)*(r[11]/dev2-r[7]/dev1)   );
	double d2 = sqrt((r[8]/dev2+r[4]/dev1)*(r[8]/dev2+r[4]/dev1) + (r[9]/dev2+r[5]/dev1)*(r[9]/dev2+r[5]/dev1) + (r[10]/dev2+r[6]/dev1)*(r[10]/dev2 + r[6]/dev1) + (r[11]/dev2+r[7]/dev1)*(r[11]/dev2+r[7]/dev1)   );
	f[12] = d1*d2/2;

	return f;
}

void gali2Indicator(double &GALI2, double t, valarray<double> &x)
{
	/*double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	dev = sqrt(x[8]*x[8]+x[9]*x[9]+x[10]*x[10]+x[11]*x[11]);
	for(int i=8; i<12; i++) x[i]/=dev;
*/			
		
	GALI2 = (t!=0) ? x[12]/t : 0;
}

valarray<double> classicGali2Force(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);

	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6]; 
	f[8] = r[10];
	f[9] = r[11];
	f[10] = D11*r[8]+D12*r[9]+D14*r[11];
	f[11] = D21*r[8]+D22*r[9]+D23*r[10]; 
	return f;
}

void classicGali2Indicator(double &GALI2, double t, valarray<double> &r)
{
	double dev1 = sqrt(r[4]*r[4]+r[5]*r[5]+r[6]*r[6]+r[7]*r[7]);
	double dev2 = sqrt(r[8]*r[8]+r[9]*r[9]+r[10]*r[10]+r[11]*r[11]);
	
	double d1 = sqrt((r[8]/dev2-r[4]/dev1)*(r[8]/dev2-r[4]/dev1) + (r[9]/dev2-r[5]/dev1)*(r[9]/dev2-r[5]/dev1) + (r[10]/dev2-r[6]/dev1)*(r[10]/dev2-r[6]/dev1) + (r[11]/dev2-r[7]/dev1)*(r[11]/dev2-r[7]/dev1)   );
	double d2 = sqrt((r[8]/dev2+r[4]/dev1)*(r[8]/dev2+r[4]/dev1) + (r[9]/dev2+r[5]/dev1)*(r[9]/dev2+r[5]/dev1) + (r[10]/dev2+r[6]/dev1)*(r[10]/dev2 + r[6]/dev1) + (r[11]/dev2+r[7]/dev1)*(r[11]/dev2+r[7]/dev1)   );
	GALI2 = d1*d2/2;	
}
valarray<double> osaliForce(double t, const valarray<double> &r)
{
	valarray<double> f(dim);
	//valarray<double>r0({0.0, 0.0, 1, 0.0});
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);		
	double D11= 1 - (1-mu)*(norm2(r0) - 3*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - 3 * (r[2]-1+mu)*(r[2]-1+mu)) /norm5(r1);
	double D12 = (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r0[3]*(r[2]-1+mu)/norm5(r1);
	double D21= (1-mu)*3*r[3]*(r[2]+mu)/norm5(r0) + mu*3*r[3]*(r[2]-1+mu)/norm5(r1);
	double D22 = 1 - (1-mu)*(norm2(r0) - 3*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - 3*r[3]*r[3])/norm5(r1);

	double D14 = 2*n;
	double D23 = -2*n;

	f[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
	f[2] = r[0];
	f[3] = r[1];
	f[4] = r[6];
	f[5] = r[7];
	f[6] = D11*r[4]+D12*r[5]+D14*r[7];
	f[7] = D21*r[4]+D22*r[5]+D23*r[6]; 
	f[8] = r[10];
	f[9] = r[11];
	f[10] = D11*r[8]+D12*r[9]+D14*r[11];
	f[11] = D21*r[8]+D22*r[9]+D23*r[10]; 
	valarray<double>ff({r[2],r[3],r[0],r[1]});
	valarray<double>d1({r[4],r[5],r[6],r[7]});
	d1=project(d1,ff);
	d1/=sqrt((d1*d1).sum());
	valarray<double>d2({r[8], r[9], r[10], r[11]});
	d2=project(d2,ff);
	d2/=sqrt((d2*d2).sum());
	double m1 = sqrt(((d1-d2)*(d1-d2)).sum());
	double m2 = sqrt(((d1+d2)*(d1+d2)).sum());
	f[12] = (m1<=m2) ? m1 : m2;

	return f;
}

void osaliIndicator(double &SALI, double t, valarray<double> &x)
{
	double dev = sqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
	for(int i=4; i<8; i++) x[i]/=dev;
	dev = sqrt(x[8]*x[8]+x[9]*x[9]+x[10]*x[10]+x[11]*x[11]);
	for(int i=8; i<12; i++) x[i]/=dev;			
		
	SALI = (t!=0) ? x[12]/t : 0;
}
double jacobiIntegral(const valarray<double> &r)
{
	return r[2]*r[2]+r[3]*r[3] + 2*(1-mu)/norm(r-rb0) + 2*mu/norm(r-rb1) - r[0]*r[0] -r[1]*r[1];
}

double secondIntegral(const valarray<double> &r)
{
	valarray<double>r0(r-rb0);
	valarray<double>r1(r-rb1);
	return 2*((r[2] - (1-mu)*(r[2]+mu)/norm3(r0) - mu*(r[2]-1+mu)/norm3(r1) )*r[4] + (r[3] - (1-mu)*r[3]/norm3(r0) - mu*r[3]/norm3(r1))*r[5] - r[0]*r[6] - r[1]*r[7]);
}

valarray<double> jacobiGradient(const valarray<double> &r)
{
	valarray<double> f(dim);
	f[0] = -2*r[0];
	f[1] = -2*r[1];
	f[2]= 2*r[2] - 2* (1-mu)* (r[2]-rb0[2])/norm3(r - rb1) - 2*mu*(r[2]-rb1[2])/norm3(r-rb0);
	f[3] = 2*r[3] - 2* (1-mu) * (r[3]-rb0[3])/norm3(r - rb1) - 2*mu*r[3]/norm3(r-rb0);
	return f;
}


void inverseProblem(const valarray<double> &x, double &a, double &e, double &u)
{
	valarray<double> r = x[ slice(2,2,1) ];
	valarray<double> v = x[ slice(0,2,1) ];
	valarray<double> I({1,0});
	valarray<double> e_vector(2);
	
	e_vector = r*(pow(module(v),2)/mu - 1/module(r)) - v*scalarProduct(r,v)/mu;
	double c = (r[0]*v[1]-r[1]*v[0]);
	double p = c*c/mu;
	e = module(e_vector);
	a = p/(1-e*e);
	u = acos(scalarProduct(e_vector, I)/e);	
}	

void findmax(vector<double> &sol)
{
	float min = numeric_limits<double>::max();
 	float max = - numeric_limits<double>::max();
		for(int j=0; j<sol.size(); j++)
		{
			float c=sol[j];
			if (c<min) min=c;
			if (c>max) max=c;
		} 
	sol.push_back(min);
	sol.push_back(max);
}


void integrate(double C0, double x0, double dt, double tmax, vector<double>& xaxis, vector<double>& yaxis, int whichIndicator, double tolerance, double &ind, double &CPU_time, bool evolution)
{
	clock_t t0=clock();
    valarray<double> *initial;
	ind=0;
	if (whichIndicator==1) {
		dim=8;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1, 0, 0, 0});
		//cout<< "FLI = ";
	}
	else if (whichIndicator==2) {
		dim=10;
		//cout<<"MEGNO = ";
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) ,x0, 0,   1,0,0,0, 0,0});
		
	}
	else if (whichIndicator==3) {
		dim=8;
		//<<"Modified MEGNO = ";
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) ,x0, 0,   1,0,0,0});
	}
	else if (whichIndicator==4) {
		dim=13;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) ,x0, 0,   1,0,0,0, 0,0,0,1, 0});
		//cout<<"SALI = ";
	}
	else if (whichIndicator==5) {
		dim=13;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0,  0,0,0,1, 0});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==6) {
		dim=12;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0, 0,0,0,1});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==7) {
		dim=8;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==8) {
		dim=10;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0, 0,0});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==9) {
		dim=13;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0, 0,1,0,0, 0});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==10) {
		dim=12;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0, 0,0,0,1});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==11) {
		dim=12;
		//initial.resize(dim);
		initial = new valarray<double>({0, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0)) - C0) , x0, 0, 1,0,0,0, 0,0, 0,0});
		//cout<< "OFLI = ";
	}
	else if (whichIndicator==0) {
		dim=4;
		//initial.resize(dim);
		double vx=0.0;
		initial = new valarray<double>({vx, sqrt(x0*x0 + 2*(1-mu)/(abs(x0+mu)) +2*mu/(abs(-1+mu+x0) - vx*vx) - C0) ,x0, 0});		
	}

	else {
		cout<<"Wrong parameter!";
		return;
	}

	valarray<double> r(dim), pred(dim);
	r = *initial;
	double initJacobiIntegral = jacobiIntegral(r);
	double initSecondIntegral = secondIntegral(r);
    cout<< "x0 = "  << x0 <<endl;
	//for (int p=0; p<dim; p++) cout << endl << initial[p]<< " ";
		
    double t=0, FLI = 0, OFLI=0, LI =0, max_FLI = 0, MEGNO1=0, MEGNO2=0, OMEGNO1=0, OMEGNO2=0, UMEGNO1=0, UMEGNO2=0, SALI=0, GALI2, min=0.01, error1=0, error2 = 0, t_pred=0, pred_ofli=1;
	int i=0;

	
	valarray<double> k[14] = valarray<double>(dim);
	//cout<<"x0 = " << x0 <<endl;
	//{
		for (i; t<tmax;  i++)
		{
			if (whichIndicator==1) {
				fastFehlberg8(t,r,dt, tolerance, &fliForce,&k[0]);
				//rkutta4(t,r,dt,&fliForce);
				fliIndicator(FLI, t,r);
				if (FLI>max_FLI) max_FLI=FLI;
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(max_FLI);}
				ind=max_FLI;
			}
			else if (whichIndicator==2) {
				fastFehlberg8(t,r,dt,  tolerance, &megnoForce,&k[0]);
				megnoIndicatorNoNormalize(MEGNO1, MEGNO2, t,r);
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(MEGNO2);}
				ind=MEGNO2;
			}
			else if (whichIndicator==3) {
				fastFehlberg8(t,r,dt,  tolerance, &fliForce,&k[0]);
                ofliIndicator(OFLI, t,r);
                if (OFLI>max_FLI) max_FLI=OFLI;
                if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(fabs(max_FLI));}

                ind=max_FLI;
			}
			else if (whichIndicator==4) {
				fastFehlberg8(t,r,dt,  tolerance, &megnoForce,&k[0]);
				omegnoIndicatorNoNormalize(MEGNO1, MEGNO2, t,r);
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(MEGNO2);}
				ind=fabs(MEGNO2);
			}
			else if (whichIndicator==5) {
				fastFehlberg8(t,r,dt,  tolerance, &modSaliForce,&k[0]);
				modSaliIndicatorNoNorm(SALI, t,r);
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(log10(SALI));}
				ind=log10(SALI);
			}
			else if (whichIndicator==6) {
				fastFehlberg8(t,r,dt, tolerance, &saliForce,&k[0]);
				//ofliIndicator(FLI, t, r);
				saliIndicator(SALI, t, r);
				//if (FLI>max_FLI) max_FLI=FLI;
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(log10(SALI));}
				//cout<<t<<" " << log10(SALI) << endl;
				ind=log10(SALI);
			}
			else if (whichIndicator==7) {
				fastFehlberg8(t,r,dt, tolerance, &fliForce,&k[0]);
				ofliIndicatorNoNorm(FLI, t,r);
				if (FLI>max_FLI) max_FLI=FLI;
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(max_FLI);}
				ind=max_FLI;
				//cout<< t << " " << FLI <<endl;
			}
			else if (whichIndicator==8) {
				fastFehlberg8(t,r,dt, tolerance, &classicMegnoForce,&k[0]);
				classicMegnoIndicator(MEGNO1, MEGNO2, t,r);
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(MEGNO2);}
				ind=MEGNO2;
			}
			else if (whichIndicator==9) {
				fastFehlberg8(t,r,dt, tolerance, &gali2Force,&k[0]);
				gali2Indicator(GALI2, t,r);
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(GALI2);}
				ind=GALI2;
			}
			else if (whichIndicator==10) {
				fastFehlberg8(t,r,dt, tolerance, &classicGali2Force,&k[0]);
				//ofliIndicator(FLI, t, r);
				classicGali2Indicator(GALI2, t, r);
				//if (FLI>max_FLI) max_FLI=FLI;
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(log10(GALI2));}
				//cout<<t<<" " << log10(SALI) << endl;
				ind=GALI2;
			}
			else if (whichIndicator==11) {
				fastFehlberg8(t,r,dt,  tolerance, &megno_omegnoForce,&k[0]);
				megno_omegnoIndicator(UMEGNO1, UMEGNO2,t,r);
				//omegnoIndicatorNoNormalize(OMEGNO1, OMEGNO2, t, r);
				//ind = ( fabs(MEGNO2 - 2)<0.01) ? MEGNO2 - (MEGNO2-fabs(OMEGNO2))*(exp(log(3)*(MEGNO2-fabs(OMEGNO2))) - 1) : MEGNO2;				
				ind = UMEGNO2;
				//cout<<t<< " " << MEGNO2<<" " << OMEGNO2 << " " << ind<<endl;
				if (evolution && t>=0.1)	{ xaxis.push_back(t); yaxis.push_back(ind);}
				
			}
			else if (whichIndicator==0) {
				pred=r;
				t_pred = t;
				fastFehlberg8(t,r,dt,  tolerance, &onlyForce,&k[0]);
				
				if (pred[3]*r[3]<0) {
					double t_interpol = t_pred - pred[3]*(t-t_pred)/(r[3]-pred[3]); 
					xaxis.push_back(pred[2] + (r[2]-pred[2]) * (t_interpol - t_pred)/(t-t_pred) ); 
				        yaxis.push_back(pred[0] + (r[0]-pred[0]) * (t_interpol - t_pred)/(t-t_pred));
					}
			}
			error1 = (jacobiIntegral(r) - initJacobiIntegral)/initJacobiIntegral;
			error2 = (secondIntegral(r) - initSecondIntegral)/initSecondIntegral;
		}
		//cout<< "Error = " << error2 << endl;
		//if (k==0)
		//{
			//r[0]=-r[0]; r[1]=-r[1];
			//n=-n;
		/*	cout<<" Half-point; " << "x=" << r[2] << " " << "y=" << r[3]<< endl;
			t=0; n=-n;
		for (i; t<tmax;  i++)
		{
			
			fehlberg8(t,r, dt,  tolerance, &force);
			double error = (jacobiIntegral(r) - initJacobiIntegral)/initJacobiIntegral;		
		}*/
		//}
	//}
	//valarray<double> diff(r-initial);
	//cout<< "Vector after integration: " << endl << "; x = " << diff[2] << "; y = " << diff[3] << endl << endl;
        cout<<"x0: " << x0 <<endl;
	cout<<"Final indicator value: " << ind<<endl;
	cout<< "Jacobi integral error: " << error1 << endl;
	
	//cout<<"SALI: " << SALI <<endl;
	//findmax(xaxis);
	//findmax(yaxis);
	clock_t t1=clock();
	CPU_time=(double)(t1-t0)/CLOCKS_PER_SEC;
	cout<<"done for "<< i << " iterations" << endl << "CPU-time: " << CPU_time << " seconds" << endl<<endl;
}
