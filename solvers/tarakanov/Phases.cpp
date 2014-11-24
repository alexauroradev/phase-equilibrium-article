#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>

#define NPH 3 //number of phases
#define NKV 4 //number of equilibrium constants
#define NC  5 //number of components

void Renorm(double *c, double *s);
void KVAL(double p, double t, double *kv);
void Phase_Equilibrium(double p, double t, double *c, double *s, double *x, double *y, double *z);
void Solver(int icase, double *c, double *s, double *kv);
void Pheq(double *c, double *s, double *kv);
void Phase_Equilibrium(double p, double t, double *c, double *s);
void Fractions(double *c, double *s, double *kv, double *x, double *y, double *z);
//void Check();

double Decision1(double *c, double *kv, double beta);
double Decision2(double *c, double *kv, double beta);

int Select(double *c, double *kv);
double StarsCorrelation(double p, double T, double A, double B, double C);

using namespace std;
int main() {
	double c[NC], s[NPH], x[NC], y[NC], z[NC];
	double norm;
	double p = 1.0e7;
	double tmin = 530.0, tmax = 690.0, dt = 0.1;
	int i;

	c[0] = 0.0;
	c[1] = 0.25;//C5H12
	c[2] = 0.25;//C10H22
	c[3] = 0.50;//H2O
	c[4] = 0.0;

	norm=0.0;

	for(i=0; i<NC; i++)
		norm += c[i];

	for(i=0; i<NC; i++)
		c[i] /= norm;

	fstream f("res.csv", ios::out);
	f << "T,s_oil,s_gas,s_wat,me,lam1,lam2" << endl;

	for(double T = tmin; T <= tmax; T += dt) {
		Phase_Equilibrium(p, T, c, s, x, y, z);
		f << T << "," << s[0] << "," << s[1] << "," << s[2] << "," << s[2] * z[3] << "," << s[0] * x[1] << "," << s[0] * x[2] << endl;
	}

	f.close();
	return 0;
}

double StarsCorrelation(double p, double T, double A, double B, double C){
	return (A / p * exp(-B/(T - C)));
}

void Renorm(double *c, double *s) {
	int i=0;
	double norm = 0.0;
	for(i=0;i<NPH;i++) {
		norm += s[i];
	}
	for(i=0;i<NPH;i++) {
		s[i] = s[i]/norm;
	}
	norm = 0.0;
	for(i=0;i<NC;i++) {
		norm += c[i];
	}
	for(i=0;i<NPH;i++) {
		s[i] = s[i]*norm;
	}
}

void KVAL(double p, double t, double *kv) {
	kv[0] = 1e-10;
	kv[1] = StarsCorrelation(p, t, 1.0026e9, 2477.07, 39.94);//C5H12
	kv[2] = StarsCorrelation(p, t, 1.1981e9, 3456.80, 78.67);//C10H22
	kv[3] = StarsCorrelation(p, t, 11.8572e9, 3816.44, 46.13);//H2O
}

void Phase_Equilibrium(double p, double t, double *c, double *s, double *x, double *y, double *z) {
	double kv[NKV];
	KVAL(p,t,kv);
	double sigma = (c[0]+c[1]+c[2]+c[3]+c[4]);
	double f0=sigma, f1=sigma, f2=sigma;
	double beta1, beta2;
	beta1 = c[4];
	beta2 = sigma;
	s[0] = 1.0;
	s[1] = 1.0;
	s[2] = 1.0;
	Solver(1,c,s,kv);

	if(c[0] != 0.0) {
		f0 = sigma-c[3]+(kv[0]+kv[3]-1.0)*s[1];
	}
	if(c[1] != 0.0) {
		f1 = sigma-c[3]+(kv[1]+kv[3]-1.0)*s[1];
	}
	if(c[2] != 0.0) {
		f2 = sigma-c[3]+(kv[2]+kv[3]-1.0)*s[1];
	}

	//    cout << f0 << " " << f1 << " " << f2 << "\n";
	if(s[0]<0.0 || s[1]<0.0 || s[2]<0.0 || f0<=0.0 || f1<=0.0 || f2<=0.0) {
		if(c[4] == 0.0 && ((kv[0]+kv[3]-1.0)*c[0]+(kv[1]+kv[3]-1.0)*c[1]+(kv[2]+kv[3]-1.0)*c[2])<=0.0) {
			Solver(3,c,s,kv);
		}
		else {
			if(c[3]/kv[3]<(c[0]/kv[0]+c[1]/kv[1]+c[2]/kv[2])) {
				Solver(4,c,s,kv);
				if(s[0]<0.0 || s[1]<0.0 || s[2]<0.0) {
					if((c[3]+c[4])==0 && ((kv[0]-1.0)*c[0]+(kv[1]-1.0)*c[1]+(kv[2]-1.0)*c[2]) <= 0.0) {
						Solver(5,c,s,kv);
					}
					else {
						Solver(6,c,s,kv);
					}
				}
			}
			else {
				Solver(2,c,s,kv);
				if(s[0]<0.0 || s[1]<0.0 || s[2]<0.0) {
					if((c[0]+c[1]+c[2]+c[4])==0.0 && kv[3]<1.0) {
						Solver(7,c,s,kv);
					}
					else {
						Solver(6,c,s,kv);
					}
				}
			}
		}
	}

	//    Solver(1,c,s,kv);
	Renorm(c,s);
	Fractions(c,s,kv,x,y,z);
}

void Fractions(double *c, double *s, double *kv, double *x, double *y, double *z) {
	int icase=0;
	int i=0;
	if(s[0] > 0.0 && s[1] > 0.0 && s[2] > 0.0) {
		icase = 1;
	}
	if(s[0] > 0.0 && s[1] > 0.0 && s[2] <= 0.0) {
		icase = 2;
	}
	if(s[0] > 0.0 && s[1] <= 0.0 && s[2] > 0.0) {
		icase = 3;
	}
	if(s[0] <= 0.0 && s[1] > 0.0 && s[2] > 0.0) {
		icase = 4;
	}
	if(s[0] <= 0.0 && s[1] <= 0.0 && s[2] > 0.0) {
		icase = 5;
	}
	if(s[0] <= 0.0 && s[1] > 0.0 && s[2] <= 0.0) {
		icase = 6;
	}
	if(s[0] > 0.0 && s[1] <= 0.0 && s[2] <= 0.0) {
		icase = 7;
	}

	for(i=0;i<NC;i++) {
		x[i]=0.0;
		y[i]=0.0;
		z[i]=0.0;
	}

	if(icase == 1) {
		x[0] = c[0]/(s[0]+kv[0]*s[1]);
		x[1] = c[1]/(s[0]+kv[1]*s[1]);
		x[2] = c[2]/(s[0]+kv[2]*s[1]);

		y[0] = kv[0]*c[0]/(s[0]+kv[0]*s[1]);
		y[1] = kv[1]*c[1]/(s[0]+kv[1]*s[1]);
		y[2] = kv[2]*c[2]/(s[0]+kv[2]*s[1]);
		y[4] = kv[3];
		y[3] = 1.0 - y[0] - y[1] - y[2] - y[4];

		z[3] = 1.0;
	}
	if(icase == 2) {
		x[0] = c[0]/(s[0]+kv[0]*s[1]);
		x[1] = c[1]/(s[0]+kv[1]*s[1]);
		x[2] = c[2]/(s[0]+kv[2]*s[1]);

		y[0] = kv[0]*c[0]/(s[0]+kv[0]*s[1]);
		y[1] = kv[1]*c[1]/(s[0]+kv[1]*s[1]);
		y[2] = kv[2]*c[2]/(s[0]+kv[2]*s[1]);
		y[3] = (1.0 - y[0] - y[1] - y[2] - y[4])*c[3]/(c[3]+c[4]);
		y[4] = (1.0 - y[0] - y[1] - y[2] - y[4])*c[4]/(c[3]+c[4]);
	}
	if(icase == 3) {
		x[0] = c[0]/(c[0]+c[1]+c[2]);
		x[1] = c[1]/(c[0]+c[1]+c[2]);
		x[2] = c[2]/(c[0]+c[1]+c[2]);

		z[3] = 1.0;
	}
	if(icase == 4) {
		y[3] = kv[3];
		y[0] = (1.0-y[3])*c[0]/(c[0]+c[1]+c[2]+c[4]);
		y[1] = (1.0-y[3])*c[1]/(c[0]+c[1]+c[2]+c[4]);
		y[2] = (1.0-y[3])*c[2]/(c[0]+c[1]+c[2]+c[4]);
		y[4] = (1.0-y[3])*c[4]/(c[0]+c[1]+c[2]+c[4]);

		z[3] = 1.0;
	}
	if(icase == 5) {
		z[3] = 1.0;
	}
	if(icase == 6) {
		y[0] = c[0]/(c[0]+c[1]+c[2]+c[3]+c[4]);
		y[1] = c[1]/(c[0]+c[1]+c[2]+c[3]+c[4]);
		y[2] = c[2]/(c[0]+c[1]+c[2]+c[3]+c[4]);
		y[3] = c[3]/(c[0]+c[1]+c[2]+c[3]+c[4]);
		y[4] = c[4]/(c[0]+c[1]+c[2]+c[3]+c[4]);
	}
	if(icase == 7) {
		x[0] = c[0]/(c[0]+c[1]+c[2]);
		x[1] = c[1]/(c[0]+c[1]+c[2]);
		x[2] = c[2]/(c[0]+c[1]+c[2]);
	}
}

double Decision1(double *c, double *kv, double beta, double &der) {
	double sigma=(c[0]+c[1]+c[2]+c[3]+c[4]);
	double res = 0;
	der = 0;

	for (int i = 0; i < 3; i++) {
		res += (kv[i] + kv[3] - 1) * c[i] / (sigma - c[3] + (kv[i] + kv[3] - 1) * beta);
		der -= pow(kv[i] + kv[3] - 1, 2) * c[i] / pow(sigma - c[3] + (kv[i] + kv[3] - 1) * beta, 2);

	}

	if(c[4] != 0.0) {
		res += c[4]/beta;
		der -= c[4] / (beta * beta);
	}

	return res;
}

double Decision2(double *c, double *kv, double beta) {
	double sigma=(c[0]+c[1]+c[2]+c[3]+c[4]);
	double res = 0.0;
	res += (kv[0]-1.0)*c[0]/(sigma + (kv[0]-1.0)*beta);
	res += (kv[1]-1.0)*c[1]/(sigma + (kv[1]-1.0)*beta);
	res += (kv[2]-1.0)*c[2]/(sigma + (kv[2]-1.0)*beta);
	if((c[3]+c[4]) != 0.0) {
		res += (c[3]+c[4])/beta;
	}
	return res;
}

void Solver(int icase, double *c, double *s, double *kv) {
	//    cout << icase << "\n";
	double sigma = (c[0]+c[1]+c[2]+c[3]+c[4]);
	int iter=0,num=100;
	double beta,beta1,beta2;
	//double pos_beta[3];
	// int i;
	beta1 = c[4];
	beta2 = sigma;
	if(icase == 1) {
		beta1 = c[4];
		beta2 = sigma;

		/*
		   for(i=0;i<3;i++) {
		   pos_beta[i]=0.0;
		   if(c[i] != 0.0) {
		   pos_beta[i] = (sigma-c[3])/(1.0-kv[i]-kv[3]);
		   }
		   }
		   beta2=pos_beta[0];
		   for(i=0;i<3;i++) {
		   if(beta2>pos_beta[i]) {
		   beta2=pos_beta[i];
		   }
		   }
		   beta1=0.0;
		   if(beta2<c[4]) {
		   beta2 = sigma;
		   beta1 = sigma;
		   }
		 */
		

/*		for(iter=0;iter<num;iter++) {
			beta = 0.5 * (beta1 + beta2);
			if(Decision1(c,kv,beta) > 0.0) {
				beta1 = beta;
			}
			else {
				beta2 = beta;
			}
		}*/

		beta = c[4];
		double v;
		int iters = 0;
		do {
			double der;
			v = Decision1(c, kv, beta, der);
			beta -= v / der;
			iters++;
			if (beta > 2. * sigma)
				break;
		} while (fabs(v) > 1e-10 && iters < 100);
		/*
		if (iters == 100) {
			cerr << "Did not converge! v = " << v << endl;
		} else 
			cout << "Converged in " << iters << " iterations with error " << v << " to value beta = " << beta << endl;
		*/
		//beta = .5 * (beta1 + beta2);
		//        cout << sigma - c[4] << "\n";
		s[1] = beta;
		s[2] = c[3] - kv[3]*beta;
		s[0] = sigma - s[1] - s[2];
	}
	if(icase == 2) {
		s[0] = 0.0;
		s[1] = (sigma-c[3])/(1.0-kv[3]);
		s[2] = sigma - s[1];
	}
	if(icase == 3) {
		s[0] = c[0]+c[1]+c[2];
		s[1] = 0.0;
		s[2] = c[3];
	}
	if(icase == 4) {
		beta1 = c[3]+c[4];
		beta2 = sigma;
		beta  = 0.5*(beta1+beta2);
		for(iter=0;iter<num;iter++) {
			if(Decision2(c,kv,beta)>0.0) {
				beta1 = beta;
			}
			else {
				beta2 = beta;
			}
			beta = 0.5*(beta1+beta2);
		}
		s[1] = beta;
		s[2] = 0.0;
		s[0] = sigma - s[1];        
	}
	if(icase == 5) {
		s[0] = 1.0;
		s[1] = 0.0;
		s[2] = 0.0;
	}
	if(icase == 6) {
		s[0] = 0.0;
		s[1] = 1.0;
		s[2] = 0.0;
	}
	if(icase == 7) {
		s[0] = 0.0;
		s[1] = 0.0;
		s[2] = 1.0;
	}
}

/*
   void Phase_Equilibrium(double p, double t, double *c, double *s) {
   double kv[NKV];
   KVAL(p,t,kv);
   Pheq(c,s,kv);
   }
 */

