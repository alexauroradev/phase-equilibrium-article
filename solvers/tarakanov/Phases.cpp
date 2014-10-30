#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "Util.h"

using namespace std;

#define NPH 3
#define NKV 4
#define NC  5

/* {{{ Forwards */
void Renorm(double *c, double *s);
void KVAL(double p, double t, double *kv);
void Solver(int icase, double *c, double *s, double *kv);
void Fractions(double *c, double *s, double *kv, double *x, double *y, double *z);
void Derivatives(double *c, double *s, double *kv, double *x, double *y, double *z, double *dkdp, double *dkdt, double vl, double vg, double vw, double *dvldx, double *dvgdy, double *dvwdz, double *dvdc);
void Simple_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z);
void PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dvdc, double *v);

void Material_Properties(double p, double t, double *kv, double *dkdp, double *dkdt, double *vl, double *vg, double *vw, double *dvldx, double *dvgdy, double *dvwdz);

double Decision1(double *c, double *kv, double beta);
double Decision2(double *c, double *kv, double beta);

int Select(double *c, double *kv);

void Matrix_Product(double *matrix1, double *matrix2, double *matrix, int r1, int c1, int r2, int c2);
void Linear_Solver(int size, double *matrix, double *rp, double *solution);
void Inverse_Matrix(int size, double *matrix, double *inv_matrix);
/* }}} */

int main() {
	double c[NC], s[NPH], x[NC], y[NC], z[NC], dvdc[NC+2], v, p;

	double norm;

	p = 1e7; /* 100 atm */

	c[0] = 0.0;
	c[1] = 0.0;
	c[2] = 0.0;
	c[3] = 1.0;
	c[4] = 1.0;

	norm = 0;
	for(int i = 0; i < NC; i++)
		norm += c[i];

	for(int i = 0; i < NC; i++)
		c[i] /= norm;

	ofstream f;
	f.open("res.csv", ios::out);
    f << "p,T,lambda,gamma,omega,v,dvdc0,dvdc1,dvdc2,dvdc3,dvdc4,dvdp,dvdT" << endl;

	for (double T = 100; T < 1200; T += 10) {
		PHEQ(p, T, c, s, x, y, z, dvdc, &v);
		f << util::format("%e,%e,%e,%e,%e,%e", p, T, s[0], s[1], s[2], v);

		/* dvdc = {dvdc0, dvdc1, ..., dvdcn, dvdp, dvdT}*/
		for (int j = 0; j < NC + 2; j++)
			f << util::format(",%e", dvdc[j]);

		f << endl;
	}

	f.close();
	return 0;
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
	double pcr[NKV], tcr[NKV];
	int i;
	pcr[0] = 1.4e6;    
	pcr[1] = 3.739e6;
	pcr[2] = 4.6e6;
	pcr[3] = 22.064e6;

	tcr[0] = 707.00;    
	tcr[1] = 425.16;
	tcr[2] = 190.55;
	tcr[3] = 647.0;

	for(i=0;i<NKV;i++) {
		kv[i] = pcr[i]/p*exp(5.37*(1.0-tcr[i]/t));
	}
	//    kv[2]=1.0e-10;
}


double Decision1(double *c, double *kv, double beta) {
	double sigma=(c[0]+c[1]+c[2]+c[3]+c[4]);
	double res = 0.0;
	res += (kv[0]+kv[3]-1.0)*c[0]/(sigma - c[3] + (kv[0]+kv[3]-1.0)*beta);
	res += (kv[1]+kv[3]-1.0)*c[1]/(sigma - c[3] + (kv[1]+kv[3]-1.0)*beta);
	res += (kv[2]+kv[3]-1.0)*c[2]/(sigma - c[3] + (kv[2]+kv[3]-1.0)*beta);
	if(c[4] != 0.0) {
		res += c[4]/beta;
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
	double sigma = (c[0]+c[1]+c[2]+c[3]+c[4]);
	int i=0,iter=0,num=100;
	double beta,beta1,beta2;
	double beta_roots[4];
	for(i=0;i<4;i++) {
		beta_roots[i] = -1.0;
	}
	if(icase == 1) {
		for(i=0;i<3;i++) {
			if(c[i] > 0.0) {
				beta_roots[i] = (sigma-c[3])/(1.0-kv[3]-kv[i]);;
			}
		}
		if(c[4] > 0.0) {
			beta_roots[3] = 0.0;
		}
		beta1 = -1.0;
		for(i=0;i<4;i++) {
			if(beta_roots[i]>=0.0 && (beta1 > beta_roots[i] || beta1 < 0.0)) {
				beta1 = beta_roots[i];
			}
		}
		beta2 = -1.0;
		for(i=0;i<4;i++) {
			if(beta1 < beta_roots[i] && (beta_roots[i] < beta2 || beta2 < beta1)) {
				beta2 = beta_roots[i];
			}
		}
		beta = 0.5*(beta1+beta2);
		if(beta1 < beta2 && beta1 >= 0.0) {    
			for(iter=0;iter<num;iter++) {
				if(Decision1(c,kv,beta)>0.0) {
					beta1 = beta;
				}
				else {
					beta2 = beta;
				}
				beta = 0.5*(beta1+beta2);
			}
			s[1] = beta;
			s[2] = c[3] - kv[3]*beta;
			s[0] = sigma - s[1] - s[2];
		}

		else {            
			s[0] = 1.0*sigma;
			s[1] = 1.0*sigma;
			s[2] = -1.0*sigma;
		}
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
		Renorm(c,s);
	}
	if(icase == 6) {
		s[0] = 0.0;
		s[1] = 1.0;
		s[2] = 0.0;
		Renorm(c,s);
	}
	if(icase == 7) {
		s[0] = 0.0;
		s[1] = 0.0;
		s[2] = 1.0;
		Renorm(c,s);
	}
}

int Select(double *c, double *kv) {
	double beta_min=0.0;
	double sigma=c[0]+c[1]+c[2]+c[3]+c[4];
	int i=0,icase=0;
	double liq,gas,wat;
	liq = c[0]+c[1]+c[2];
	gas = c[4];
	wat = c[3];
	if(liq != 0.0 && gas != 0.0 && wat != 0.0 ) {
		icase = 1;
	}
	if(liq != 0.0 && gas != 0.0 && wat == 0.0 ) {
		icase = 2;
	}
	if(liq != 0.0 && gas == 0.0 && wat != 0.0 ) {
		icase = 3;
	}
	if(liq == 0.0 && gas != 0.0 && wat != 0.0 ) {
		icase = 4;
	}
	if(liq != 0.0 && gas == 0.0 && wat == 0.0 ) {
		icase = 5;
	}
	if(liq == 0.0 && gas != 0.0 && wat == 0.0 ) {
		icase = 6;
	}
	if(liq == 0.0 && gas == 0.0 && wat != 0.0 ) {
		icase = 7;
	}

	if(icase == 1) {
		if((c[0]/kv[0]+c[1]/kv[1]+c[2]/kv[2]) < c[3]/kv[3]) {
			beta_min = (sigma-c[3])/(1.0-kv[3]);
			if(beta_min > 0.0 && (c[3]-kv[3]*beta_min) > 0.0) {
				return 2;
			}
			else {
				return 6;
			}
		}
		else {
			if(Decision2(c,kv,c[3]+c[4]) > 0.0 && Decision2(c,kv,sigma) < 0.0) {
				return 4;
			}
			else {
				return 6;
			}
		}
	}
	if(icase == 2) {
		if(Decision2(c,kv,c[3]+c[4]) > 0.0 && Decision2(c,kv,sigma) < 0.0) {
			return 4;
		}
		else {
			return 6;
		}
	}
	if(icase == 3) {
		if((kv[0]*c[0]+kv[1]*c[1]+kv[2]*c[2])/(c[0]+c[1]+c[2])+kv[3] < 1.0) {
			return 3;
		}
		else {
			if((c[0]/kv[0]+c[1]/kv[1]+c[2]/kv[2]) < c[3]/kv[3]) {
				beta_min = (sigma-c[3])/(1.0-kv[3]);
				if(beta_min > 0.0 && (c[3]-kv[3]*beta_min) > 0.0) {
					return 2;
				}
				else {
					return 6;
				}
			}
			else {
				if(Decision2(c,kv,c[3]+c[4]) > 0.0 && Decision2(c,kv,sigma) < 0.0) {
					return 4;
				}
				else {
					return 6;
				}
			}
		}
	}
	if(icase == 4) {
		beta_min = (sigma-c[3])/(1.0-kv[3]);
		if(beta_min > 0.0 && (c[3]-kv[3]*beta_min) > 0.0) {
			return 2;
		}
		else {
			return 6;
		}
	}
	if(icase == 5) {
		if(Decision2(c,kv,c[3]+c[4]) > 0.0 && Decision2(c,kv,sigma) < 0.0) {
			return 4;
		}
		else {
			if(Decision2(c,kv,c[3]+c[4]) < 0.0) {
				return 5;
			}
			if(Decision2(c,kv,sigma) > 0.0) {
				return 6;
			}
		}
	}
	if(icase == 6) {
		return 6;
	}
	if(icase == 7) {
		if(kv[3] < 1.0) {
			return 7;
		}
		else {
			return 6;
		}
	}

	((void (*)())0)();
	return -1; /* BUG */
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

void Simple_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z) {
	double kv[NKV];
	int i=0;
	KVAL(p,t,kv);
	Solver(1,c,s,kv);
	if(s[0] < 0.0 || s[1] < 0.0 || s[2] < 0.0) {
		i=Select(c,kv);
		Solver(i,c,s,kv);
	}
	Fractions(c,s,kv,x,y,z);
}

void Matrix_Product(double *matrix1, double *matrix2, double *matrix, int r1, int c1, int r2, int c2) {
	double *inner_matrix;
	int i=0,j=0,k=0;
	if(c1 == r2) {
		inner_matrix = new double[r1*c2];
		for(i=0;i<r1;i++) {
			for(j=0;j<c2;j++) {
				inner_matrix[i*c2+j] = 0.0;
				for(k=0;k<c1;k++) {
					inner_matrix[i*c2+j] += matrix1[i*c1+k]*matrix2[k*c2+j];
				}
			}
		}
		for(i=0;i<r1;i++) {
			for(j=0;j<c2;j++) {
				matrix[i*c2+j]=inner_matrix[i*c2+j];
			}
		}
		delete inner_matrix;
	}
}

void Linear_Solver(int size, double *matrix, double *rp, double *solution) {
	double *matrix1, *rp1;
	double mult=0.0;
	matrix1 = new double[size*size];
	rp1 = new double[size];
	int i=0,j=0,k=0;
	for(i=0;i<size;i++) {
		for(j=0;j<size;j++) {
			matrix1[i*size+j] = matrix[i*size+j];
		}
		rp1[i] = rp[i];
	}
	for(i=0;i<size-1;i++) {
		for(j=i+1;j<size;j++) {
			if(matrix1[i*size+i] == 0.0) {
				cout << "Fatal error in Linear_Solver()\n";
			}
			mult = matrix1[j*size+i]/matrix1[i*size+i];
			for(k=0;k<size;k++) {
				matrix1[j*size+k] -= mult*matrix1[i*size+k];
			}
			rp1[j] -= mult*rp1[i];
		}
	}
	for(i=size-1;i>0;i--) {
		for(j=i-1;j>=0;j--) {
			if(matrix1[i*size+i] == 0.0) {
				cout << "Fatal error in Linear_Solver()\n";
			}
			mult = matrix1[j*size+i]/matrix1[i*size+i];
			for(k=0;k<size;k++) {
				matrix1[j*size+k] -= mult*matrix1[i*size+k];
			}
			rp1[j] -= mult*rp1[i];
		}
	}
	for(i=0;i<size;i++) {
		solution[i] = rp1[i]/matrix1[i*size+i];
	}
	delete matrix1,rp1;
}

void Inverse_Matrix(int size, double *matrix, double *inv_matrix) {
	int i=0,j=0,k=0;
	double *rp,*inv_matrix1;
	rp = new double[size];
	inv_matrix1 = new double[size*size];
	for(i=0;i<size;i++) {
		for(j=0;j<size;j++) {
			rp[j] = 0.0;
			if(j==i) {
				rp[j]=1.0;
			}
		}
		Linear_Solver(size,matrix,rp,rp);
		for(j=0;j<size;j++) {
			inv_matrix1[j*size+i] = rp[j];
		}
	}
	for(i=0;i<size;i++) {
		for(j=0;j<size;j++) {
			inv_matrix[i*size+j] = inv_matrix1[i*size+j];
		}
	}
	delete inv_matrix1,rp;
}





void Derivatives(double *c, double *s, double *kv, double *x, double *y, double *z, double *dkdp, double *dkdt, double vl, double vg, double vw, double *dvldx, double *dvgdy, double *dvwdz, double *dvdc) {
	int icase=0;
	int i=0,j=0,k=0;
	double *dxdc,*dydc,*dzdc;
	double *dxds,*dyds,*dzds;
	double *dsdc, *dsds;
	double *dvldc, *dvgdc, *dvwdc;
	dxdc = new double[NC*(NC+2)];
	dydc = new double[NC*(NC+2)];
	dzdc = new double[NC*(NC+2)];
	dvldc = new double[NC+2];
	dvgdc = new double[NC+2];
	dvwdc = new double[NC+2];

	for(i=0;i<NC;i++) {
		for(j=0;j<(NC+2);j++) {
			dxdc[i*(NC+2)+j] = 0.0;
			dydc[i*(NC+2)+j] = 0.0;
			dzdc[i*(NC+2)+j] = 0.0;
		}
	}
	for(i=0;i<(NC+2);i++) {
		dvdc[i] = 0.0;
		dvldc[i] = 0.0;
		dvgdc[i] = 0.0;
		dvwdc[i] = 0.0;
	}
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

	if(icase == 1) {
		dxds = new double[NC*3];
		dyds = new double[NC*3];
		dzds = new double[NC*3];
		dsds = new double[3*3];
		dsdc = new double[3*(NC+2)];

		for(i=0;i<NC;i++) {
			for(j=0;j<3;j++) {
				dxds[i*3+j] = 0.0;
				dyds[i*3+j] = 0.0;
				dzds[i*3+j] = 0.0;
			}
		}
		for(i=0;i<3;i++) {
			for(j=0;j<3;j++) {
				dsds[i*3+j] = 0.0;
			}
			for(j=0;j<(NC+2);j++) {
				dsdc[i*(NC+2)+j] = 0.0;
			}
		}
		dxds[0*3+0] = -      c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dxds[0*3+1] = -kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dxds[1*3+0] = -      c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dxds[1*3+1] = -kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dxds[2*3+0] = -      c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dxds[2*3+1] = -kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);

		dyds[0*3+0] = -      kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dyds[0*3+1] = -kv[0]*kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dyds[1*3+0] = -      kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dyds[1*3+1] = -kv[1]*kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dyds[2*3+0] = -      kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dyds[2*3+1] = -kv[2]*kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dyds[3*3+1] = -kv[3]*kv[3]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);
		dyds[3*3+2] = -      kv[3]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);
		dydc[4*3+1] = -      c[4]/s[1]/s[1];

		dzds[3*3+1] = -kv[3]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);
		dzds[3*3+2] = -      c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);

		dsds[0*3+0]  = c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[0*3+0] += c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[0*3+0] += c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dsds[0*3+1]  = kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[0*3+1] += kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[0*3+1] += kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);

		dsds[1*3+0]  = kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[1*3+0] += kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[1*3+0] += kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dsds[1*3+1]  = kv[0]*kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[1*3+1] += kv[1]*kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[1*3+1] += kv[2]*kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dsds[1*3+1] += kv[3]*kv[3]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);
		dsds[1*3+1] += c[4]/s[1]/s[1];
		dsds[1*3+2]  = kv[3]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);

		dsds[2*3+1]  = kv[3]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);
		dsds[2*3+2]  = c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1]);

		dsdc[0*(NC+2)+0]  = 1.0/(s[0]+kv[0]*s[1]);
		dsdc[0*(NC+2)+1]  = 1.0/(s[0]+kv[1]*s[1]);
		dsdc[0*(NC+2)+2]  = 1.0/(s[0]+kv[2]*s[1]);
		dsdc[0*(NC+2)+5]  = -s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dsdc[0*(NC+2)+5]  -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dsdc[0*(NC+2)+5]  -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dsdc[0*(NC+2)+6]  = -s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dsdc[0*(NC+2)+6]  -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dsdc[0*(NC+2)+6]  -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];

		dsdc[1*(NC+2)+0] = kv[0]/(s[0]+kv[0]*s[1]);
		dsdc[1*(NC+2)+1] = kv[1]/(s[0]+kv[1]*s[1]);
		dsdc[1*(NC+2)+2] = kv[2]/(s[0]+kv[2]*s[1]);
		dsdc[1*(NC+2)+3] = kv[3]/(s[2]+kv[3]*s[1]);
		dsdc[1*(NC+2)+4] = 1.0/s[1];
		dsdc[1*(NC+2)+5]  = s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dsdc[1*(NC+2)+5] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dsdc[1*(NC+2)+5] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dsdc[1*(NC+2)+5] += s[2]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdp[3];
		dsdc[1*(NC+2)+6]  = s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dsdc[1*(NC+2)+6] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dsdc[1*(NC+2)+6] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];
		dsdc[1*(NC+2)+6] += s[2]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdt[3];

		dsdc[2*(NC+2)+3] = 1.0/(s[2]+kv[3]*s[1]);
		dsdc[2*(NC+2)+5] = -s[1]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdp[3];
		dsdc[2*(NC+2)+6] = -s[1]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdt[3];

		Inverse_Matrix(3,dsds,dsds);
		Matrix_Product(dsds,dsdc,dsdc,3,3,3,NC+2);

		Matrix_Product(dxds,dsdc,dxdc,NC,3,3,NC+2);
		Matrix_Product(dyds,dsdc,dydc,NC,3,3,NC+2);
		Matrix_Product(dzds,dsdc,dzdc,NC,3,3,NC+2);

		dxdc[0*(NC+2)+0] += 1.0/(s[0]+kv[0]*s[1]);
		dxdc[1*(NC+2)+1] += 1.0/(s[0]+kv[1]*s[1]);
		dxdc[2*(NC+2)+2] += 1.0/(s[0]+kv[2]*s[1]);
		dxdc[0*(NC+2)+5] -= s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dxdc[0*(NC+2)+6] -= s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dxdc[1*(NC+2)+5] -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dxdc[1*(NC+2)+6] -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dxdc[2*(NC+2)+5] -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dxdc[2*(NC+2)+6] -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];

		dydc[0*(NC+2)+0] += kv[0]/(s[0]+kv[0]*s[1]);
		dydc[1*(NC+2)+1] += kv[1]/(s[0]+kv[1]*s[1]);
		dydc[2*(NC+2)+2] += kv[2]/(s[0]+kv[2]*s[1]);
		dydc[3*(NC+2)+3] += kv[3]/(s[2]+kv[3]*s[1]);
		dydc[4*(NC+2)+4] += 1.0/s[1];
		dydc[0*(NC+2)+5] += s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dydc[0*(NC+2)+6] += s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dydc[1*(NC+2)+5] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dydc[1*(NC+2)+6] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dydc[2*(NC+2)+5] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dydc[2*(NC+2)+6] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];
		dydc[3*(NC+2)+5] += s[2]*c[3]/(s[0]+kv[3]*s[1])/(s[0]+kv[3]*s[1])*dkdp[3];
		dydc[3*(NC+2)+6] += s[2]*c[3]/(s[0]+kv[3]*s[1])/(s[0]+kv[3]*s[1])*dkdt[3];

		dzdc[3*(NC+2)+3] += 1.0/(s[2]+kv[3]*s[1]);
		dzdc[3*(NC+2)+5] -= s[1]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdp[3];
		dzdc[3*(NC+2)+6] -= s[1]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdt[3];

		for(i=0;i<(NC+2);i++) {
			dvdc[i] = dsdc[0*(NC+2)+i]*vl + dsdc[1*(NC+2)+i]*vg + dsdc[2*(NC+2)+i]*vw;
		}
		Matrix_Product(dvldx,dxdc,dvldc,1,NC,NC,NC+2);
		Matrix_Product(dvgdy,dydc,dvgdc,1,NC,NC,NC+2);
		Matrix_Product(dvwdz,dzdc,dvwdc,1,NC,NC,NC+2);
		for(i=0;i<(NC+2);i++) {
			dvdc[i] += s[0]*dvldc[i] + s[1]*dvgdc[i] + s[2]*dvwdc[i];
			if(i==5) {
				dvdc[i] += s[0]*dvldx[5] + s[1]*dvgdy[5] + s[2]*dvwdz[5];
			}
			if(i==6) {
				dvdc[i] += s[0]*dvldx[6] + s[1]*dvgdy[6] + s[2]*dvwdz[6];
			}
		}
		for(i=0;i<(NC+2);i++) {
			//            dvdc[i] = dxdc[0*(NC+2)+i];
		}        
		delete dxds,dyds,dzds,dsds,dsdc;
	}
	if(icase == 2) {
		dxds = new double[NC*2];
		dyds = new double[NC*2];
		dzds = new double[NC*2];
		dsds = new double[2*2];
		dsdc = new double[2*(NC+2)];

		for(i=0;i<NC;i++) {
			for(j=0;j<2;j++) {
				dxds[i*2+j] = 0.0;
				dyds[i*2+j] = 0.0;
				dzds[i*2+j] = 0.0;
			}
		}

		for(i=0;i<2;i++) {
			for(j=0;j<2;j++) {
				dsds[i*2+j] = 0.0;
			}
			for(j=0;j<(NC+2);j++) {
				dsdc[i*(NC+2)+j] = 0.0;
			}
		}

		dxds[0*2+0] = -      c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dxds[0*2+1] = -kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dxds[1*2+0] = -      c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dxds[1*2+1] = -kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dxds[2*2+0] = -      c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dxds[2*2+1] = -kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);

		dyds[0*2+0] = -      kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dyds[0*2+1] = -kv[0]*kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dyds[1*2+0] = -      kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dyds[1*2+1] = -kv[1]*kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dyds[2*2+0] = -      kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dyds[2*2+1] = -kv[2]*kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dyds[3*2+1] = -c[3]/s[1]/s[1];
		dyds[4*2+1] = -c[4]/s[1]/s[1];        

		dsds[0*2+0]  = c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[0*2+0] += c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[0*2+0] += c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dsds[0*2+1]  = kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[0*2+1] += kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[0*2+1] += kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);

		dsds[1*2+0]  = kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[1*2+0] += kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[1*2+0] += kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dsds[1*2+1]  = kv[0]*kv[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1]);
		dsds[1*2+1] += kv[1]*kv[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1]);
		dsds[1*2+1] += kv[2]*kv[2]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1]);
		dsds[1*2+1] += c[3]/s[1]/s[1];
		dsds[1*2+1] += c[4]/s[1]/s[1];

		dsdc[0*(NC+2)+0] = 1.0/(s[0]+kv[0]*s[1]);
		dsdc[0*(NC+2)+1] = 1.0/(s[0]+kv[1]*s[1]);
		dsdc[0*(NC+2)+2] = 1.0/(s[0]+kv[2]*s[1]);
		dsdc[0*(NC+2)+5]  = -s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dsdc[0*(NC+2)+5]  -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dsdc[0*(NC+2)+5]  -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dsdc[0*(NC+2)+6]  = -s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dsdc[0*(NC+2)+6]  -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dsdc[0*(NC+2)+6]  -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];

		dsdc[1*(NC+2)+0] = kv[0]/(s[0]+kv[0]*s[1]);
		dsdc[1*(NC+2)+1] = kv[1]/(s[0]+kv[1]*s[1]);
		dsdc[1*(NC+2)+2] = kv[2]/(s[0]+kv[2]*s[1]);
		dsdc[1*(NC+2)+3] = 1.0/s[1];
		dsdc[1*(NC+2)+4] = 1.0/s[1];
		dsdc[1*(NC+2)+5]  = s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dsdc[1*(NC+2)+5] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dsdc[1*(NC+2)+5] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dsdc[1*(NC+2)+6]  = s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dsdc[1*(NC+2)+6] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dsdc[1*(NC+2)+6] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];

		Inverse_Matrix(2,dsds,dsds);
		Matrix_Product(dsds,dsdc,dsdc,2,2,2,NC+2);

		Matrix_Product(dxds,dsdc,dxdc,NC,2,2,NC+2);
		Matrix_Product(dyds,dsdc,dydc,NC,2,2,NC+2);
		Matrix_Product(dzds,dsdc,dzdc,NC,2,2,NC+2);

		dxdc[0*(NC+2)+0] += 1.0/(s[0]+kv[0]*s[1]);
		dxdc[1*(NC+2)+1] += 1.0/(s[0]+kv[1]*s[1]);
		dxdc[2*(NC+2)+2] += 1.0/(s[0]+kv[2]*s[1]);
		dxdc[0*(NC+2)+5] -= s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dxdc[0*(NC+2)+6] -= s[1]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dxdc[1*(NC+2)+5] -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dxdc[1*(NC+2)+6] -= s[1]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dxdc[2*(NC+2)+5] -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dxdc[2*(NC+2)+6] -= s[1]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];

		dydc[0*(NC+2)+0] += kv[0]/(s[0]+kv[0]*s[1]);
		dydc[1*(NC+2)+1] += kv[1]/(s[0]+kv[1]*s[1]);
		dydc[2*(NC+2)+2] += kv[2]/(s[0]+kv[2]*s[1]);
		dydc[3*(NC+2)+3] += 1.0/s[1];
		dydc[4*(NC+2)+4] += 1.0/s[1];
		dydc[0*(NC+2)+5] += s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdp[0];
		dydc[0*(NC+2)+6] += s[0]*c[0]/(s[0]+kv[0]*s[1])/(s[0]+kv[0]*s[1])*dkdt[0];
		dydc[1*(NC+2)+5] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdp[1];
		dydc[1*(NC+2)+6] += s[0]*c[1]/(s[0]+kv[1]*s[1])/(s[0]+kv[1]*s[1])*dkdt[1];
		dydc[2*(NC+2)+5] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdp[2];
		dydc[2*(NC+2)+6] += s[0]*c[2]/(s[0]+kv[2]*s[1])/(s[0]+kv[2]*s[1])*dkdt[2];

		for(i=0;i<(NC+2);i++) {
			dvdc[i] = dsdc[0*(NC+2)+i]*vl + dsdc[1*(NC+2)+i]*vg;
		}

		Matrix_Product(dvldx,dxdc,dvldc,1,NC,NC,NC+2);
		Matrix_Product(dvgdy,dydc,dvgdc,1,NC,NC,NC+2);
		Matrix_Product(dvwdz,dzdc,dvwdc,1,NC,NC,NC+2);
		for(i=0;i<(NC+2);i++) {
			dvdc[i] += s[0]*dvldc[i] + s[1]*dvgdc[i] + s[2]*dvwdc[i];
			if(i==5) {
				dvdc[i] += s[0]*dvldx[5] + s[1]*dvgdy[5] + s[2]*dvwdz[5];
			}
			if(i==6) {
				dvdc[i] += s[0]*dvldx[6] + s[1]*dvgdy[6] + s[2]*dvwdz[6];
			}
		}        
		delete dxds,dyds,dzds,dsds,dsdc;
	}
	if(icase == 3) {
		dsdc = new double[2*(NC+2)];
		double dadc4=0.0,dbdc4=0.0,dgdc4=0.0;
		for(i=0;i<2;i++) {
			for(j=0;j<(NC+2);j++) {
				dsdc[i*(NC+2)+j]=0.0;
			}
		}        
		dsdc[0*(NC+2)+0] = 1.0;
		dsdc[0*(NC+2)+1] = 1.0;
		dsdc[0*(NC+2)+2] = 1.0;
		dsdc[1*(NC+2)+3] = 1.0;
		for(i=0;i<3;i++) {
			for(j=0;j<3;j++) {
				dxdc[i*(NC+2)+j] = -c[j]/(c[0]+c[1]+c[2])/(c[0]+c[1]+c[2]);
				if(i==j) {
					dxdc[i*(NC+2)+j] += 1.0/(c[0]+c[1]+c[2]);
				}
			}
		}
		for(i=0;i<(NC+2);i++) {
			dvdc[i] = dsdc[0*(NC+2)+i]*vl + dsdc[1*(NC+2)+i]*vw;
		}
		dbdc4 = 1.0/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]-kv[3]);
		dgdc4 = -kv[3]*dbdc4;
		dadc4 = 1.0 - dbdc4 - dgdc4;

		dxdc[0*(NC+2)+4] = -x[0]*s[0]*(dadc4+kv[0]*dbdc4);
		dxdc[1*(NC+2)+4] -= x[1]*s[0]*(dadc4+kv[1]*dbdc4);
		dxdc[2*(NC+2)+4] -= x[2]*s[0]*(dadc4+kv[2]*dbdc4);
		dydc[0*(NC+2)+4] = -kv[0]*x[0]*s[0]*(dadc4+kv[0]*dbdc4);
		dydc[1*(NC+2)+4] -= kv[1]*x[1]*s[0]*(dadc4+kv[1]*dbdc4);
		dydc[2*(NC+2)+4] -= kv[2]*x[2]*s[0]*(dadc4+kv[2]*dbdc4);
		dydc[4*(NC+2)+4] = -dydc[0*(NC+2)+4] - dydc[1*(NC+2)+4] - dydc[2*(NC+2)+4] ;
		dvdc[4] = dadc4*vl + dbdc4*vg + dgdc4*vw;

		Matrix_Product(dvldx,dxdc,dvldc,1,NC,NC,NC+2);
		Matrix_Product(dvgdy,dydc,dvgdc,1,NC,NC,NC+2);
		Matrix_Product(dvwdz,dzdc,dvwdc,1,NC,NC,NC+2);

		for(i=0;i<(NC+2);i++) {
			dvdc[i] += s[0]*dvldc[i] + s[1]*dvgdc[i] + s[2]*dvwdc[i];
			if(i==5) {
				dvdc[i] += s[0]*dvldx[5] + s[1]*dvgdy[5] + s[2]*dvwdz[5];
			}
			if(i==6) {
				dvdc[i] += s[0]*dvldx[6] + s[1]*dvgdy[6] + s[2]*dvwdz[6];
			}
		}
		delete dsdc;
	}
	if(icase==4) {
		dsdc = new double[2*(NC+2)];
		for(i=0;i<2;i++) {
			for(j=0;j<(NC+2);j++) {
				dsdc[i*(NC+2)+j]=0.0;
			}
		}
		dsdc[0*(NC+2)+0] = 1.0/(1.0-kv[3]);
		dsdc[0*(NC+2)+1] = 1.0/(1.0-kv[3]);
		dsdc[0*(NC+2)+2] = 1.0/(1.0-kv[3]);
		dsdc[0*(NC+2)+4] = 1.0/(1.0-kv[3]);
		dsdc[0*(NC+2)+5] = (c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdp[3];
		dsdc[0*(NC+2)+6] = (c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdt[3];

		dsdc[1*(NC+2)+0] = -kv[3]/(1.0-kv[3]);
		dsdc[1*(NC+2)+1] = -kv[3]/(1.0-kv[3]);
		dsdc[1*(NC+2)+2] = -kv[3]/(1.0-kv[3]);
		dsdc[1*(NC+2)+3] = 1.0;
		dsdc[1*(NC+2)+4] = -kv[3]/(1.0-kv[3]);
		dsdc[1*(NC+2)+5] = (c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdp[3];
		dsdc[1*(NC+2)+6] = (c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdt[3];

		for(i=0;i<NC;i++) {
			for(j=0;j<NC;j++) {
				if(i!=3) {
					dydc[i*(NC+2)+j] = -(1.0-kv[3])*c[i]/(c[0]+c[1]+c[2]+c[4])/(c[0]+c[1]+c[2]+c[4]);
					if(i==j) {
						dydc[i*(NC+2)+j] += (1.0-kv[3])/(c[0]+c[1]+c[2]+c[4]);
					}
				}
			}
			if(i!=3) {
				dydc[i*(NC+2)+5] = -c[i]/(c[0]+c[1]+c[2]+c[4])*dkdp[3];
				dydc[i*(NC+2)+6] = -c[i]/(c[0]+c[1]+c[2]+c[4])*dkdt[3];
			}
		}
		for(i=0;i<(NC+2);i++) {
			dvdc[i] = dsdc[0*(NC+2)+i]*vg + dsdc[1*(NC+2)+i]*vw;
		}
		Matrix_Product(dvldx,dxdc,dvldc,1,NC,NC,NC+2);
		Matrix_Product(dvgdy,dydc,dvgdc,1,NC,NC,NC+2);
		Matrix_Product(dvwdz,dzdc,dvwdc,1,NC,NC,NC+2);
		for(i=0;i<(NC+2);i++) {
			dvdc[i] += s[0]*dvldc[i] + s[1]*dvgdc[i] + s[2]*dvwdc[i];
			if(i==5) {
				dvdc[i] += s[0]*dvldx[5] + s[1]*dvgdy[5] + s[2]*dvwdz[5];
			}
			if(i==6) {
				dvdc[i] += s[0]*dvldx[6] + s[1]*dvgdy[6] + s[2]*dvwdz[6];
			}
		}                        
		delete dsdc;        
	}
	if(icase==5) {
		dvdc[4] = 1.0/(1.0-kv[3])*vg - kv[3]/(1.0-kv[3])*vw;
		for(i=0;i<3;i++) {
			dvdc[i] = dvdc[4];
			if(kv[i]+kv[3]<1.0) {
				dvdc[i]=vw;
			}
		}
		dvdc[3] = vw;
		dvdc[5] = s[2]*dvwdz[5];
		dvdc[6] = s[2]*dvwdz[6];        
	}

	if(icase==6) {
		for(i=0;i<(NC+2);i++) {
			dvdc[i] = vg;
		}
		for(i=0;i<NC;i++) {
			for(j=0;j<NC;j++) {
				dydc[i*(NC+2)+j] = -c[i]/(c[0]+c[1]+c[2]+c[3]+c[4])/(c[0]+c[1]+c[2]+c[3]+c[4]);
				if(i==j) {
					dydc[i*(NC+2)+j] += 1.0/(c[0]+c[1]+c[2]+c[3]+c[4]);
				}
			}
		}
		Matrix_Product(dvgdy,dydc,dvgdc,1,NC,NC,NC+2);
		for(i=0;i<(NC+2);i++) {
			dvdc[i] += s[0]*dvldc[i] + s[1]*dvgdc[i] + s[2]*dvwdc[i];
			if(i==5) {
				dvdc[i] += s[0]*dvldx[5] + s[1]*dvgdy[5] + s[2]*dvwdz[5];
			}
			if(i==6) {
				dvdc[i] += s[0]*dvldx[6] + s[1]*dvgdy[6] + s[2]*dvwdz[6];
			}
		}
	}
	if(icase==7) {        
		dvdc[0] = vl + dvldx[0]/(c[0]+c[1]+c[2])*(1.0-x[0]) - s[0]*dvldx[0]/(c[0]+c[1]+c[2])*x[0] - s[0]*dvldx[0]/(c[0]+c[1]+c[2])*x[0];
		dvdc[1] = vl + dvldx[1]/(c[0]+c[1]+c[2])*(1.0-x[1]) - s[0]*dvldx[1]/(c[0]+c[1]+c[2])*x[1] - s[0]*dvldx[1]/(c[0]+c[1]+c[2])*x[1];
		dvdc[2] = vl + dvldx[2]/(c[0]+c[1]+c[2])*(1.0-x[2]) - s[0]*dvldx[2]/(c[0]+c[1]+c[2])*x[2] - s[0]*dvldx[2]/(c[0]+c[1]+c[2])*x[2];
		dvdc[4]  = (1.0-1.0/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]))*vl;
		dvdc[4] += dvldx[0]*x[0]*(1.0-(kv[0]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]));
		dvdc[4] += dvldx[1]*x[1]*(1.0-(kv[1]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]));
		dvdc[4] += dvldx[2]*x[2]*(1.0-(kv[2]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]));
		dvdc[4] += 1.0/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2])*vg;
		dvdc[4] += s[1]/s[0]*dvgdy[0]*kv[0]*x[0]*(1.0-(kv[0]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]));
		dvdc[4] += s[1]/s[0]*dvgdy[1]*kv[1]*x[1]*(1.0-(kv[1]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]));
		dvdc[4] += s[1]/s[0]*dvgdy[2]*kv[2]*x[2]*(1.0-(kv[2]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]));
		dvdc[4] -= s[1]/s[0]*dvgdy[4]*(kv[0]*x[0]*(1.0-(kv[0]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]) ) + kv[1]*x[1]*(1.0-(kv[1]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2]) ) + kv[2]*x[2]*(1.0-(kv[2]+1)/(1.0-kv[0]*x[0]-kv[1]*x[1]-kv[2]*x[2])));
		dvdc[3] = dvdc[4];
		if(kv[3]+kv[0]*x[0]+kv[1]*x[1]+kv[2]*x[2] < 1.0) {
			dvdc[3] = vw;
		}
		dvdc[i] += s[0]*dvldx[5] + s[1]*dvgdy[5] + s[2]*dvwdz[5];
		dvdc[i] += s[0]*dvldx[6] + s[1]*dvgdy[6] + s[2]*dvwdz[6];
	}        
	delete[] dxdc;
	delete[] dydc;
	delete[] dzdc;
	delete[] dvldc;
	delete[] dvgdc;
	delete[] dvwdc;
}


void Material_Properties(double p, double t, double *kv, double *dkdp, double *dkdt, double *vl, double *vg, double *vw, double *dvldx, double *dvgdy, double *dvwdz) {
	int i=0;
	double pcr[NKV], tcr[NKV];
	pcr[0] = 1.4e6;    
	pcr[1] = 3.739e6;
	pcr[2] = 4.6e6;
	pcr[3] = 22.064e6;

	tcr[0] = 707.00;    
	tcr[1] = 425.16;
	tcr[2] = 190.55;
	tcr[3] = 647.0;

	for(i=0;i<NKV;i++) {
		kv[i] = pcr[i]/p*exp(5.37*(1.0-tcr[i]/t));
	}

	for(i=0;i<NKV;i++) {
		dkdp[i] = -kv[i]/p;
		dkdt[i] = 5.37*tcr[i]/t/t*kv[i];
	}
	vl[0] = 1.5;
	vg[0] = 20.0;
	vw[0] = 2.0;
	for(i=0;i<(NC+2);i++) {
		dvldx[i] = 0.0;
		dvgdy[i] = 0.0;
		dvwdz[i] = 0.0;
	}
}

void PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dvdc, double *v) {
	double kv[NKV],dkdp[NKV],dkdt[NKV];
	double vl,vg,vw;
	double dvldx[NC+2],dvgdy[NC+2],dvwdz[NC+2];
	int i=0;
	double eps = 1.0e-8;
	double norm;
	double c1[NC];
	Material_Properties(p,t,kv,dkdp,dkdt,&vl,&vg,&vw,dvldx,dvgdy,dvwdz);
	Solver(1,c,s,kv);
	if(s[0] < 0.0 || s[1] < 0.0 || s[2] < 0.0) {
		i=Select(c,kv);
		Solver(i,c,s,kv);
	}
	Fractions(c,s,kv,x,y,z);
	Derivatives(c,s,kv,x,y,z,dkdp,dkdt,vl,vg,vw,dvldx,dvgdy,dvwdz,dvdc);
	v[0] = s[0]*vl + s[1]*vg + s[2]*vw;
}
