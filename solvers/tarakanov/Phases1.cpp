#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

#define NPH 3
#define NKV 4
#define NC  5

void Renorm(double *c, double *s);
void KVAL(double p, double t, double *kv);
void Solver(int icase, double *c, double *s, double *kv);
void Fractions(double *c, double *s, double *kv, double *x, double *y, double *z);
void Derivatives(double *c, double *s, double *kv, double *x, double *y, double *z, double *dkdp, double *dkdt, double vl, double vg, double vw, double *dvldx, double *dvgdy, double *dvwdz, double *dvdc);
void Simple_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z);
void PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dvdc, double *v);
void Fraction_Derivatives(double p, double t, double *c, double *s, double *kv, double *dkdp, double *dkdt, double *dXdC,  double *dYdC,  double *dZdC);
void Fraction_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dXdC, double *dYdC, double *dZdC);
void Lasy_Derivatives(double *c, double *s, double *x, double *y, double *z,
                      double *kv, double *dkdp, double *dkdt,
                      double *dXdC, double *dYdC, double *dZdC,
                      double *dvldx, double *dvgdy, double *dvwdz,
                      double *dvldp, double *dvgdp, double *dvwdp,
                      double *dvldt, double *dvgdt, double *dvwdt,
                      double *v, double *dvdc);

void Lasy_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dXdC, double *dYdC, double *dZdC, double *v, double *dvdc);
     
void Material_Properties(double p, double t, double *kv, double *dkdp, double *dkdt, double *vl, double *vg, double *vw, double *dvldx, double *dvgdy, double *dvwdz);
void Fraction_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dXdC, double *dYdC, double *dZdC);
double Decision1(double *c, double *kv, double beta);
double Decision2(double *c, double *kv, double beta);

int Select(double *c, double *kv);

void Matrix_Product(double *matrix1, double *matrix2, double *matrix, int r1, int c1, int r2, int c2);
void Linear_Solver(int size, double *matrix, double *rp, double *solution);
void Inverse_Matrix(int size, double *matrix, double *inv_matrix);

void Check();
/*
int main() {
    Check();
    return 0;
    }
*/
int main() {
    FILE *fp;
    fp = fopen("Diagram.txt","w");
//	fp = stdout;
    double c[NC],s[NPH],x[NC],y[NC],z[NC], dvdc[NC+2], dXdC[NC*(NC+2)], dYdC[NC*(NC+2)], dZdC[NC*(NC+2)], v[0];
    double norm=0.0, p=1.0e7, tmin=300.0, tmax = 1200.0, dt=0.0;
    double pmin = 1.0e7, pmax = 1.0e7, dp=0.0;
    double vnext=0.0,vcurrent=0.0;
    int i=0,num=1000;
    dt = (tmax-tmin)/num;
    dp = (pmax-pmin)/num;
    for(i=0;i<NC;i++) {
        c[i] = 1.0;
        }
    c[0] = 1.0;
    c[1] = 1.0;
    c[2] = 1.0;
    c[3] = 6.0;
    c[4] = 1.0;
    norm=0.0;
    for(i=0;i<NC;i++) {
        norm += c[i];
        }
    for(i=0;i<NC;i++) {
        c[i] = c[i]/norm;
        }
    dp = 0.0;
    dt = 0.0;

	for (double T = 200; T < 1200; T += 1) {
		Lasy_PHEQ(1e7, T, c, s, x, y, z, dXdC, dYdC, dZdC, v, dvdc);
		fprintf(fp, "%e %e %e %e %e\n", T, s[0] * x[0], s[0] * x[1], s[0] * x[2], s[2]);
	}

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
//    cout << icase << "\n";
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
        beta1 = 1.0;
        for(i=0;i<3;i++) {
            if(c[i] <= 0.0) {
                beta_roots[i] = 1.0;
                }
            }
        if(c[4] <= 0.0) {
            beta_roots[3] = 1.0;
            }
        for(i=0;i<4;i++) {
            if(beta_roots[i]<=0.0 && (beta1 < beta_roots[i] || beta1 > 0.0)) {
                beta1 = beta_roots[i];
                }
            }
        beta2 = -1.0;
        for(i=0;i<3;i++) {
            if(c[i] <= 0.0) {
                beta_roots[i] = -1.0;
                }
            }
        if(c[4] <= 0.0) {
            beta_roots[3] = -1.0;
            }

        for(i=0;i<4;i++) {
            if((0.0 < beta_roots[i] && beta_roots[i] < beta2) || (beta2 < 0.0)) {
                beta2 = beta_roots[i];
                }
            }
//        cout << beta1 << " " << beta2 << "\n";
//        cout << beta_roots[0] << "\n";
//        cout << beta_roots[1] << "\n";
//        cout << beta_roots[2] << "\n";
//        cout << beta_roots[3] << "\n";

//        cout << beta1 << " " << beta2 << "\n";
        beta = 0.5*(beta1+beta2);
        if(beta1 <= 0.0 && beta2 > 0.0) {    
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
//            cout << Decision1(c,kv,beta) << "\n";
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
//        cout << "I am here\n";
//        cout << ((kv[0]*c[0]+kv[1]*c[1]+kv[2]*c[2])/(c[0]+c[1]+c[2])+kv[3]) << "\n";
        if(((kv[0]*c[0]+kv[1]*c[1]+kv[2]*c[2])/(c[0]+c[1]+c[2])+kv[3]) < 1.0) {
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
        if(Decision2(c,kv,0.0) > 0.0 && Decision2(c,kv,sigma) < 0.0) {
            return 4;
            }
        else {
            if(Decision2(c,kv,0.0) < 0.0) {
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
    }

void Fractions(double *c, double *s, double *kv, double *x, double *y, double *z) {
    int icase=0;
    int i=0;
    double eps = 1.0e-6;
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
        y[3] = kv[3];
        y[4] = 1.0 - y[0] - y[1] - y[2] - y[3];
        
        z[3] = 1.0;
        if((x[0]+x[1]+x[2]+x[3]+x[4]-1.0) > eps || (1.0-(x[0]+x[1]+x[2]+x[3]+x[4])) > eps) {
            cout << "An error in x[i] has appeared \n";
            }
        if((y[0]+y[1]+y[2]+y[3]+y[4]-1.0) > eps || (1.0-(y[0]+y[1]+y[2]+y[3]+y[4])) > eps) {
            cout << "An error in y[i] has appeared \n";
            }
        if((z[0]+z[1]+z[2]+z[3]+z[4]-1.0) > eps || (1.0-(z[0]+z[1]+z[2]+z[3]+z[4])) > eps) {
            cout << "An error in z[i] has appeared \n";
            }
        }
    if(icase == 2) {
        x[0] = c[0]/(s[0]+kv[0]*s[1]);
        x[1] = c[1]/(s[0]+kv[1]*s[1]);
        x[2] = c[2]/(s[0]+kv[2]*s[1]);
        
        y[0] = kv[0]*c[0]/(s[0]+kv[0]*s[1]);
        y[1] = kv[1]*c[1]/(s[0]+kv[1]*s[1]);
        y[2] = kv[2]*c[2]/(s[0]+kv[2]*s[1]);
        if((c[3]+c[4]) > 0.0) {
            y[3] = (1.0 - y[0] - y[1] - y[2])*c[3]/(c[3]+c[4]);
            y[4] = (1.0 - y[0] - y[1] - y[2])*c[4]/(c[3]+c[4]);
            }
        if((x[0]+x[1]+x[2]+x[3]+x[4]-1.0) > eps || (1.0-(x[0]+x[1]+x[2]+x[3]+x[4])) > eps) {
            cout << "An error in x[i] has appeared \n";
            }
        if((y[0]+y[1]+y[2]+y[3]+y[4]-1.0) > eps || (1.0-(y[0]+y[1]+y[2]+y[3]+y[4])) > eps) {
            cout << "An error in y[i] has appeared \n";
            }
        }
    if(icase == 3) {
        if((c[0]+c[1]+c[2]) > 0.0) {
            x[0] = c[0]/(c[0]+c[1]+c[2]);
            x[1] = c[1]/(c[0]+c[1]+c[2]);
            x[2] = c[2]/(c[0]+c[1]+c[2]);
            }
        if((x[0]+x[1]+x[2]+x[3]+x[4]-1.0) > eps || (1.0-(x[0]+x[1]+x[2]+x[3]+x[4])) > eps) {
            cout << "An error in x[i] has appeared \n";
            }

        z[3] = 1.0;
        if((z[0]+z[1]+z[2]+z[3]+z[4]-1.0) > eps || (1.0-(z[0]+z[1]+z[2]+z[3]+z[4])) > eps) {
            cout << "An error in z[i] has appeared \n";
            }
        }
    if(icase == 4) {
        y[3] = kv[3];
        if((c[0]+c[1]+c[2]+c[4])>0.0) {
            y[0] = (1.0-y[3])*c[0]/(c[0]+c[1]+c[2]+c[4]);
            y[1] = (1.0-y[3])*c[1]/(c[0]+c[1]+c[2]+c[4]);
            y[2] = (1.0-y[3])*c[2]/(c[0]+c[1]+c[2]+c[4]);
            y[4] = (1.0-y[3])*c[4]/(c[0]+c[1]+c[2]+c[4]);
            }
        z[3] = 1.0;
        if((y[0]+y[1]+y[2]+y[3]+y[4]-1.0) > eps || (1.0-(y[0]+y[1]+y[2]+y[3]+y[4])) > eps) {
            cout << "An error in y[i] has appeared \n";
            }
        if((z[0]+z[1]+z[2]+z[3]+z[4]-1.0) > eps || (1.0-(z[0]+z[1]+z[2]+z[3]+z[4])) > eps) {
            cout << "An error in z[i] has appeared \n";
            }
        }
    if(icase == 5) {
        z[3] = 1.0;
        if((z[0]+z[1]+z[2]+z[3]+z[4]-1.0) > eps || (1.0-(z[0]+z[1]+z[2]+z[3]+z[4])) > eps) {
            cout << "An error in z[i] has appeared \n";
            }
        }
    if(icase == 6) {
        if((c[0]+c[1]+c[2]+c[3]+c[4])>0.0) {
            y[0] = c[0]/(c[0]+c[1]+c[2]+c[3]+c[4]);
            y[1] = c[1]/(c[0]+c[1]+c[2]+c[3]+c[4]);
            y[2] = c[2]/(c[0]+c[1]+c[2]+c[3]+c[4]);
            y[3] = c[3]/(c[0]+c[1]+c[2]+c[3]+c[4]);
            y[4] = c[4]/(c[0]+c[1]+c[2]+c[3]+c[4]);
            }
        if((y[0]+y[1]+y[2]+y[3]+y[4]-1.0) > eps || (1.0-(y[0]+y[1]+y[2]+y[3]+y[4])) > eps) {
            cout << "An error in y[i] has appeared \n";
            }
        }
    if(icase == 7) {
        if((c[0]+c[1]+c[2])>0.0) {
            x[0] = c[0]/(c[0]+c[1]+c[2]);
            x[1] = c[1]/(c[0]+c[1]+c[2]);
            x[2] = c[2]/(c[0]+c[1]+c[2]);
            }
        if((x[0]+x[1]+x[2]+x[3]+x[4]-1.0) > eps || (1.0-(x[0]+x[1]+x[2]+x[3]+x[4])) > eps) {
            cout << "An error in x[i] has appeared \n";
            }
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


void Material_Properties(double p, double t,
                         double *kv, double *dkdp, double *dkdt,
                         double *dvldx, double *dvgdy, double *dvwdz,
                         double *dvldp, double *dvgdp, double *dvwdp,
                         double *dvldt, double *dvgdt, double *dvwdt) {
    int i=0;
    double pcr[NKV], tcr[NKV];
    pcr[0] = 15e5;    
    pcr[1] = 20.8e5;
    pcr[2] = 33.3e5;
    pcr[3] = 22.064e6;

    tcr[0] = 707.00;    
    tcr[1] = 617.6;
    tcr[2] = 469.6;
    tcr[3] = 647.0;

    for(i=0;i<NKV;i++) {
        kv[i] = pcr[i]/p*exp(5.37*(1.0-tcr[i]/t));
        }

    for(i=0;i<NKV;i++) {
        dkdp[i] = -kv[i]/p;
        dkdt[i] = 5.37*tcr[i]/t/t*kv[i];
        }
    for(i=0;i<NC;i++) {
        dvldx[i] = 1.5;
        dvgdy[i] = 20.0;
        dvwdz[i] = 2.0;
        dvldp[i] = 0.0;
        dvgdp[i] = 0.0;
        dvwdp[i] = 0.0;
        dvldt[i] = 0.0;
        dvgdt[i] = 0.0;
        dvwdt[i] = 0.0;
        }    
    }


void Lasy_Derivatives(double *c, double *s, double *x, double *y, double *z,
                      double *kv, double *dkdp, double *dkdt,
                      double *dXdC, double *dYdC, double *dZdC,
                      double *dvldx, double *dvgdy, double *dvwdz,
                      double *dvldp, double *dvgdp, double *dvwdp,
                      double *dvldt, double *dvgdt, double *dvwdt,
                      double *v, double *dvdc) {
    int icase=0;
    int i=0,j=0,k=0;
    double dxdc[NC*(NC+2)],dydc[NC*(NC+2)],dzdc[NC*(NC+2)];

    for(i=0;i<NC;i++) {
        for(j=0;j<(NC+2);j++) {
            dxdc[i*(NC+2)+j] = 0.0;
            dydc[i*(NC+2)+j] = 0.0;
            dzdc[i*(NC+2)+j] = 0.0;
            dXdC[i*(NC+2)+j] = 0.0;
            dYdC[i*(NC+2)+j] = 0.0;
            dZdC[i*(NC+2)+j] = 0.0;
            }
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
    if(s[0] > 0.0 && s[1] <= 0.0 && s[2] <= 0.0) {
        icase = 5;
        }
    if(s[0] <= 0.0 && s[1] > 0.0 && s[2] <= 0.0) {
        icase = 6;
        }
    if(s[0] <= 0.0 && s[1] <= 0.0 && s[2] > 0.0) {
        icase = 7;
        }
//    cout << icase << "\n";
    if(icase == 1) {
        double dxds[NC*3], dyds[NC*3], dzds[NC*3], dsds[3*3], dsdc[3*(NC+2)];
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
        dyds[4*3+1] = -      c[4]/s[1]/s[1];

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
        Matrix_Product(dsds,dsdc,dsdc, 3,3,3,NC+2);
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
        dydc[3*(NC+2)+5] += s[2]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdp[3];
        dydc[3*(NC+2)+6] += s[2]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdt[3];

        dzdc[3*(NC+2)+3] += 1.0/(s[2]+kv[3]*s[1]);
        dzdc[3*(NC+2)+5] -= s[1]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdp[3];
        dzdc[3*(NC+2)+6] -= s[1]*c[3]/(s[2]+kv[3]*s[1])/(s[2]+kv[3]*s[1])*dkdt[3];

        for(i=0;i<NC;i++) {
            for(j=0;j<(NC+2);j++) {
                dXdC[i*(NC+2)+j] = x[i]*dsdc[0*(NC+2)+j] + s[0]*dxdc[i*(NC+2)+j];
                dYdC[i*(NC+2)+j] = y[i]*dsdc[1*(NC+2)+j] + s[1]*dydc[i*(NC+2)+j];
                dZdC[i*(NC+2)+j] = z[i]*dsdc[2*(NC+2)+j] + s[2]*dzdc[i*(NC+2)+j];
                }
            }
        }
    if(icase == 2) {
        double dxds[NC*2], dyds[NC*2], dzds[NC*2], dsds[2*2], dsdc[2*(NC+2)];
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

        for(i=0;i<NC;i++) {
            for(j=0;j<(NC+2);j++) {
                dXdC[i*(NC+2)+j] = x[i]*dsdc[0*(NC+2)+j] + s[0]*dxdc[i*(NC+2)+j];
                dYdC[i*(NC+2)+j] = y[i]*dsdc[1*(NC+2)+j] + s[1]*dydc[i*(NC+2)+j];
                dZdC[i*(NC+2)+j] = 0.0;
                }
            }
        }

    if(icase == 3) {
        dXdC[0*(NC+2)+0] = 1.0;
        dXdC[1*(NC+2)+1] = 1.0;
        dXdC[2*(NC+2)+2] = 1.0;
        dZdC[3*(NC+2)+3] = 1.0;

        double dbdc4=0.0;
        dbdc4 = (c[0]+c[1]+c[2])/((1.0-kv[3])*(c[0]+c[1]+c[2])-kv[0]*c[0]-kv[1]*c[1]-kv[2]*c[2]);
        dXdC[0*(NC+2)+4] = -kv[0]*c[0]/(c[0]+c[1]+c[2])*dbdc4;
        dXdC[1*(NC+2)+4] = -kv[1]*c[1]/(c[0]+c[1]+c[2])*dbdc4;
        dXdC[2*(NC+2)+4] = -kv[2]*c[2]/(c[0]+c[1]+c[2])*dbdc4;
        dZdC[3*(NC+2)+4] = -kv[3]*dbdc4;
        dYdC[0*(NC+2)+4] = kv[0]*c[0]/(c[0]+c[1]+c[2])*dbdc4;
        dYdC[1*(NC+2)+4] = kv[1]*c[1]/(c[0]+c[1]+c[2])*dbdc4;
        dYdC[2*(NC+2)+4] = kv[2]*c[2]/(c[0]+c[1]+c[2])*dbdc4;
        dYdC[3*(NC+2)+4] = kv[3]*dbdc4;
        dYdC[4*(NC+2)+4] = 1.0;
        }
    if(icase == 4) {
        dYdC[0*(NC+2)+0] = 1.0;
        dYdC[1*(NC+2)+1] = 1.0;
        dYdC[2*(NC+2)+2] = 1.0;
        dYdC[4*(NC+2)+4] = 1.0;

        dYdC[3*(NC+2)+0] = kv[3]/(1.0-kv[3]);
        dYdC[3*(NC+2)+1] = kv[3]/(1.0-kv[3]);
        dYdC[3*(NC+2)+2] = kv[3]/(1.0-kv[3]);
        dYdC[3*(NC+2)+4] = kv[3]/(1.0-kv[3]);
        dYdC[3*(NC+2)+3] = 1.0;
        dYdC[3*(NC+2)+5] = (c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdp[3];
        dYdC[3*(NC+2)+6] = (c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdt[3];

        dZdC[3*(NC+2)+3] = 1.0;
        dZdC[3*(NC+2)+0] = -kv[3]/(1.0-kv[3]);
        dZdC[3*(NC+2)+1] = -kv[3]/(1.0-kv[3]);
        dZdC[3*(NC+2)+2] = -kv[3]/(1.0-kv[3]);
        dZdC[3*(NC+2)+4] = -kv[3]/(1.0-kv[3]);
        dZdC[3*(NC+2)+5] = -(c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdp[3];
        dZdC[3*(NC+2)+6] = -(c[0]+c[1]+c[2]+c[4])/(1.0-kv[3])/(1.0-kv[3])*dkdt[3];
        }
    if(icase == 5) {
        dXdC[0*(NC+2)+0] = 1.0;
        dXdC[1*(NC+2)+1] = 1.0;
        dXdC[2*(NC+2)+2] = 1.0;

        double dadc4=0.0, dbdc4=0.0;
        dbdc4 = (c[0]+c[1]+c[2])/(c[0]+c[1]+c[2]-kv[0]*c[0]-kv[1]*c[1]-kv[2]*c[2]);
        dadc4 = 1.0-dbdc4;
        dXdC[0*(NC+2)+4] = kv[0]*c[0]/(c[0]+c[1]+c[2])*dbdc4;
        dXdC[1*(NC+2)+4] = kv[1]*c[1]/(c[0]+c[1]+c[2])*dbdc4;
        dXdC[2*(NC+2)+4] = kv[2]*c[2]/(c[0]+c[1]+c[2])*dbdc4;
        dYdC[0*(NC+2)+4] = kv[0]*c[0]/(c[0]+c[1]+c[2])*dbdc4;
        dYdC[1*(NC+2)+4] = kv[1]*c[1]/(c[0]+c[1]+c[2])*dbdc4;
        dYdC[2*(NC+2)+4] = kv[2]*c[2]/(c[0]+c[1]+c[2])*dbdc4;
        if((kv[0]*c[0]+kv[1]*c[1]+kv[2]*c[2])/(c[0]+c[1]+c[2])+kv[3] < 1.0) {
            dZdC[3*(NC+2)+3] = 1.0;
            }
        else {
            dXdC[0*(NC+2)+3] = kv[0]*c[0]/(c[0]+c[1]+c[2])*dbdc4;
            dXdC[1*(NC+2)+3] = kv[1]*c[1]/(c[0]+c[1]+c[2])*dbdc4;
            dXdC[2*(NC+2)+3] = kv[2]*c[2]/(c[0]+c[1]+c[2])*dbdc4;
            dYdC[0*(NC+2)+3] = kv[0]*c[0]/(c[0]+c[1]+c[2])*dbdc4;
            dYdC[1*(NC+2)+3] = kv[1]*c[1]/(c[0]+c[1]+c[2])*dbdc4;
            dYdC[2*(NC+2)+3] = kv[2]*c[2]/(c[0]+c[1]+c[2])*dbdc4;
            }
        }
    if(icase == 6) {
        dYdC[0*(NC+2)+0] = 1.0;
        dYdC[1*(NC+2)+1] = 1.0;
        dYdC[2*(NC+2)+2] = 1.0;
        dYdC[3*(NC+2)+3] = 1.0;
        dYdC[4*(NC+2)+4] = 1.0;
        }
    if(icase == 7) {
        dZdC[3*(NC+2)+3] = 1.0;

        dZdC[3*(NC+2)+4] = -kv[3]/(1.0-kv[3]);
        dYdC[3*(NC+2)+4] = kv[3]/(1.0-kv[3]);
        if(kv[0]+kv[3]<1.0) {
            dXdC[0*(NC+2)+0] = 1.0;
            }
        else {
            dYdC[0*(NC+2)+0] = 1.0;
            dZdC[3*(NC+2)+0] = -kv[3]/(1.0-kv[3]);
            dYdC[3*(NC+2)+0] = kv[3]/(1.0-kv[3]);
            }
        if(kv[1]+kv[3]<1.0) {
            dXdC[1*(NC+2)+1] = 1.0;
            }
        else {
            dYdC[1*(NC+2)+1] = 1.0;
            dZdC[3*(NC+2)+1] = -kv[3]/(1.0-kv[3]);
            dYdC[3*(NC+2)+1] = kv[3]/(1.0-kv[3]);
            }
        if(kv[2]+kv[3]<1.0) {
            dXdC[2*(NC+2)+2] = 1.0;
            }
        else {
            dYdC[2*(NC+2)+2] = 1.0;
            dZdC[3*(NC+2)+2] = -kv[3]/(1.0-kv[3]);
            dYdC[3*(NC+2)+2] = kv[3]/(1.0-kv[3]);
            }
        }
    *v = 0.0;
    for(i=0;i<NC;i++) {
        *v += dvldx[i]*s[0]*x[i] + dvgdy[i]*s[1]*y[i] + dvwdz[i]*s[2]*z[i];
        }
    for(i=0;i<(NC+2);i++) {
        dvdc[i] = 0.0;
        }
    for(i=0;i<(NC+2);i++) {
        for(j=0;j<NC;j++) {
            dvdc[i] += dvldx[j]*dXdC[j*(NC+2)+i];
            dvdc[i] += dvgdy[j]*dYdC[j*(NC+2)+i];
            dvdc[i] += dvwdz[j]*dZdC[j*(NC+2)+i];
            }
        if(i==NC) {
            for(j=0;j<NC;j++) {
                dvdc[i] += dvldp[j]*s[0]*x[j];
                dvdc[i] += dvgdp[j]*s[1]*y[j];
                dvdc[i] += dvwdp[j]*s[2]*z[j];
                }
            }
        if(i==(NC+1)) {
            for(j=0;j<NC;j++) {
                dvdc[i] += dvldt[j]*s[0]*x[j];
                dvdc[i] += dvgdt[j]*s[1]*y[j];
                dvdc[i] += dvwdt[j]*s[2]*z[j];
                }
            }
        }
    }

void Lasy_PHEQ(double p, double t, double *c, double *s, double *x, double *y, double *z, double *dXdC, double *dYdC, double *dZdC, double *v, double *dvdc) {
    double kv[NKV],dkdp[NKV],dkdt[NKV];
    double dvldx[NC],dvgdy[NC],dvwdz[NC];
    double dvldp[NC], dvgdp[NC], dvwdp[NC];
    double dvldt[NC], dvgdt[NC], dvwdt[NC];
    int i=0;
    Material_Properties(p,t,kv,dkdp,dkdt,dvldx,dvgdy,dvwdz,dvldp,dvgdp,dvwdp,dvldt,dvgdt,dvwdt);
    Solver(1,c,s,kv);
    if(s[0] < 0.0 || s[1] < 0.0 || s[2] < 0.0) {
        i=Select(c,kv);
        Solver(i,c,s,kv);
        }
    Fractions(c,s,kv,x,y,z);    
    Lasy_Derivatives(c,s,x,y,z,kv,dkdp,dkdt,dXdC,dYdC,dZdC,dvldx,dvgdy,dvwdz,dvldp,dvgdp,dvwdp,dvldt,dvgdt,dvwdt,v,dvdc);    
    }
