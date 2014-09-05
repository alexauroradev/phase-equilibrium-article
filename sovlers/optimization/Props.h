#ifndef __PROPS_H__
#define __PROPS_H__

#include <cmath>

struct ValueAndDerivatives {
	double value;
	double dp;
	double dT;
	ValueAndDerivatives(double v, double dp, double dT) :
		value(v), dp(dp), dT(dT)
	{ }
	operator const double() const {
		return value;
	}
};

struct Function {
	virtual ValueAndDerivatives operator()(double p, double T) const = 0;
	virtual ~Function() { }
};

template<int np, int nT>
struct Derivative{
	Derivative(const Function &f) : f(f) { }
	double operator()(double p, double T) const;
private:
	const Function &f;
};

template<>
inline double Derivative<0, 1>::operator()(double p, double T) const {
	return f(p, T).dT;
}

template<>
inline double Derivative<1, 0>::operator()(double p, double T) const {
	return f(p, T).dp;
}

struct Kappa : public Function {
	double A, B, C;
	Kappa(double A, double B, double C) : A(A), B(B), C(C) { }
	virtual ValueAndDerivatives operator()(double p, double T) const {
		double lK = log(A) - log(p) - B / (T - C);
		return ValueAndDerivatives(lK, -1 / p, B / pow(C - T, 2));
	}
};

struct VaporizationEnthalpy : public Function {
	double D;
	double Tc;
	VaporizationEnthalpy(double D, double Tc) : D(D), Tc(Tc) { }
	virtual ValueAndDerivatives operator()(double p, double T) const {
		const double gamma = .38;
		const double eps = 1e-6;
		double DT, Dh;
		DT = Tc - T;
		double adj = sqrt(DT * DT + eps * Tc * Tc);
		DT = .5 * (DT + adj);
		Dh = D * pow(DT, gamma);
		double dlogDhdT = -gamma / adj;
		double dDhdT = dlogDhdT * Dh;

		return ValueAndDerivatives(Dh, 0, dDhdT);
	}
};

struct GasEnthalpy : public Function {
	double A, B, Qadd;
	GasEnthalpy(double A, double B, double Qadd = 0) : A(A), B(B), Qadd(Qadd) { }
	virtual ValueAndDerivatives operator()(double p, double T) const {
		double h = A * T + 0.5 * B * T * T + Qadd;
		double dhdT = A + B * T;

		return ValueAndDerivatives(h, 0, dhdT);
	}
};

struct SolidEnthalpy : public Function {
	double A, B, Qadd;
	SolidEnthalpy(double A, double B, double Qadd = 0) : A(A), B(B), Qadd(Qadd) { }
	virtual ValueAndDerivatives operator()(double p, double T) const {
		double h = A * T + 0.5 * B * T * T + Qadd;
		double dhdT = A + B * T;

		return ValueAndDerivatives(h, 0, dhdT);
	}
};

#endif
