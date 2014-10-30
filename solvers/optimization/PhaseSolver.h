#ifndef __PHASE_SOLVER_H__
#define __PHASE_SOLVER_H__

#include "LinearAlgebra.h"
#include "Util.h"

#include "Props.h"

enum Status {
	Ok = 0,
	MaximumNumberOfIterationsExceeded = 1,
	StepBecameTooSmall = 2,
	StepFactorIsZero = 3
};

class PhaseSolver {
	const int nLO, nHO, nG, nS;

	const double eps;
	double p;
	double eta;

	la::vector<double> c;
	la::matrix<double> H;
	la::vector<double> G;
	la::matrix<double> S;

	la::vector<double> x;
	la::vector<double> dx;

	Status status;
	double residual;

	double Tmin;
	double Tmax;
	double hRef;
	double sumc;

	double getAlpha() const;
	double computePhi() const;
	void computeF();
	void computeFx();
	void computeFa();

	static inline double L(const double x);
	static inline double Lx(const double x);
	static inline double Lxx(const double x);

	static inline double B(const double x, const double a, const double b, const double eps);
	static inline double Bx(const double x, const double a, const double b, const double eps);
	static inline double Bxx(const double x, const double a, const double b, const double eps);

	util::ptr_vector<Function> h;
	util::ptr_vector<Function> hGL;
	util::ptr_vector<Function> kappa;

public:

	PhaseSolver(const int nLO, const int nHO, const int nG, const int nS,
		const double eps = 1e-8);
	void setTempRange(double Tmin, double Tmax);
	void setStandartEnthalpy(const int id, Function *f);
	void setVaporizationEnthalpy(const int id, Function *f);
	void setEquilibriumConstantLog(const int id, Function *f);
	void setParams(const la::vector<double> &c, const double p, const double Tapprox,
		const double eta);
	const la::vector<double> &solve();
	const la::matrix<double> &getDerivatives();
	Status getStatus(double &res) const {
		res = residual;
		return status;
	}
	bool checkF();
	bool checkFx();
	bool checkFa();
};

#endif
