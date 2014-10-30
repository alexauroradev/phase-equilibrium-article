#include "PhaseSolver.h"
#include <cmath>

#include <iostream>

PhaseSolver::PhaseSolver(const int nLO, const int nHO, const int nG, const int nS,
	const double eps)
:
	nLO(nLO), nHO(nHO), nG(nG), nS(nS),
	eps(eps),
	c(nLO + nHO + nG + nS + 1), H(nLO + 2, nLO + 2), G(nLO + 2), S(nLO + 2, 2 + c.N()),
	x(nLO + 2), dx(nLO + 2),
	h(nLO + nHO + nG + nS + 1), hGL(nLO + 1), kappa(nLO + 1)
{
}

void PhaseSolver::setTempRange(double Tmin, double Tmax) {
	this->Tmin = Tmin;
	this->Tmax = Tmax;
}

void PhaseSolver::setStandartEnthalpy(const int id, Function *f) {
	h.put(id, f);
}

void PhaseSolver::setVaporizationEnthalpy(const int id, Function *f) {
	hGL.put(id, f);
}

void PhaseSolver::setEquilibriumConstantLog(const int id, Function *f) {
	kappa.put(id, f);
}

void PhaseSolver::setParams(const la::vector<double> &c, const double p, const double Tapprox,
	const double eta)
{
	this->p = p;
	this->eta = eta;

	const int n = nLO + nHO + nG + nS + 1;
	assert(c.N() == n);

	sumc = 0;
	for (int i = 0; i < c.N(); i++)
		sumc += c[i];
	sumc += (nLO + 1) * eps;

	for (int i = 0; i <= nLO; i++)
		this->c[i] = (c[i] + eps) / sumc;
	for (int i = nLO + 1; i < c.N(); i++)
		this->c[i] = c[i] / sumc;

	for (int i = 0; i <= nLO; i++)
		x[i] = .5 * this->c[i];
	x[nLO + 1] = Tapprox;

	double C = 0;
	for (int i = 0; i < c.N(); i++) {
		double v = fabs(Derivative<0, 1>(h[i])(p, Tapprox));
		if (v > C)
			C = v;
	}

	hRef = C * Tapprox;
}

double PhaseSolver::getAlpha() const {
	double safe = 0.9;
	double amax = 1 / safe;

	for (int i = 0; i <= nLO; i++) {
		if (x[i] + amax * dx[i] > c[i])
			amax = (c[i] - x[i]) / dx[i];
		if (x[i] + amax * dx[i] < 0)
			amax = (0 - x[i]) / dx[i];
	}

	const int i = nLO + 1;
	if (x[i] + amax * dx[i] > Tmax)
		amax = (Tmax - x[i]) / dx[i];
	if (x[i] + amax * dx[i] < Tmin)
		amax = (Tmin - x[i]) / dx[i];

	return std::min(1., amax * safe);
}

const la::vector<double> &PhaseSolver::solve() {
	int itmax = 50;
	do {
		computeF();
		computeFx();

		H.lu_factorize();
		dx = -G;
		H.lu_solve(dx);
		double alpha = getAlpha();

		dx *= alpha;
		x += dx;

		G[nLO + 1] /= hRef;

		residual = G.norm();

		if (residual < 1e-6) {
			status = Ok;
			break;
		}

		double dz = 0;
		for (int i = 0; i <= nLO; i++)
			dz += fabs(dx[i]);
		dz /= sumc;
		dz += dx[nLO + 1] / x[nLO + 1];

		if (dz < 1e-14) {
			if (alpha < 1e-6)
				status = StepFactorIsZero;
			else
				status = StepBecameTooSmall;
		}

		itmax--;

		if (!itmax) {
			status = MaximumNumberOfIterationsExceeded;
			break;
		}
	} while (true);

	return x;
}

const la::matrix<double> &PhaseSolver::getDerivatives() {
	computeFx();
	computeFa();

	H.lu_factorize();
	S.negate();
	H.lu_solve(S);

	return S;
}

double PhaseSolver::L(const double x) { return x * log(x); }
double PhaseSolver::Lx(const double x) { return 1 + log(x); }
double PhaseSolver::Lxx(const double x) { return 1. / x; }

double PhaseSolver::B(const double x, const double a, const double b, const double eps) {
	return -eps * (log(x - a) + log(b - x));
}
double PhaseSolver::Bx(const double x, const double a, const double b, const double eps) {
	return -eps * (1 / (x - a) - 1 / (b - x));
}
double PhaseSolver::Bxx(const double x, const double a, const double b, const double eps) {
	return eps * (1 / ((x - a) * (x - a)) + 1 / ((b - x) * (b - x)));
}

double PhaseSolver::computePhi() const {
	double phi = 0;
	double T = x[nLO + 1];

	double lam = 0;
	for (int i = 1; i <= nLO; i++)
		lam += x[i];
	for (int i = nLO + 1; i <= nLO + nHO; i++)
		lam += c[i];

	double gam = c[0] - x[0];
	for (int i = 1; i <= nLO; i++)
		gam += c[i] - x[i];
	for (int i = nLO + nHO + 1; i <= nLO + nHO + nG; i++)
		gam += c[i];

	for (int i = 0; i <= nLO; i++)
		phi += x[i] * kappa[i](p, T);

	for (int i = 1; i <= nLO; i++)
		phi += L(x[i]);
	phi -= L(lam);

	for (int i = 0; i <= nLO; i++)
		phi += L(c[i] - x[i]);
	phi -= L(gam);

	phi += B(x[0], 0, c[0], eps);
	for (int i = 1; i <= nLO; i++)
		phi += B(x[i], 0, c[i], eps);

	return phi;
}

void PhaseSolver::computeF() {
	double T = x[nLO + 1];

	double lam = 0;
	for (int i = 1; i <= nLO; i++)
		lam += x[i];
	for (int i = nLO + 1; i <= nLO + nHO; i++)
		lam += c[i];

	double gam = c[0] - x[0];
	for (int i = 1; i <= nLO; i++)
		gam += c[i] - x[i];

	for (int i = nLO + nHO + 1; i <= nLO + nHO + nG; i++)
		gam += c[i];

	double Lx_gam = Lx(gam);
	double Lx_lam = Lx(lam);

	G[0] = kappa[0](p, T) - Lx(c[0] - x[0]) + Lx_gam + Bx(x[0], 0, c[0], eps);
	for (int i = 1; i <= nLO; i++)
		G[i] = kappa[i](p, T) - Lx(c[i] - x[i]) + Lx_gam - Lx_lam + Lx(x[i]) + Bx(x[i], 0, c[i], eps);

	double ent = 0;
	for (int i = 0; i <= nLO + nHO + nG + nS; i++)
		ent += c[i] * h[i](p, T);
	for (int i = 0; i <= nLO; i++)
		ent -= x[i] * hGL[i](p, T);

	G[1 + nLO] = eta - ent;
}

void PhaseSolver::computeFx() {
	double T = x[nLO + 1];

	double lam = 0;
	for (int i = 1; i <= nLO; i++)
		lam += x[i];
	for (int i = nLO + 1; i <= nLO + nHO; i++)
		lam += c[i];

	double gam = c[0] - x[0];
	for (int i = 1; i <= nLO; i++)
		gam += c[i] - x[i];

	for (int i = nLO + nHO + 1; i <= nLO + nHO + nG; i++)
		gam += c[i];

	double Lxx_gam = Lxx(gam);
	double Lxx_lam = Lxx(lam);

	for (int i = 0; i <= nLO; i++)
		H(0, i) = H(i, 0) = -Lxx_gam;
	H(0, 0) += Lxx(c[0] - x[0]) + Bxx(x[0], 0, c[0], eps);

	for (int i = 1; i <= nLO; i++)
		for (int j = 1; j <= nLO; j++)
			H(i, j) = -Lxx_lam - Lxx_gam;

	for (int i = 1; i <= nLO; i++)
		H(i, i) += Lxx(x[i]) + Lxx(c[i] - x[i]) + Bxx(x[i], 0, c[i], eps);

	for (int i = 0; i <= nLO; i++)
		H(i, nLO + 1) = Derivative<0, 1>(kappa[i])(p, T);

	for (int i = 0; i <= nLO; i++)
		H(nLO + 1, i) = hGL[i](p, T);

	double cap = 0;
	for (int i = 0; i <= nLO + nHO + nG + nS; i++)
		cap += c[i] * Derivative<0, 1>(h[i])(p, T);
	for (int i = 0; i <= nLO; i++)
		cap -= x[i] * Derivative<0, 1>(hGL[i])(p, T);

	H(nLO + 1, nLO + 1) = -cap;
}

void PhaseSolver::computeFa() {
	double T = x[nLO + 1];

	double lam = 0;
	for (int i = 1; i <= nLO; i++)
		lam += x[i];
	for (int i = nLO + 1; i <= nLO + nHO; i++)
		lam += c[i];

	double gam = c[0] - x[0];
	for (int i = 1; i <= nLO; i++)
		gam += c[i] - x[i];

	for (int i = nLO + nHO + 1; i <= nLO + nHO + nG; i++)
		gam += c[i];

	double Lxx_gam = Lxx(gam);
	double Lxx_lam = Lxx(lam);

	/* wrt p */
	for (int i = 0; i <= nLO; i++)
		S(i, 0) = Derivative<1, 0>(kappa[i])(p, T);
	double comp = 0;
	for (int i = 0; i <= nLO + nHO + nG + nS; i++)
		comp += c[i] * Derivative<1, 0>(h[i])(p, T);
	for (int i = 0; i <= nLO; i++)
		comp -= x[i] * Derivative<1, 0>(hGL[i])(p, T);
	S(nLO + 1, 0) = -comp;

	/* wrt eta */
	for (int i = 0; i <= nLO; i++)
		S(i, 1) = 0;
	S(nLO + 1, 1) = 1;

	/* wrt c_0*/
	S(0, 2) = -Lxx(c[0] - x[0]) + Lxx_gam - eps / ((c[0] - x[0]) * (c[0] - x[0]));
	for (int i = 1; i <= nLO; i++)
		S(i, 2) = Lxx_gam;
	S(nLO + 1, 2) = -h[0](p, T);

	/* wrt c_j, j is LO*/
	for (int j = 1; j <= nLO; j++) {
		S(0, 2 + j) = Lxx_gam;
		for (int i = 1; i <= nLO; i++)
			S(i, 2 + j) = Lxx_gam;
		S(j, 2 + j) -= Lxx(c[j] - x[j]) + eps / ((c[j] - x[j]) * (c[j] - x[j]));
		S(nLO + 1, 2 + j) = -h[j](p, T);
	}

	/* wrt c_j, j is HO*/
	for (int j = nLO + 1; j <= nLO + nHO; j++) {
		S(0, 2 + j) = 0;
		for (int i = 1; i <= nLO; i++)
			S(i, 2 + j) = -Lxx_lam;
		S(nLO + 1, 2 + j) = -h[j](p, T);
	}

	/* wrt c_j, j is G*/
	for (int j = nLO + nHO + 1; j <= nLO + nHO + nG; j++) {
		for (int i = 0; i <= nLO; i++)
			S(i, 2 + j) = Lxx_gam;
		S(nLO + 1, 2 + j) = -h[j](p, T);
	}

	/* wrt c_j, j is S*/
	for (int j = nLO + nHO + nG + 1; j <= nLO + nHO + nG + nS; j++) {
		for (int i = 0; i <= nLO; i++)
			S(i, 2 + j) = 0;
		S(nLO + 1, 2 + j) = -h[j](p, T);
	}
}

bool PhaseSolver::checkF() {
	computeF();

	double dh = 1e-6;
	bool good = true;

	for (int i = 0; i <= nLO; i++) {
		x[i] += dh;
		double fp = computePhi();
		x[i] -= 2 * dh;
		double fm = computePhi();
		x[i] += dh;

		double num = 0.5 * (fp - fm) / dh;
		double ex = G[i];

		std::cerr << "G[" << i << "]: num = " << num << " vs ex = " << ex << ", [" << ex - num << "]" << std::endl;
		if (fabs(num - ex) > 1e-8 * (fabs(num) + fabs(ex))) {
			std::cerr << "^^^ bad" << std::endl;
			good = false;
		}
	}
	return good;
}

bool PhaseSolver::checkFx() {
	computeFx();

	double dh = 1e-6;
	bool good = true;

	for (int j = 0; j <= nLO; j++) {
		x[j] += dh;
		computeF();
		la::vector<double> gp(G);
		x[j] -= 2 * dh;
		computeF();
		la::vector<double> gm(G);
		x[j] += dh;

		for (int i = 0; i <= nLO + 1; i++) {
			double ex = H(i, j);
			double num = (0.5 / dh) * (gp[i] - gm[i]);
			std::cerr << "H[" << i << ", " << j << "]: num = " << num << " vs ex = " << ex << ", [" << ex - num << "]" << std::endl;
			if (fabs(num - ex) > 1e-8 * (fabs(num) + fabs(ex))) {
				std::cerr << "^^^ bad" << std::endl;
				good = false;
			}
		}
	}

	/* wrt T*/ {
		const int j = nLO + 1;
		const double dT = dh * x[j];
		x[j] += dT;
		computeF();
		la::vector<double> gp(G);
		x[j] -= 2 * dT;
		computeF();
		la::vector<double> gm(G);
		x[j] += dT;

		for (int i = 0; i <= nLO + 1; i++) {
			double ex = H(i, j);
			double num = (0.5 / dT) * (gp[i] - gm[i]);
			std::cerr << "H[" << i << ", " << j << "]: num = " << num << " vs ex = " << ex << ", [" << ex - num << "]" << std::endl;
			if (fabs(num - ex) > 1e-8 * (fabs(num) + fabs(ex))) {
				std::cerr << "^^^ bad" << std::endl;
				good = false;
			}
		}
	}
	return good;
}

bool PhaseSolver::checkFa() {
	computeFa();

	double dh = 1e-6;
	bool good = true;

	/* wrt p */ {
		double dp = dh * p;
		p += dp;
		computeF();
		la::vector<double> gp(G);
		p -= 2 * dp;
		computeF();
		la::vector<double> gm(G);
		p += dp;

		for (int i = 0; i <= nLO + 1; i++) {
			double ex = S(i, 0);
			double num = (0.5 / dp) * (gp[i] - gm[i]);
			std::cerr << "S[" << i << ", " << 0 << "]: num = " << num << " vs ex = " << ex << ", [" << ex - num << "]" << std::endl;
			if (fabs(num - ex) > 1e-8 * (fabs(num) + fabs(ex))) {
				std::cerr << "^^^ bad" << std::endl;
				good = false;
			}
		}
	}

	/* wrt eta */ {
		double de = dh * eta;
		eta += de;
		computeF();
		la::vector<double> gp(G);
		eta -= 2 * de;
		computeF();
		la::vector<double> gm(G);
		eta += de;

		for (int i = 0; i <= nLO + 1; i++) {
			double ex = S(i, 1);
			double num = (0.5 / de) * (gp[i] - gm[i]);
			std::cerr << "S[" << i << ", " << 1 << "]: num = " << num << " vs ex = " << ex << ", [" << ex - num << "]" << std::endl;
			if (fabs(num - ex) > 1e-8 * (fabs(num) + fabs(ex))) {
				std::cerr << "^^^ bad" << std::endl;
				good = false;
			}
		}
	}

	for (int j = 0; j < c.N(); j++) {
		c[j] += dh;
		computeF();
		la::vector<double> gp(G);
		c[j] -= 2 * dh;
		computeF();
		la::vector<double> gm(G);
		c[j] += dh;

		for (int i = 0; i <= nLO + 1; i++) {
			double ex = S(i, j + 2);
			double num = (0.5 / dh) * (gp[i] - gm[i]);
			std::cerr << "S[" << i << ", " << j + 2 << "]: num = " << num << " vs ex = " << ex << ", [" << ex - num << "]" << std::endl;
			if (fabs(num - ex) > 1e-7 * (fabs(num) + fabs(ex))) {
				std::cerr << "^^^ bad" << std::endl;
				good = false;
			}
		}
	}
	return good;
}
