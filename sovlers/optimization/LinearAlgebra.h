#ifndef __LINEARALGEBRA_H__
#define __LINEARALGEBRA_H__

#include <algorithm>
#include <functional>
#include <cassert>
#include <cmath>
#include <ostream>

#include "Util.h"

namespace la {

template <typename R>
class vector {
	R *v;
	int sz;
	template <class UnaryFunctor>
	vector<R> &transform(UnaryFunctor f) {
		std::transform(v, v + sz, v, f);
		return *this;
	}
	template <class BinaryFunctor>
	vector<R> &transform(vector<R> v2, BinaryFunctor f) {
		assert(sz == v2.sz);
		std::transform(v, v + sz, v2.v, v, f);
		return *this;
	}
	struct mul_by : std::unary_function<R,R> {
		R w;
		mul_by(R w) : w(w) {}
		R operator()(const R &v) {
			return w * v;
		}
	};
	struct set {
		R w;
		set(R w) : w(w) {}
		R operator()() { return w; }
	};
public:
	~vector() {
		delete[] v;
	}
	vector() : sz(0) {
		v = 0;
	}
	void resize(int _sz) {
		delete[] v;
		sz = _sz;
		v = new R[sz];
	}
	explicit vector(int sz) : sz(sz) {
		v = new R[sz];
	}
	vector(int sz, const R o) : sz(sz) {
		v = new R[sz];
		std::generate(v, v + sz, set(o));
	}
	vector(const vector<R> &o) {
		sz = o.sz;
		v = new R[sz];
		std::copy(o.v, o.v + sz, v);
	}
	vector<R> &operator=(const vector<R> &o) {
		assert(sz == o.sz);
		for (int i = 0; i < sz; i++)
			v[i] = o.v[i];
		return *this;
	}
	template <class UnaryFunctor>
	vector(const vector<R> &o, UnaryFunctor f) {
		sz = o.sz;
		v = new R[sz];
		std::transform(o.v, o.v + sz, v, f);
	}
	template <class BinaryFunctor>
	vector(const vector<R> &v1, const vector<R> &v2, BinaryFunctor f) {
		assert(v1.sz == v2.sz);
		sz = v1.sz;
		v = new R[sz];
		std::transform(v1.v, v1.v + sz, v2.v, v, f);
	}
	int N() const {
		return sz;
	}
	R &operator[](int n) {
		return v[n];
	}
	const R &operator[](int n) const {
		return v[n];
	}
	R &operator()(int n) {
		return v[n];
	}
	const R &operator()(int n) const {
		return v[n];
	}
	friend vector<R> operator +(const vector<R> &o1, const vector<R> &o) {
		return vector<R>(o1, o, std::plus<R>());
	}
	friend vector<R> operator -(const vector<R> &o1, const vector<R> &o) {
		return vector<R>(o1, o, std::minus<R>());
	}
	vector<R> operator -() const {
		return vector<R>(*this, std::negate<R>());
	}
	vector<R> operator *(const vector<R> &o) const {
		return vector<R>(*this, o, std::multiplies<R>());
	}
	vector<R> operator /(const vector<R> &o) const {
		return vector<R>(*this, o, std::divides<R>());
	}
	vector<R> &operator +=(const vector<R> &o) {
		return this->transform(o, std::plus<R>());
	}
	vector<R> &operator -=(const vector<R> &o) {
		return this->transform(o, std::minus<R>());
	}
	vector<R> &operator *=(const vector<R> &o) {
		return this->transform(o, std::multiplies<R>());
	}
	vector<R> &operator /=(const vector<R> &o) {
		return this->transform(o, std::divides<R>());
	}
	friend vector<R> operator *(const R &o, const vector<R> &v) {
		return vector(v, mul_by(o));
	}
	vector<R> &operator *=(const R &o) {
		return this->transform(mul_by(o));
	}
	R dot(const vector<R> &o) const {
		R sum = 0;
		for (int i = 0; i < sz; i++)
			sum += v[i] * o[i];
		return sum;
	}
	static R dot(const vector<R> &o1, const vector &o2) {
		return o1.dot(o2);
	}
	R norm() const {
		return std::sqrt(this->dot(*this));
	}
	friend std::ostream &operator<<(std::ostream &o, const vector<R> &vec) {
		o << "[";
		for (int i = 0; i < vec.N() - 1; i++)
			o << util::format("% 8.4e ", vec(i));
		o << util::format("% 8.4e]", vec(vec.N() - 1));
		return o;
	}
};

template <typename R> class matrix;

template <typename R>
struct qrinfo {
	int k;
	vector<int> pvt;
	vector<R> rdiag;
	qrinfo(const matrix<R> &a) : k(std::min(a.N(), a.M())), pvt(k), rdiag(k) { }
};

template <typename R>
class matrix {
	R *v;
	int m, n;
	template <class UnaryFunctor>
	matrix<R> &transform(UnaryFunctor f) {
		std::transform(v, v + n * m, v, f);
		return *this;
	}
	template <class BinaryFunctor>
	matrix<R> &transform(matrix<R> v2, BinaryFunctor f) {
		assert(n == v2.n && m == v2.m);
		std::transform(v, v + n * m, v2.v, v, f);
		return *this;
	}
	struct mul_by : std::unary_function<R,R> {
		R w;
		mul_by(R w) : w(w) {}
		R operator()(const R &v) {
			return w * v;
		}
	};
	struct set {
		R w;
		set(R w) : w(w) {}
		R operator()() { return w; }
	};
public:
	~matrix() {
		delete[] v;
	}
	matrix() : m(0), n(0) {
		v = 0;
	}
	matrix(int m, int n) : m(m), n(n) {
		v = new R[n * m];
	}
	matrix(const matrix<R> &o) {
		n = o.n;
		m = o.m;
		v = new R[n * m];
		std::copy(o.v, o.v + n * m, v);
	}
	int N() const {
		return n;
	}
	int M() const {
		return m;
	}
	void resize(int _m, int _n) {
		delete[] v;
		m = _m;
		n = _n;
		v = new R[n * m];
	}
	matrix<R> &operator=(const matrix<R> &o) {
		assert(n == o.n && m == o.m);
		for (int i = 0; i < m * n; i++)
			v[i] = o.v[i];
		return *this;
	}
	template <class UnaryFunctor>
	matrix(const matrix<R> &m1, UnaryFunctor f) {
		n = m1.n;
		m = m1.m;
		v = new R[n * m];
		std::transform(m1.v, m1.v + n * m, v, f);
	}
	template <class BinaryFunctor>
	matrix(const matrix<R> &m1, const matrix<R> &m2, BinaryFunctor f) {
		assert(m1.n == m2.n && m1.m == m2.m);
		n = m1.n;
		m = m1.m;
		v = new R[n * m];
		std::transform(m1.v, m1.v + n * m, m2.v, v, f);
	}
	void setIdentitiy() {
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				(*this)(i, j) = static_cast<R>((i == j) ? 1 : 0);
	}
	R &operator()(int i, int j) {
		return v[i * n + j];
	}
	const R &operator()(int i, int j) const {
		return v[i * n + j];
	}
	matrix<R> operator+(const matrix<R> &o) const {
		return matrix(*this, o, std::plus<R>());
	}
	matrix<R> &operator+=(const matrix<R> &o) {
		return transform(o, std::plus<R>());
	}
	matrix<R> operator-(const matrix<R> &o) const {
		return matrix(*this, o, std::minus<R>());
	}
	matrix<R> &operator-=(const matrix<R> &o) {
		return transform(o, std::minus<R>());
	}
	matrix<R> operator*(const R &o) const {
		return matrix(*this, mul_by(o));
	}
	friend matrix<R> operator*(const R &m, const matrix<R> &o) {
		return matrix(m, mul_by(o));
	}
	matrix<R> &negate() {
		return transform(std::negate<R>());
	}
	void lu_factorize() {
		assert(n == m);
		R p;
		matrix<R> &a = *this;
		for (int i = 0; i < m; i++) {
			a(i, i) = p = static_cast<R>(1) / a(i, i);
			for (int k = i+1; k < m; k++) {
				a(k, i) *= p;
				for (int j = i+1; j < n; j++)
					a(k, j) -= a(k, i) * a(i, j);
			}
		}
	}
	void lu_solve(vector<R> &b) const {
		assert(m == n && n == b.N());
		const matrix<R> &a = *this;
		for (int i = 1; i < m; i++)
			for (int j = 0; j < i; j++)
				b[i] -= a(i, j) * b[j];
		for (int i = m - 1; i >= 0; i--) {
			for (int j = i + 1; j < n; j++)
				b[i] -= a(i, j) * b[j];
			b[i] *= a(i, i);
		}
	}
	void lu_solve(matrix<R> &b) const {
		assert(m == n && n == b.M());
		const matrix<R> &a = *this;
		for (int k = 0; k < b.N(); k++) {
			for (int i = 1; i < m; i++)
				for (int j = 0; j < i; j++)
					b(i, k) -= a(i, j) * b(j, k);
			for (int i = m - 1; i >= 0; i--) {
				for (int j = i + 1; j < n; j++)
					b(i, k) -= a(i, j) * b(j, k);
				b(i, k) *= a(i, i);
			}
		}
	}
	void qr_factorize(qrinfo<R> &info) {
		int k = std::min(n, m);
		matrix<R> &a = *this;
		assert(info.pvt.N() == k && info.rdiag.N() == k);
		for (int j = 0; j < k; j++) {
			R alpha_max = 0;
			int q = j;
			for (int p = j; p < n; p++) {
				R alpha = 0;
				for (int i = j; i < m; i++)
					alpha += a(i, p) * a(i, p);
				if (alpha > alpha_max) {
					q = p;
					alpha_max = alpha;
				}
			}
			info.pvt[j] = q;
			if (q != j)
				for (int i = 0; i < m; i++)
					std::swap(a(i, j), a(i, q));
			alpha_max = std::sqrt(alpha_max);
			if (a(j, j) < 0)
				alpha_max = -alpha_max;
			info.rdiag[j] = alpha_max;
			a(j, j) -= alpha_max;
			R beta = 0;
			for (int i = j; i < m; i++)
				beta += a(i, j) * a(i, j);
			if (beta > 1e-20) {
				beta = 2 / beta;
				for (int p = j + 1; p < n; p++) {
					R gamma = 0;
					for (int i = j; i < n; i++)
						gamma += a(i, j) * a(i, p);
					gamma *= beta;

					for (int i = j; i < n; i++)
						a(i, p) -= gamma * a(i, j);
				}
			}
		}
	}
	void qr_factorize_noinfo() {
		int k = std::min(n, m);
		matrix<R> &a = *this;
		for (int j = 0; j < k; j++) {
			R alpha_max = 0;
			int q = j;
			for (int p = j; p < n; p++) {
				R alpha = 0;
				for (int i = j; i < m; i++)
					alpha += a(i, p) * a(i, p);
				if (alpha > alpha_max) {
					q = p;
					alpha_max = alpha;
				}
			}
			if (q != j)
				for (int i = 0; i < m; i++)
					std::swap(a(i, j), a(i, q));
			alpha_max = std::sqrt(alpha_max);
			if (a(j, j) < 0)
				alpha_max = -alpha_max;
			a(j, j) -= alpha_max;
			R beta = 0;
			for (int i = j; i < m; i++)
				beta += a(i, j) * a(i, j);
			if (beta > 1e-20) {
				beta = 2 / beta;
				for (int p = j + 1; p < n; p++) {
					R gamma = 0;
					for (int i = j; i < n; i++)
						gamma += a(i, j) * a(i, p);
					gamma *= beta;

					for (int i = j; i < n; i++)
						a(i, p) -= gamma * a(i, j);
				}
			}
		}
	}
	vector<R> annulator() {
		assert(n == m);
		int k = n;
		matrix<R> &a = *this;
		qr_factorize_noinfo();
		vector<R> ans(m, 0);
		vector<R> u(m, 0);
		ans(m - 1) = 1;
		for (int j = k - 1; j >= 0; j--) {
			for (int i = j; i < m; i++)
				u[i] = a(i, j);
			R beta = u.dot(u);
			if (beta > 1e-20) {
				beta = 2 / beta;
				ans -= (beta * u.dot(ans)) * u;
			}
		}
		return ans;
	}
	friend std::ostream &operator<<(std::ostream &o, const matrix<R> &A) {
		for (int i = 0; i < A.M(); i++) {
			o << "[";
			for (int j = 0; j < A.N() - 1; j++)
				o << util::format("% 8.4e ", A(i, j));
			o << util::format("% 8.4e]", A(i, A.N() - 1));
			if (i < A.M() - 1)
				o << std::endl;
		}
		return o;
	}
};

}; /* namespace la */

#endif
