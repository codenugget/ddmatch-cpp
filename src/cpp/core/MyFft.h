#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>

#include "MyArrays.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384
#endif

inline std::vector<std::complex<double>> DFT(const std::vector<std::complex<double>>& x) {
	int N = static_cast<int>(x.size());
	TGrid<std::complex<double>> grid(N,N,{0.0,0.0});
	for(int r = 0; r < N; ++r) {
		for(int c = 0; c < N; ++c) {
			double Im = -2.0 * M_PI * static_cast<double>(r * c) / static_cast<double>(N);
			grid[r][c] = std::exp(std::complex(0.0, Im));
		}
	}
	std::vector<std::complex<double>> ret(N, {0.0,0.0});
	for(int i = 0; i < N; ++i) {
		for(int j = 0; j < N; ++j) {
			ret[i] += grid[i][j] * x[j];
		}
	}
	return ret;
}


inline std::vector<std::complex<double>> IDFT(const std::vector<std::complex<double>>& x) {
	int N = static_cast<int>(x.size());
	TGrid<std::complex<double>> grid(N,N,{0.0,0.0});
	for(int r = 0; r < N; ++r) {
		for(int c = 0; c < N; ++c) {
			double Im = 2.0 * M_PI * static_cast<double>(r * c) / static_cast<double>(N);
			grid[r][c] = std::exp(std::complex(0.0, Im));
		}
	}
	std::vector<std::complex<double>> ret(N, {0.0, 0.0});
	for(int n = 0; n < N; ++n) {
		for(int k=0; k < N; ++k)
			ret[n] += x[k] * grid[n][k];
		ret[n] /= static_cast<double>(N);
	}
	return ret;
}

// Algorithm S2
inline std::vector<std::complex<double>> ZeroPad(const std::vector<std::complex<double>>& x, int n) {
	int m = static_cast<int>(x.size());
	assert(m <= n);
	std::vector<std::complex<double>> x_hat = x;
	x_hat.resize(n, {0.0, 0.0});
	return x_hat;
}

// Algorithm S5: Runs in O(n log n) time. NOTE: only works when length(x) is a power of two
inline std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>>& x) {
	int n = static_cast<int>(x.size());
	if (n <= 1)
		return x;
	int nodd = n / 2 + ((n&1) ? 1 : 0);
	std::vector<std::complex<double>> xe, xo;
	xe.reserve(n / 2);
	xo.reserve(nodd);
	bool even = true;
	for(auto e : x) {
		if (even)
			xe.push_back(e);
		else
			xo.push_back(e);
		even = !even;
	}

	auto y_p = FFT(xe);
	auto y_pp = FFT(xo);
	std::vector<std::complex<double>> y(n);
	for(int k = 0; k < n/2; ++k) {
		std::complex w = std::exp(std::complex(0.0, -2*M_PI*double(k)/double(n)));
		y[k]       = y_p[k] + w * y_pp[k];
		y[k+(n/2)] = y_p[k] - w * y_pp[k];
	}
	return y;
}

// Algorithm S6: Runs in O(n log n) time.
inline std::vector<std::complex<double>> IFFT(const std::vector<std::complex<double>>& y) {
	int n = static_cast<int>(y.size());
	if (n <= 1)
		return y;
	int nodd = n / 2 + ((n&1) ? 1 : 0);
	std::vector<std::complex<double>> ye, yo;
	ye.reserve(n / 2);
	yo.reserve(nodd);
	bool even = true;
	for(auto e : y) {
		if (even)
			ye.push_back(e);
		else
			yo.push_back(e);
		even = !even;
	}
	auto x_p = IFFT(ye);
	auto x_pp = IFFT(yo);
	std::vector<std::complex<double>> x(n);
	for(int k = 0; k < n/2; ++k) {
		std::complex w = std::exp(std::complex(0.0, 2*M_PI*double(k)/double(n)));
		x[k]       = (x_p[k] + w * x_pp[k])/2.0;
		x[k+(n/2)] = (x_p[k] - w * x_pp[k])/2.0;
	}
	return x;
}

// Algorithm S4: Multiply a circulant matrix by a vector
inline std::vector<std::complex<double>> CirculantMultiply(
	const std::vector<std::complex<double>>& c,
	const std::vector<std::complex<double>>& x) {
	// Compute the product y = Gx of a circulant matrix G
	// and a vector x, where G is generated by its first column
	// c = (c[0], c[1], . . . , c[n−1]).
	int n = static_cast<int>(c.size());
	assert(static_cast<int>(x.size()) == n);

	auto C = FFT(c);
	auto X = FFT(x);
	std::vector<std::complex<double>> Y(n);
	for (int k = 0; k < n; ++k)
		Y[k] = C[k] * X[k];
	auto y = IFFT(Y);
	return y;
}

// Algorithm S1: Multiply a Toeplitz matrix by a vector using circulant embedding. Runs in O(n log n) time.
inline std::vector<std::complex<double>> ToeplitzMultiplyE(
	const std::vector<std::complex<double>>& r,
	const std::vector<std::complex<double>>& c,
	const std::vector<std::complex<double>>& x) {
	// Compute the product y = Tx of a Toeplitz matrix T
	// and a vector x, where T is specified by its first row
	// r = (r[0], r[1], r[2], . . . , r[N−1]) and its first column
	// c = (c[0], c[1], c[2], . . . , c[M−1]), where r[0] = c[0].
	int N = static_cast<int>(r.size());
	int M = static_cast<int>(c.size());
	assert(r[0] == c[0]);
	assert(static_cast<int>(x.size()) == N);

	// n=2^ceil(log2(M+N-1));
	int raise = static_cast<int>(ceil(log2(M + N - 1)));
	int n = 1 << raise;
	// Form an array c_hat by concatenating c, n − (M+N−1)
	// zeros, and the reverse of the last N−1 elements of r.
	//c_hat = ZeroArray(n);
	std::vector<std::complex<double>> c_hat(n, {0.0, 0.0});
	for (int k = 0; k < M; ++k)
		c_hat[k] = c[k];
	for (int k = 1; k < N; ++k)
		c_hat[n-k] = r[k];
	// c_hat= (c[0], c[1],..., c[M−1], 0,..., 0,r[N−1],...,r[2],r[1]);
	auto x_hat = ZeroPad(x, n);										// call Algorithm S2
	auto y_hat = CirculantMultiply(c_hat, x_hat); // call Algorithm S4

 	// The result is the first M elements of y_hat.
	std::vector<std::complex<double>> y(M);
	for(int k = 0; k < M; ++k)
		y[k] = y_hat[k];

	//y_hat.resize(M);
	//return y_hat;
	return y;
}

// Algorithm S7. Multiply a skew-circulant matrix by a vector. Runs in O(n log n) time.
std::vector<std::complex<double>> SkewCirculantMultiply(
	const std::vector<std::complex<double>>& c,
	const std::vector<std::complex<double>>& x) {
	int n = static_cast<int>(c.size());
	assert(static_cast<int>(x.size() == n));
	std::vector<std::complex<double>> c_hat(n);
	std::vector<std::complex<double>> x_hat(n);
	for(int k = 0; k < n; ++k) {
		std::complex<double> H = std::exp(std::complex<double>(0.0, -M_PI * double(k) / double(n)));
		c_hat[k] = c[k] * H; // c_hat = H * c
		x_hat[k] = x[k] * H; // x_hat = H * x
	}
	auto y = CirculantMultiply(c_hat, x_hat); // call algorithm S4
	for (int k = 0; k < n; ++k)
		y[k] = y[k] * std::exp(std::complex<double>(0.0, M_PI * double(k) / double(n)));
	return y;
}

// Algorithm S3: Multiply a Toeplitz matrix by a vector using Pustylnikov's decomposition. Runs in O(n log n) time.
inline std::vector<std::complex<double>> ToeplitzMultiplyP(
	const std::vector<std::complex<double>>& r,
	std::vector<std::complex<double>>& c, // NOTE: see if we can change to const reference later
	const std::vector<std::complex<double>>& x) {
	// Compute the product y = T x of a Toeplitz matrix T
	// and a vector x where T is specified by its first row
	// r = (r[0], r[1], r[2], . . . , r[N−1]) and its first column
	// c = (c[0], c[1], c[2], . . . , c[M−1]), where r[0] = c[0].
	int N = static_cast<int>(r.size());
	int M = static_cast<int>(c.size());
	assert(N > 0);
	assert(M > 0);
	assert(r[0] == c[0]);
	assert(static_cast<int>(x.size()) == N);
	int raise = static_cast<int>(ceil(log2(std::max(M, N))));
	int n = 1 << raise;
	if (N != n) // NOTE: when n differs it is always larger than N due to max statement above
		c = ZeroPad(c, n);
	std::vector<std::complex<double>> c_p(n);
	std::vector<std::complex<double>> c_pp(n);
	c_p[0] = 0.5 * c[0];
	c_pp[0] = 0.5 * c[0];
	for(int k = 1; k < n; ++k) {
		c_p[k]  = 0.5 * (c[k] + r[n-k]);
		c_pp[k] = 0.5 * (c[k] - r[n-k]);
	}
	auto y_p  = CirculantMultiply(c_p,  x);
	auto y_pp = SkewCirculantMultiply(c_pp, x);
	std::vector<std::complex<double>> y(M);
	for(int k = 0; k < M; ++k)
		y[k] = y_p[k] + y_pp[k];
	return y;
}

// Algorithm 1. CZT Algorithm. Runs in O(n log n) time.
// Te complex numbers A and W are parameters of the transform that defne the logarithmic spiral contour
// and the locations of the samples on it
inline std::vector<std::complex<double>> CZT(
	const std::vector<std::complex<double>>& x,
	const int M,
	const std::complex<double> W,
	const std::complex<double> A) {
	// A = diag(A^-0, A^-1, ..., A^-(N-1))
	int N = static_cast<int>(x.size());
	std::vector<std::complex<double>> X(N);
	std::vector<std::complex<double>> r(N);
	std::vector<std::complex<double>> c(M);
	for (int k = 0; k < N; ++k) {
		//X[k] = W^(k^2/2) * A^(-k) * x[k];
		X[k] = std::pow(W, k*k/2.0) * std::pow(A, -k) * x[k];
		//r[k] = W^(-k^2/2);
		r[k] = std::pow(W, -k*k/2.0); // NOTE: this expression should be possible to evaluated instead of storing it...
	}
	for (int k = 0; k < M; ++k) {
		//c[k] = W^(-k^2/2);
		c[k] = std::pow(W, -k*k/2.0); // NOTE: this expression should be possible to evaluated instead of storing it...
	}
	// After next line, Length(X) = M
	X = ToeplitzMultiplyE(r, c, X);
	for (int k = 0; k < M; ++k) {
		//X[k] = W^(k^2/2) * X[k];
		X[k] = std::pow(W, k*k/2.0) * X[k];
	}
	return X;
}

std::vector<std::complex<double>> ICZT(
	const std::vector<std::complex<double>>& X,
	const int N,
	const std::complex<double> W,
	const std::complex<double> A) {
	assert(static_cast<int>(X.size()) == N);
	int n = N;
	std::vector<std::complex<double>> x(n);
	for(int k = 0; k < n; ++k)
		x[k] = std::pow(W, -k*k/2.0) * X[k];
	// precompute the necessary polynomial products
	std::vector<std::complex<double>> p(n);
	p[0] = std::complex<double>(1.0, 0.0);
	for(int k = 1; k < n; ++k)
		p[k] = p[k-1] * (std::pow(W, k) - std::complex(1.0, 0.0));
	// copute the genarating vector u
	std::vector<std::complex<double>> u(n);
	for(int k = 0; k < n; ++k)
		u[k] = std::pow(-1, k) * std::pow(W, (2.0*k*k-(2.0*n-1)*k+n*(n-1))/2.0) / (p[n-k-1] * p[k]);
	std::vector<std::complex<double>> z(n, {0.0, 0.0});

	std::vector<std::complex<double>> u_hat(n);
	u_hat[0] = {0.0, 0.0};
	for (int k = 1; k < n; ++k)
		u_hat[n-k] = u[k];
	std::vector<std::complex<double>> u_tilde(n, {0.0, 0.0});
	u_tilde[0] = u[0];

	auto x_p = ToeplitzMultiplyE(u_hat, z, x);
	x_p = ToeplitzMultiplyE(z, u_hat, x_p);
	auto x_pp = ToeplitzMultiplyE(u, u_tilde, x);
	x_pp = ToeplitzMultiplyE(u_tilde, u, x_pp);
	for (int k = 0; k < n; ++k) {
		x[k] = (x_pp[k] - x_p[k]) / u[0]; // subtract and divide by u0
	}

	for (int k = 0; k < n; ++k) {
		x[k] = std::pow(A, k) * std::pow(W, -k*k/2.0) * x[k]; // multiply by A^-1 * Q^-1
	}
	return x;
}
