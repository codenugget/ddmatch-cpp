#pragma once

#include <cmath>
#include <complex>

#include "MyArrays.h"

inline std::vector<std::complex<double>> dft(const std::vector<std::complex<double>>& x) {
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


inline std::vector<std::complex<double>> idft(const std::vector<std::complex<double>>& x) {
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

// requires the size of x to be a power of 2
inline std::vector<std::complex<double>> fft_pow2(const std::vector<std::complex<double>>& x) {
	return {};
}

// 8-point FFT, with 0-padding
inline std::vector<std::complex<double>> fft_8pt0(const std::vector<std::complex<double>>& x) {
	return {};
}

/*
ZeroPad(x, n) {
	m = length(x);
	assert( m <= n );
	x_hat = EmptyArray(n);
	for (int k = 0; k < M; ++k) {
		x_hat[k] = x[k];
	}
	for (int k = M; k < N; ++k) {
		x_hat[k] = 0;
	}
	return x_hat;
}



ToeplitzMultiplyE(r,c,x){
	// compute the product y = Tx of a Toeplitz matrix T
	// and a vector x, where T is specified by its first row
	// r = (r[0], r[1], ..., r[N-1]) and its first column
	// c = (c[0], c[1], ..., c[M-1]) where r[0] == c[0]
	int N = static_cast<int>(r.size());
	int M = static_cast<int>(c.size());
	assert(r[0] == c[0]);
	assert(x.size() == N);
	n = 2^(log2(M+N-1));
	// Form an array c_hat by concatenating c, n - (M + N - 1)
	// zeros, and the reverseof the last N-1 elements of r
	c_hat = ZeroArray(n);
	for(int k = 0; k < M; ++k) {
		c_hat[k] = c[k];
	}
	for(int k = 1; k < N; ++k) {
		c_hat[n-k] = r[k];
	}
	// c_hat = (c[0], c[1], ..., c[M-1], 0, ..., 0, r[N-1], ..., r[2], r[1])

	x_hat = ZeroPad(x, n); // call Algorithm S2
	y_hat = CirculantMultiply(c_hat, x_hat); // call Algorithm S4
	// The result is the first M elements of y_hat
	y = EmptyArray(M);
	for(int k = 0; k < M; ++k) {
		y[k] = y_hat[k];
	}
	return y;
}

// chirp-z transform
inline std::vector<std::complex<double>> fft_czt(const std::vector<std::complex<double>>& x, int M) {
	int N = static_cast<int>(x.size()); //length(x);
	std::vector<std::complex<double>> X(N, {0.0, 0.0}); // X = emptyarray(N);
	std::vector<std::complex<double>> r(N, {0.0, 0.0}); // r = emptyarray(N);
	std::vector<std::complex<double>> c(N, {0.0, 0.0}); // c = emptyarray(N);
	for(int k = 0; k < N; ++i) {
		X[k] = W^(k^2/2) * A^(-k) * x[k];
		r[k] = W^(-k^2/2);
	}
	for (int k = 0; k < M; ++k) {
		c[k] = W^(-k^2/2);
	}
	// after this length(X) = M
	X = ToeplitzMultiplyE(r, c, X);
	for (int k = 0; k < M; ++k) {
		X[k] =W^(k^2/2) * X[k];
	}
	return X;
}
*/

inline std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& x) {
	// NOTES:
	//   How should non-power of 2 resolutions be handled?
	//   * evaluate dft directly
	//   * 8-point FFT, with 0-padding
	//   * chirp-z transform

	// placeholder implementation is to use the dft until the proper fft is fixed
	return dft(x);
}

// inverse chirp-z transform
inline std::vector<std::complex<double>> ifft_czt(const std::vector<std::complex<double>>& x) {
	return {};
}

inline std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& x) {
	// placeholder implementation is to use the idft until the proper fft is fixed
	return idft(x);
}
