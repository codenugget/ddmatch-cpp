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


inline std::vector<std::complex<double>> idft(
	const std::vector<std::complex<double>>& x
	) {
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
