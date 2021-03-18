#include <gtest/gtest.h>

#include "core/MyFft.h"


TEST(fft_Test, dft_idft_real) {
	const double cTolerance = 1e-10;
	{
		std::vector<std::complex<double>> orig_values{1,4,3,2};
		std::vector<std::complex<double>> freq_vals = dft(orig_values);
		ASSERT_EQ(freq_vals.size(), 4);
		EXPECT_NEAR(freq_vals[0].real(), 10, cTolerance);
		EXPECT_NEAR(freq_vals[0].imag(), 0, cTolerance);

		EXPECT_NEAR(freq_vals[1].real(), -2, cTolerance);
		EXPECT_NEAR(freq_vals[1].imag(), -2, cTolerance);

		EXPECT_NEAR(freq_vals[2].real(), -2, cTolerance);
		EXPECT_NEAR(freq_vals[2].imag(),  0, cTolerance);

		EXPECT_NEAR(freq_vals[3].real(), -2, cTolerance);
		EXPECT_NEAR(freq_vals[3].imag(),  2, cTolerance);

		std::vector<std::complex<double>> restored_values = idft(freq_vals);
		ASSERT_EQ(restored_values.size(), 4);

		EXPECT_NEAR(restored_values[0].real(), 1, cTolerance);
		EXPECT_NEAR(restored_values[0].imag(), 0, cTolerance);

		EXPECT_NEAR(restored_values[1].real(), 4, cTolerance);
		EXPECT_NEAR(restored_values[1].imag(), 0, cTolerance);

		EXPECT_NEAR(restored_values[2].real(), 3, cTolerance);
		EXPECT_NEAR(restored_values[2].imag(), 0, cTolerance);

		EXPECT_NEAR(restored_values[3].real(), 2, cTolerance);
		EXPECT_NEAR(restored_values[3].imag(), 0, cTolerance);
	}
}

TEST(fft_Test, dft_idft_complex) {
	const double cTolerance = 1e-10;
	{
		std::vector<std::complex<double>> orig_values{{1,-1},{4,4},{-3,3},{-2,-2}};
		std::vector<std::complex<double>> freq_vals = dft(orig_values);
		ASSERT_EQ(freq_vals.size(), 4);
		EXPECT_NEAR(freq_vals[0].real(), 0, cTolerance);
		EXPECT_NEAR(freq_vals[0].imag(), 4, cTolerance);

		EXPECT_NEAR(freq_vals[1].real(), 10, cTolerance);
		EXPECT_NEAR(freq_vals[1].imag(), -10, cTolerance);

		EXPECT_NEAR(freq_vals[2].real(), -4, cTolerance);
		EXPECT_NEAR(freq_vals[2].imag(),  0, cTolerance);

		EXPECT_NEAR(freq_vals[3].real(), -2, cTolerance);
		EXPECT_NEAR(freq_vals[3].imag(),  2, cTolerance);

		std::vector<std::complex<double>> restored_values = idft(freq_vals);
		ASSERT_EQ(restored_values.size(), 4);

		EXPECT_NEAR(restored_values[0].real(), 1, cTolerance);
		EXPECT_NEAR(restored_values[0].imag(),-1, cTolerance);

		EXPECT_NEAR(restored_values[1].real(), 4, cTolerance);
		EXPECT_NEAR(restored_values[1].imag(), 4, cTolerance);

		EXPECT_NEAR(restored_values[2].real(),-3, cTolerance);
		EXPECT_NEAR(restored_values[2].imag(), 3, cTolerance);

		EXPECT_NEAR(restored_values[3].real(),-2, cTolerance);
		EXPECT_NEAR(restored_values[3].imag(),-2, cTolerance);
	}
}
