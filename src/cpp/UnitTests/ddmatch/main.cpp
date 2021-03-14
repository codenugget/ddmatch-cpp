#include <gtest/gtest.h>

#include "ddmatch/Diffeo_functions.h"


TEST(Diffeo_functionsTest, boundary_conditions1) {
  const auto [i0, i1, sh0, sh1, frac] = boundary_conditions1(-1e-16, 64);
  EXPECT_EQ(i0, 63);
  EXPECT_EQ(i1, 0);
  EXPECT_NEAR(sh0, -64.0, 0.01);
  EXPECT_NEAR(sh1, 0.0, 0.01);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
