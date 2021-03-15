#include <gtest/gtest.h>

#include "core/MyArrays.h"

TEST(TGrid_Test, DefaultConstructor) {
  {
    TGrid<int> def;
    EXPECT_EQ(def.rows(), 0);
    EXPECT_EQ(def.cols(), 0);
    EXPECT_EQ(def.size(), 0);
    EXPECT_EQ(def.size_bytes(), 0);
    EXPECT_EQ(def.data(), nullptr);
  }
  {
    TGrid<float> def;
    EXPECT_EQ(def.rows(), 0);
    EXPECT_EQ(def.cols(), 0);
    EXPECT_EQ(def.size(), 0);
    EXPECT_EQ(def.size_bytes(), 0);
    EXPECT_EQ(def.data(), nullptr);
  }
  {
    TGrid<double> def;
    EXPECT_EQ(def.rows(), 0);
    EXPECT_EQ(def.cols(), 0);
    EXPECT_EQ(def.size(), 0);
    EXPECT_EQ(def.size_bytes(), 0);
    EXPECT_EQ(def.data(), nullptr);
  }
}

TEST(TGrid_Test, Constructor) {
  {
    TGrid<int> g1(16,8, 1);
    EXPECT_EQ(g1.rows(), 16);
    EXPECT_EQ(g1.cols(), 8);
    EXPECT_EQ(g1.size(), 16*8);
    EXPECT_EQ(g1.size_bytes(), 16*8*sizeof(int));
    ASSERT_NE(g1.data(), nullptr);
    EXPECT_EQ(g1[0][0], 1);
    {
      int* v = g1.data();
      for(int i = 0; i < g1.size(); ++i)
        ASSERT_EQ(v[i], 1);
    }
  }
  {
    TGrid<float> g2(8,16, 2.0f);
    EXPECT_EQ(g2.rows(), 8);
    EXPECT_EQ(g2.cols(), 16);
    EXPECT_EQ(g2.size(), 8*16);
    EXPECT_EQ(g2.size_bytes(), 8*16*sizeof(float));
    ASSERT_NE(g2.data(), nullptr);
    EXPECT_NEAR(g2[0][0], 2.0f, 1e-20f);
    {
      float* v = g2.data();
      for(int i = 0; i < g2.size(); ++i)
        ASSERT_NEAR(v[i], 2.0f, 1e-20f);
    }
  }
  {
    TGrid<double> g3(16,16, 3.0);
    EXPECT_EQ(g3.rows(), 16);
    EXPECT_EQ(g3.cols(), 16);
    EXPECT_EQ(g3.size(), 16*16);
    EXPECT_EQ(g3.size_bytes(), 16*16*sizeof(double));
    ASSERT_NE(g3.data(), nullptr);
    EXPECT_NEAR(g3[0][0], 3.0f, 1e-20);
    {
      double* v = g3.data();
      for(int i = 0; i < g3.size(); ++i)
        ASSERT_NEAR(v[i], 3.0f, 1e-20);
    }
  }
}

TEST(TGrid_Test, OutOfBounds) {
  TGrid<int> grid(2,2, 0);
  EXPECT_DEATH(grid[-1][0], "");
  EXPECT_DEATH(grid[0][-1], "");
  EXPECT_DEATH(grid[-1][-1], "");
  EXPECT_DEATH(grid[2][0], "");
  EXPECT_DEATH(grid[0][2], "");
  EXPECT_DEATH(grid[2][2], "");
}


TEST(TGrid_Test, values_like) {
  {
    TGrid<int> grid_orig(3,4, -1);
    TGrid<int> grid_test = values_like(grid_orig, 2);
    EXPECT_EQ(grid_test.rows(), 3);
    EXPECT_EQ(grid_test.cols(), 4);
    EXPECT_EQ(grid_test.size(), 12);
    EXPECT_EQ(grid_test.size_bytes(), 12 * sizeof(int));
    int* p = grid_test.data();
    ASSERT_NE(p, nullptr);
    for(int i = 0; i < grid_test.size(); ++i)
      EXPECT_EQ(p[i], 2) << " element " << i << " has unexpected value";
  }
}

TEST(TGrid_Test, zeros_like) {
  {
    TGrid<int> grid_orig(4,3, -1);
    TGrid<int> grid_test = zeros_like(grid_orig);
    EXPECT_EQ(grid_test.rows(), 4);
    EXPECT_EQ(grid_test.cols(), 3);
    EXPECT_EQ(grid_test.size(), 12);
    EXPECT_EQ(grid_test.size_bytes(), 12 * sizeof(int));
    int* p = grid_test.data();
    ASSERT_NE(p, nullptr);
    for(int i = 0; i < grid_test.size(); ++i)
      EXPECT_EQ(p[i], 0) << " element " << i << " has unexpected value";
  }
}

TEST(TGrid_Test, ones_like) {
  {
    TGrid<int> grid_orig(2,2, -1);
    TGrid<int> grid_test = ones_like(grid_orig);
    EXPECT_EQ(grid_test.rows(), 2);
    EXPECT_EQ(grid_test.cols(), 2);
    EXPECT_EQ(grid_test.size(), 4);
    EXPECT_EQ(grid_test.size_bytes(), 4 * sizeof(int));
    int* p = grid_test.data();
    ASSERT_NE(p, nullptr);
    for(int i = 0; i < grid_test.size(); ++i)
      EXPECT_EQ(p[i], 1) << " element " << i << " has unexpected value";
  }
}

TEST(TGrid_Test, copyto_TGrid_TGrid) {
  {
    TGrid<int> grid_a(2,2, -1);
    TGrid<int> grid_b(2,3, -2);
    EXPECT_FALSE(copyto(grid_a, grid_b));
    EXPECT_EQ(grid_a.rows(), 2);
    EXPECT_EQ(grid_a.cols(), 2);
    int* p = grid_a.data();
    ASSERT_NE(p, nullptr);
    for(int i = 0; i < grid_a.size(); ++i)
      EXPECT_EQ(p[i], -1) << " element " << i << " has unexpected value";
  }

  {
    TGrid<int> grid_a(2,2, -1);
    TGrid<int> grid_b(2,2, 4);
    EXPECT_TRUE(copyto(grid_a, grid_b));
    EXPECT_EQ(grid_a.rows(), 2);
    EXPECT_EQ(grid_a.cols(), 2);
    int* p = grid_a.data();
    ASSERT_NE(p, nullptr);
    for(int i = 0; i < grid_a.size(); ++i)
      EXPECT_EQ(p[i], 4) << " element " << i << " has unexpected value";
  }
}

TEST(TGrid_Test, copyto_TVecRef_stdvec) {
  {
    TGrid<int> grid_a(2,2, -1);
    std::vector<int> b{1,2,3};
    EXPECT_FALSE(copyto(grid_a[0], b));
    EXPECT_FALSE(copyto(grid_a[1], b));
    EXPECT_EQ(grid_a.rows(), 2);
    EXPECT_EQ(grid_a.cols(), 2);
    int* p = grid_a.data();
    ASSERT_NE(p, nullptr);
    for(int i = 0; i < grid_a.size(); ++i)
      EXPECT_EQ(p[i], -1) << " element " << i << " has unexpected value";
  }

  {
    TGrid<int> grid_a(2,2, -1);
    std::vector<int> b{1,2};
    std::vector<int> c{3,4};
    EXPECT_TRUE(copyto(grid_a[0], b));
    EXPECT_TRUE(copyto(grid_a[1], c));
    EXPECT_EQ(grid_a.rows(), 2);
    EXPECT_EQ(grid_a.cols(), 2);
    EXPECT_EQ(grid_a[0][0], 1);
    EXPECT_EQ(grid_a[0][1], 2);
    EXPECT_EQ(grid_a[1][0], 3);
    EXPECT_EQ(grid_a[1][1], 4);
  }
}

// elem_func
