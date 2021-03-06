#pragma once

#include <memory>
#include <string>
#include <tuple>

#include "image/Image.h"

using ImageLib::dImage;

template<typename T>
class VecT {
public:
  VecT(int cols, T val) : m_cols(cols) {
    m_array.resize(cols, val);
  }
  VecT(int cols, T* ptr) : m_cols(cols), m_ptr(ptr) {
  }

  T& operator[](int i) {
    assert(i >= 0);
    assert(i < m_cols);
    T *p = m_ptr ? m_ptr : m_array.data();
    return p[i];
  }

private:
  int m_cols = 0;
  T* m_ptr = nullptr;
  std::vector<T> m_array;
};

template<typename T>
class GridT {
public:
  GridT(int rows, int cols, T val) : m_rows(rows), m_cols(cols) {
    m_array.resize(rows*cols, val);
  }
  GridT(int rows, int cols, T* ptr) : m_rows(rows), m_cols(cols), m_ptr(ptr) {
  }

  VecT<T> operator[](int i) {
    assert(i >= 0);
    assert(i < m_rows);
    T *p = m_ptr ? m_ptr : m_array.data();
    return VecT<T>(m_cols, &p[i * m_cols]);
  }

private:
  int m_rows = 0;
  int m_cols = 0;
  T* m_ptr = nullptr;
  std::vector<T> m_array;
};

template<typename T>
class CubeT {
public:
  CubeT(int depths, int rows, int cols, T val) : m_depths(depths), m_rows(rows), m_cols(cols) {
    m_array.resize(depths*rows*cols, val);
  }
  CubeT(int depths, int rows, int cols, T* ptr) : m_depths(depths), m_rows(rows), m_cols(cols), m_ptr(ptr) {
  }

  GridT<T> operator[](int depth) {
    assert(depth >= 0);
    assert(depth < m_depths);
    T *p = m_ptr ? m_ptr : m_array.data();
    return GridT<T>(m_rows, m_cols, &p[depth*m_rows*m_cols]);
  }

private:
  int m_depths = 0;
  int m_rows = 0;
  int m_cols = 0;
  T* m_ptr = nullptr;
  std::vector<T> m_array;
};

typedef GridT<double> dGrid;
typedef CubeT<double> dCube;

typedef std::unique_ptr<ImageLib::dImage> ImgPtr;
typedef std::vector<double> VecDbl;
typedef std::vector<VecDbl> GridDbl;
typedef std::vector<ImgPtr> VecImg;
typedef std::vector<VecImg> GridImg;
typedef std::vector<GridImg> CubeImg;

typedef std::vector<int> VecInt;
typedef std::vector<VecInt> GridInt;
// Implements to algorithm in the paper by Modin and Karlsson (to be published).
class DiffeoFunctionMatching final {
public:
  //Parameters
  //----------
  // source      : array_like: Numpy array (float64) for the source image.
  // target      : array_like: Numpy array (float64) for the target image.
  //   target Must be of the same shape as `source`.
  // alpha       : float     : Parameter for ?
  // beta        : float     : Parameter for ?
  // sigma       : float     : Parameter for penalizing change of volume (divergence).
  // compute_phi : bool      : Whether to compute the forward phi mapping or not.
  //Returns
  //-------
  // null for invalid input or object
  // second in tuple is a message (usually descriptive error state)
  static std::tuple<std::unique_ptr<DiffeoFunctionMatching>, std::string> create(
    const dImage* source, dImage* target,
    double alpha, double beta, double sigma,
    bool compute_phi);

  // niter=300, epsilon=0.1
  // niter   - Number of iterations to take.
  // epsilon - The stepsize in the gradient descent method.
  void run(int niter, double epsilon);

private:
  DiffeoFunctionMatching(const dImage* source, dImage* target,
    double alpha, double beta, double sigma,
    bool compute_phi) :
    m_source(source), m_target(target), m_alpha(alpha), m_beta(beta), m_sigma(sigma),
    m_compute_phi(compute_phi)
  {
  }

  void setup();

  const dImage* m_source;
  dImage* m_target;
  double m_alpha;
  double m_beta;
  double m_sigma;
  bool m_compute_phi;

  int m_s = 0; // TODO: fix width, height

  VecDbl m_E;

  ImgPtr m_I0;
  ImgPtr m_I1;
  ImgPtr m_I;

  ImgPtr m_dIdx;
  ImgPtr m_dIdy;
  ImgPtr m_vx;
  ImgPtr m_vy;
  ImgPtr m_divv;

  GridDbl m_idx;
  GridDbl m_idy;
  GridDbl m_phiinvx;
  GridDbl m_phiinvy;
  GridDbl m_psiinvx;
  GridDbl m_psiinvy;

  GridDbl m_phix;
  GridDbl m_phiy;
  GridDbl m_psix;
  GridDbl m_psiy;

  GridDbl m_tmpx;
  GridDbl m_tmpy;

  GridImg m_g;
  GridImg m_h;

  ImgPtr m_hdet;
  ImgPtr m_dhaadx;
  ImgPtr m_dhbadx;
  ImgPtr m_dhabdx;
  ImgPtr m_dhbbdx;
  ImgPtr m_dhaady;
  ImgPtr m_dhbady;
  ImgPtr m_dhabdy;
  ImgPtr m_dhbbdy;

  ImgPtr m_yddy;
  ImgPtr m_yddx;
  ImgPtr m_xddy;
  ImgPtr m_xddx;

  CubeImg m_G;
  VecImg m_Jmap;

  GridDbl m_multipliers;
  GridDbl m_Linv;
};
