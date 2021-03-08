#pragma once

#include <memory>
#include <string>
#include <tuple>

#include "core/MyArrays.h"
//#include "image/Image.h"

//using ImageLib::dImage;

//typedef std::unique_ptr<ImageLib::dImage> ImgPtr;
typedef std::vector<double> VecDbl;
//typedef std::vector<VecDbl> GridDbl;
typedef std::vector<dGrid> VecImg;
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
    const dGrid* source, dGrid* target,
    double alpha, double beta, double sigma,
    bool compute_phi);

  // niter=300, epsilon=0.1
  // niter   - Number of iterations to take.
  // epsilon - The stepsize in the gradient descent method.
  void run(int niter, double epsilon);

private:
  DiffeoFunctionMatching(const dGrid* source, dGrid* target,
    double alpha, double beta, double sigma,
    bool compute_phi) :
    m_source(source), m_target(target), m_alpha(alpha), m_beta(beta), m_sigma(sigma),
    m_compute_phi(compute_phi)
  {
  }

  void setup();

  const dGrid* m_source;
  dGrid* m_target;
  double m_alpha;
  double m_beta;
  double m_sigma;
  bool m_compute_phi;

  int m_s = 0; // TODO: fix width, height

  VecDbl m_E;

  dGrid m_I0;
  dGrid m_I1;
  dGrid m_I;

  dGrid m_dIdx;
  dGrid m_dIdy;
  dGrid m_vx;
  dGrid m_vy;
  dGrid m_divv;

  dGrid m_idx;
  dGrid m_idy;
  dGrid m_phiinvx;
  dGrid m_phiinvy;
  dGrid m_psiinvx;
  dGrid m_psiinvy;

  dGrid m_phix;
  dGrid m_phiy;
  dGrid m_psix;
  dGrid m_psiy;

  dGrid m_tmpx;
  dGrid m_tmpy;

  GridImg m_g;
  GridImg m_h;

  dGrid m_hdet;
  dGrid m_dhaadx;
  dGrid m_dhbadx;
  dGrid m_dhabdx;
  dGrid m_dhbbdx;
  dGrid m_dhaady;
  dGrid m_dhbady;
  dGrid m_dhabdy;
  dGrid m_dhbbdy;

  dGrid m_yddy;
  dGrid m_yddx;
  dGrid m_xddy;
  dGrid m_xddx;

  CubeImg m_G;
  VecImg m_Jmap;

  dGrid m_multipliers;
  dGrid m_Linv;
};
