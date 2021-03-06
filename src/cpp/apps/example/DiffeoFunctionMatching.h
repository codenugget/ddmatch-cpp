#pragma once

#include <memory>
#include <string>
#include <tuple>

#include "image/Image.h"

using ImageLib::dImage;

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
