#include <array>
#include <cmath>

#include "DiffeoFunctionMatching.h"
#include "Diffeo_functions.h"
#include "image/Image_funcs.h"

using ImageLib::dImage;

template<typename T>
std::vector<T> my_linspace(T start, T stop, int num, bool endpoint) {
  std::vector<T> ret(num, T(0));
  double step = 0;
  if (endpoint)
    step = double(stop-start) / (num-1);
  else
    step = double(stop-start) / num;

  for(int i = 0; i < num; ++i)
    ret[i] = T(start + i * step);
  return ret;
}

enum class Indexing {
  xy,
  ij
};
std::tuple<dGrid,dGrid> my_meshgrid(
  const std::vector<double>& x,
  const std::vector<double>& y,
  Indexing i = Indexing::xy) {
  switch(i) {
    case Indexing::xy: {
      int nx = (int)x.size();
      int ny = (int)y.size();
      dGrid xx(ny, nx, 0.0);
      dGrid yy(ny, nx, 0.0);

      for(int iy = 0; iy < ny; ++iy) {
        for(int ix = 0; ix < nx; ++ix) {
          xx[iy][ix] = x[ix];
          yy[iy][ix] = y[iy];
        }
      }
      return std::make_tuple(xx, yy);
    }
    case Indexing::ij: {
      // TODO: verify this case
      int nx = (int)x.size();
      int ny = (int)y.size();
      dGrid xx(ny, nx, 0.0);
      dGrid yy(ny, nx, 0.0);

      for(int iy = 0; iy < ny; ++iy) {
        for(int ix = 0; ix < nx; ++ix) {
          xx[iy][ix] = y[ix];
          yy[iy][ix] = x[iy];
        }
      }
      return std::make_tuple(xx, yy);
    }
    default: {
      // unknown / unsupported case
      return std::make_tuple(dGrid{}, dGrid{});
    }
  }
}

std::tuple<std::unique_ptr<DiffeoFunctionMatching>, std::string> DiffeoFunctionMatching::create(
    const dGrid& source, const dGrid& target,
    double alpha, double beta, double sigma,
    bool compute_phi) {
  // Check input
  if (!source.is_same_shape(target))
    return std::make_tuple(nullptr, "source and target are not of the same shape");

  if (sigma < 0)
    return std::make_tuple(nullptr, "Paramter sigma must be positive");

  if (source.rows() != source.cols())
    return std::make_tuple(nullptr, "Only square images allowed so far.");

  auto ret = std::unique_ptr<DiffeoFunctionMatching>(new DiffeoFunctionMatching(source, target, alpha, beta, sigma, compute_phi));
  ret->setup();
  return std::make_tuple(std::move(ret), "");
}

void DiffeoFunctionMatching::setup() {
  using ImageLib::zeros_like;
  using ImageLib::ones_like;
  dGrid& I0 = m_target;
  const dGrid& I1 = m_source;

  // Create optimized algorithm functions
  //self.image_compose = image_compose_2d;//generate_optimized_image_composition(I1)
  //self.diffeo_compose = diffeo_compose_2d;//generate_optimized_diffeo_composition(I1)
  //self.image_gradient = image_gradient_2d;//generate_optimized_image_gradient(I1)
  //self.diffeo_gradient_y = diffeo_gradient_y_2d;//generate_optimized_diffeo_gradient_y(I1)
  //self.diffeo_gradient_x = diffeo_gradient_x_2d;//generate_optimized_diffeo_gradient_x(I1)
  //self.evaluate = eval_diffeo_2d;//generate_optimized_diffeo_evaluation(I1)

  // Allocate and initialize variables
  m_s = I1.cols();
  //E = []
  m_E.resize(0);

  //self.I0 = np.zeros_like(I0)
  //np.copyto(self.I0,I0)
  m_I0 = I0;
  //self.I1 = np.zeros_like(I1)
  //np.copyto(self.I1,I1)
  m_I1 = I1;
  //self.I = np.zeros_like(I1)
  //np.copyto(self.I,I1)
  m_I = I1;

  //self.dIdx = np.zeros_like(I1)
  //self.dIdy = np.zeros_like(I1)
  //self.vx = np.zeros_like(I1)
  //self.vy = np.zeros_like(I1)
  //self.divv = np.zeros_like(I1)
  m_dIdx = zeros_like(I1);
  m_dIdy = zeros_like(I1);
  m_vx = zeros_like(I1);
  m_vy = zeros_like(I1);
  m_divv = zeros_like(I1);

  // Allocate and initialize the diffeos
  //x = np.linspace(0, self.s, self.s, endpoint=False)
  //[self.idx, self.idy] = np.meshgrid(x, x)
  auto x = my_linspace<double>(0, m_s, m_s, false);
  std::tie(m_idx, m_idy) = my_meshgrid(x,x);

  //self.phiinvx = self.idx.copy()
  //self.phiinvy = self.idy.copy()
  //self.psiinvx = self.idx.copy()
  //self.psiinvy = self.idy.copy()
  m_phiinvx = m_idx;
  m_phiinvy = m_idy;
  m_psiinvx = m_idx;
  m_psiinvy = m_idy;

  if (m_compute_phi) {
    //self.phix = self.idx.copy()
    //self.phiy = self.idy.copy() 
    //self.psix = self.idx.copy()
    //self.psiy = self.idy.copy() 
    m_phix = m_idx;
    m_phiy = m_idy;
    m_psix = m_idx;
    m_psiy = m_idy;
  }

  //self.tmpx = self.idx.copy()
  //self.tmpy = self.idy.copy()
  m_tmpx = m_idx;
  m_tmpy = m_idy;

  //# test case
  //#self.phiinvy += 5.e-8*self.phiinvy**2*(self.s-1-self.phiinvy)**2 + 5.e-8*self.phiinvx**2*(self.s-1-self.phiinvx)**2# compare with += 3.e-7*(...)
  //#self.phiinvx += 1.e-7*self.phiinvx**2*(self.s-1-self.phiinvx)**2

  // Allocate and initialize the metrics
  m_g.resize(2);
  //self.g = np.array([[np.ones_like(I1),np.zeros_like(I1)],[np.zeros_like(I1),np.ones_like(I1)]])
  m_g[0].emplace_back( ones_like(I1));
  m_g[0].emplace_back(zeros_like(I1));
  m_g[1].emplace_back(zeros_like(I1));
  m_g[1].emplace_back( ones_like(I1));

  m_h.resize(2);
  //self.h = np.array([[np.ones_like(I1),np.zeros_like(I1)],[np.zeros_like(I1),np.ones_like(I1)]])
  m_h[0].emplace_back(zeros_like(I1));
  m_h[0].emplace_back(zeros_like(I1));
  m_h[1].emplace_back(zeros_like(I1));
  m_h[1].emplace_back(zeros_like(I1));

  //self.hdet = np.zeros_like(I1)
  //self.dhaadx = np.zeros_like(I1)
  //self.dhbadx = np.zeros_like(I1)
  //self.dhabdx = np.zeros_like(I1)
  //self.dhbbdx = np.zeros_like(I1)
  //self.dhaady = np.zeros_like(I1)
  //self.dhbady = np.zeros_like(I1)
  //self.dhabdy = np.zeros_like(I1)
  //self.dhbbdy = np.zeros_like(I1)
  //self.yddy = np.zeros_like(I1)
  //self.yddx = np.zeros_like(I1)
  //self.xddy = np.zeros_like(I1)
  //self.xddx = np.zeros_like(I1)
  m_hdet   = zeros_like(I1);
  m_dhaadx = zeros_like(I1);
  m_dhbadx = zeros_like(I1);
  m_dhabdx = zeros_like(I1);
  m_dhbbdx = zeros_like(I1);
  m_dhaady = zeros_like(I1);
  m_dhbady = zeros_like(I1);
  m_dhabdy = zeros_like(I1);
  m_dhbbdy = zeros_like(I1);
  m_yddy   = zeros_like(I1);
  m_yddx   = zeros_like(I1);
  m_xddy   = zeros_like(I1);
  m_xddx   = zeros_like(I1);

  // self.G = np.zeros_like(np.array([self.g,self.g]))
  m_G.resize(2);
  m_G[0].resize(2);
  m_G[1].resize(2);
  m_G[0][0].emplace_back(zeros_like(I1));
  m_G[0][0].emplace_back(zeros_like(I1));
  m_G[0][1].emplace_back(zeros_like(I1));
  m_G[0][1].emplace_back(zeros_like(I1));
  m_G[1][0].emplace_back(zeros_like(I1));
  m_G[1][0].emplace_back(zeros_like(I1));
  m_G[1][1].emplace_back(zeros_like(I1));
  m_G[1][1].emplace_back(zeros_like(I1));

  //self.Jmap = np.zeros_like(np.array([I1,I1]))
  m_Jmap.emplace_back(zeros_like(I1));
  m_Jmap.emplace_back(zeros_like(I1));

  // Create wavenumber vectors
  //k = [np.hstack((np.arange(n//2),np.arange(-n//2,0))) for n in self.I0.shape]
  dGrid k(2, I1.rows(), 0.0);

  auto to_double = [](const VecInt& v) -> VecDbl {
    VecDbl ret;
    ret.reserve(v.size());
    for(auto i : v)
      ret.push_back(double(i));
    return ret;
  };

  //k.resize(2);
  auto v1 = to_double(my_linspace<int>(0, I1.cols(), I1.rows(), false));
  copyto(k[0], v1);
  auto v2 = to_double(my_linspace<int>(-I1.cols(), 0, I1.rows(), false));
  copyto(k[1], v2);

  /* not completed yet
  # Create wavenumber tensors
  K = np.meshgrid(*k, sparse=False, indexing='ij')
  */
  dGrid Kx, Ky;
  std::tie(Kx, Ky) = my_meshgrid(v1, v2, Indexing::ij);

  // Create Fourier multiplicator
  //self.multipliers = np.ones_like(K[0])
  //self.multipliers = self.multipliers*self.alpha
  //for Ki in K:
  //  Ki = Ki*self.beta
  //  self.multipliers = self.multipliers+Ki**2
  const auto mul_sq = [&](const double v) -> double {
    double m = v * m_beta;
    return m*m;
  };

  m_multipliers = values_like(Kx, m_alpha) +
    elem_func(Kx, mul_sq) + elem_func(Ky, mul_sq);
  //m_multipliers = values_like(Kx, m_alpha) + mul_pow(Kx, m_beta, 2) + mul_pow(Ky, m_beta, 2);

  //if self.alpha == 0:
  //  self.multipliers[0,0]=1.0#self.multipliers[(0 for _ in self.s)] = 1. # Avoid division by zero
  //  self.Linv = 1./self.multipliers
  //  self.multipliers[0,0]=0.
  //else:
  //  self.Linv = 1./self.multipliers
  const auto inv_f = [&](const double v) -> double {
    return 1.0 / v;
  };
  if (m_alpha == 0) {
    m_multipliers[0][0] = 1.0;
    m_Linv = elem_func(m_multipliers, inv_f);
    m_multipliers[0][0] = 0.0;
  }
  else {
    m_Linv = elem_func(m_multipliers, inv_f);
  }
}

/*bool copyto(GridDbl& dst, dImage* I) {
  int h = dst.size();
  if (h == 0)
    return false;
  int w = dst[0].size();
  if (!I or I->height() != h or I->width() != w)
    return false;

  double* v = I->data();
  for(int i = 0; i < h; ++i)
    memcpy(dst[i].data(), v, w * sizeof(double));
  return true;
}

double sum(const GridDbl& vals) {
  double s = 0;
  for(const auto& row : vals) {
    for(const auto v : row)
      s += v;
  }
  return s;
}

void set_zero(GridDbl& vals) {
  for(auto& row : vals) {
    for(auto& v : row)
      v = 0;
  }
}

  template<typename Pred>
  bool elem_set(GridDbl& res, const dImage* a, const dImage* b, Pred f) {
    int w = a->width();
    int h = a->height();
    for(int y = 0; y < h; ++y) {
      for(int x = 0; x < w; ++x) {
        res[y][x] = f(a->get(x,y,0), b->get(x,y,0));
      }
    }
    return true;
  }
  template<typename Pred>
  bool elem_add(GridDbl& res, const dImage* a, const dImage* b, Pred f) {
    int w = a->width();
    int h = a->height();
    for(int y = 0; y < h; ++y) {
      for(int x = 0; x < w; ++x) {
        res[y][x] += f(a->get(x,y,0), b->get(x,y,0));
      }
    }
    return true;
  }
*/

void DiffeoFunctionMatching::run(int niter, double epsilon) {
  // Carry out the matching process.
  // Implements to algorithm in the paper by Modin and Karlsson
  // kE = len(self.E)
  // self.E = np.hstack((self.E,np.zeros(niter)))
  int kE = (int) m_E.size();
  m_E.resize(kE+niter, 0);
  for(int k = 0; k < niter; ++k) {
    // OUTPUT
    //np.copyto(self.tmpx, self.I)
    //np.copyto(self.tmpy, self.I0)
    copyto(m_tmpx, m_I);
    copyto(m_tmpy, m_I0);

    //self.tmpx = self.tmpx-self.tmpy
    //self.tmpx **= 2
    //self.E[k+kE] = self.tmpx.sum()
    m_tmpx = m_tmpx - m_tmpy;
    m_tmpx = elem_func(m_tmpx, [](double v){ return v*v; });
    m_E[k+kE] = sum(m_tmpx);

    //np.copyto(self.tmpx, (self.h[0,0]-self.g[0,0])**2+(self.h[1,0]-self.g[1,0])**2+\
    //  (self.h[0,1]-self.g[0,1])**2+(self.h[1,1]-self.g[1,1])**2)
    const auto diff_sq = [](const double v1, const double v2){
      double v = (v1-v2);
      return v*v;
    };
    elem_set(m_tmpx, m_h[0][0], m_g[0][0], diff_sq);
    elem_add(m_tmpx, m_h[1][0], m_g[1][0], diff_sq);
    elem_add(m_tmpx, m_h[0][1], m_g[0][1], diff_sq);
    elem_add(m_tmpx, m_h[1][1], m_g[1][1], diff_sq);

    //self.E[k+kE] += self.sigma*self.tmpx.sum()
    m_E[k+kE] += m_sigma * sum(m_tmpx);

/*
    //self.image_compose(self.I1, self.phiinvx, self.phiinvy, self.I)
    image_compose_2d(m_I1.get(), m_phiinvx, m_phiinvy, m_I.get());
    //self.diffeo_gradient_y(self.phiinvy, self.yddx, self.yddy)
    //self.diffeo_gradient_x(self.phiinvx, self.xddx, self.xddy)
    diffeo_gradient_y_2d(m_phiinvy, m_yddx.get(), m_yddy.get());
    diffeo_gradient_x_2d(m_phiinvx, m_xddx.get(), m_xddy.get());
    */
  }
}

/*
      np.copyto(self.h[0,0], self.yddy*self.yddy+self.xddy*self.xddy)
      np.copyto(self.h[1,0], self.yddx*self.yddy+self.xddx*self.xddy)
      np.copyto(self.h[0,1], self.yddy*self.yddx+self.xddy*self.xddx)
      np.copyto(self.h[1,1], self.yddx*self.yddx+self.xddx*self.xddx)

      self.image_gradient(self.h[0,0], self.dhaadx, self.dhaady)
      self.image_gradient(self.h[0,1], self.dhabdx, self.dhabdy)
      self.image_gradient(self.h[1,0], self.dhbadx, self.dhbady)
      self.image_gradient(self.h[1,1], self.dhbbdx, self.dhbbdy)

      np.copyto(self.Jmap[0], -(self.h[0,0]-self.g[0,0])*self.dhaady -(self.h[0,1]-self.g[0,1])*self.dhabdy -(self.h[1,0]-self.g[1,0])*self.dhbady -(self.h[1,1]-self.g[1,1])*self.dhbbdy +\
        2*self.dhaady*self.h[0,0] + 2*self.dhabdx*self.h[0,0] + 2*self.dhbady*self.h[1,0] + 2*self.dhbbdx*self.h[1,0] +\
        2*(self.h[0,0]-self.g[0,0])*self.dhaady + 2*(self.h[1,0]-self.g[1,0])*self.dhbady + 2*(self.h[0,1]-self.g[0,1])*self.dhaadx + 2*(self.h[1,1]-self.g[1,1])*self.dhbadx)
      
      np.copyto(self.Jmap[1], -(self.h[0,0]-self.g[0,0])*self.dhaadx -(self.h[0,1]-self.g[0,1])*self.dhabdx -(self.h[1,0]-self.g[1,0])*self.dhbadx -(self.h[1,1]-self.g[1,1])*self.dhbbdx +\
        2*self.dhaady*self.h[0,1] + 2*self.dhabdx*self.h[0,1] + 2*self.dhbady*self.h[1,1] + 2*self.dhbbdx*self.h[1,1] +\
        2*(self.h[0,0]-self.g[0,0])*self.dhabdy + 2*(self.h[1,0]-self.g[1,0])*self.dhbbdy + 2*(self.h[0,1]-self.g[0,1])*self.dhabdx + 2*(self.h[1,1]-self.g[1,1])*self.dhbbdx)

      self.image_gradient(self.I, self.dIdx, self.dIdy)
      self.vx = -(self.I-self.I0)*self.dIdx + 2*self.sigma*self.Jmap[1]# axis: [1]
      self.vy = -(self.I-self.I0)*self.dIdy + 2*self.sigma*self.Jmap[0]# axis: [0]
      fftx = np.fft.fftn(self.vx)
      ffty = np.fft.fftn(self.vy)
      fftx *= self.Linv
      ffty *= self.Linv
      self.vx[:] = -np.fft.ifftn(fftx).real # vx[:]=smth will copy while vx=smth directs a pointer
      self.vy[:] = -np.fft.ifftn(ffty).real

      # STEP 4 (v = -grad E, so to compute the inverse we solve \psiinv' = -epsilon*v o \psiinv)
      np.copyto(self.tmpx, self.vx)
      self.tmpx *= epsilon
      np.copyto(self.psiinvx, self.idx)
      self.psiinvx -= self.tmpx
      if self.compute_phi: # Compute forward phi also (only for output purposes)
        np.copyto(self.psix, self.idx)
        self.psix += self.tmpx

      np.copyto(self.tmpy, self.vy)
      self.tmpy *= epsilon
      np.copyto(self.psiinvy, self.idy)
      self.psiinvy -= self.tmpy
      if self.compute_phi: # Compute forward phi also (only for output purposes)
        np.copyto(self.psiy, self.idy)
        self.psiy += self.tmpy

      self.diffeo_compose(self.phiinvx, self.phiinvy, self.psiinvx, self.psiinvy, \
                self.tmpx, self.tmpy) # Compute composition phi o psi = phi o (1-eps*v)
      np.copyto(self.phiinvx, self.tmpx)
      np.copyto(self.phiinvy, self.tmpy)
      if self.compute_phi: # Compute forward phi also (only for output purposes)
        self.diffeo_compose(self.psix, self.psiy, \
                  self.phix, self.phiy, \
                  self.tmpx, self.tmpy)
        np.copyto(self.phix, self.tmpx)
        np.copyto(self.phiy, self.tmpy)
*/


/*

class DiffeoFunctionMatching(object):
  """
  Implementation of the two component function matching algorithm.

  The computations are accelerated using the `numba` library.
  """

  def __init__(self, source, target, alpha=0.001, beta=0.03, sigma=0.5, compute_phi=True):
    """
    Initialize the matching process.

    Implements to algorithm in the paper by Modin and Karlsson (to be published).

    Parameters
    ----------
    source : array_like
      Numpy array (float64) for the source image.
    target : array_like
      Numpy array (float64) for the target image.
      Must be of the same shape as `source`.
    sigma : float
      Parameter for penalizing change of volume (divergence).
    compute_phi : bool
      Whether to compute the forward phi mapping or not.

    Returns
    -------
    None or Python generator

    See also
    --------
    N/A

    Examples
    --------
    N/A
    """ 
    
    self.source = source
    self.target = target
    self.compute_phi = compute_phi
    I0 = target
    I1 = source

    # Check input
    if (I0.shape != I1.shape):
      raise(TypeError('Source and target images must have the same shape.'))
    if (I0.dtype != I1.dtype):
      raise(TypeError('Source and target images must have the same dtype.'))
    if (sigma < 0):
      raise(ValueError('Paramter sigma must be positive.'))
    for d in I1.shape:
      if (d != I1.shape[0]):
        raise(NotImplementedError('Only square images allowed so far.'))
    if (len(I1.shape) != 2):
      raise(NotImplementedError('Only 2d images allowed so far.'))

    # Create optimized algorithm functions
    self.image_compose = generate_optimized_image_composition(I1)
    self.diffeo_compose = generate_optimized_diffeo_composition(I1)
    self.image_gradient = generate_optimized_image_gradient(I1)
    self.diffeo_gradient_y = generate_optimized_diffeo_gradient_y(I1)
    self.diffeo_gradient_x = generate_optimized_diffeo_gradient_x(I1)
    self.evaluate = generate_optimized_diffeo_evaluation(I1)

    # Allocate and initialize variables
    self.alpha = alpha
    self.beta = beta
    self.sigma = sigma
    self.s = I1.shape[0]
    self.E = []
    self.I0 = np.zeros_like(I0)
    np.copyto(self.I0,I0)
    self.I1 = np.zeros_like(I1)
    np.copyto(self.I1,I1)
    self.I = np.zeros_like(I1)
    np.copyto(self.I,I1)
    self.dIdx = np.zeros_like(I1)
    self.dIdy = np.zeros_like(I1)
    self.vx = np.zeros_like(I1)
    self.vy = np.zeros_like(I1)
    self.divv = np.zeros_like(I1)
        
    # Allocate and initialize the diffeos
    x = np.linspace(0, self.s, self.s, endpoint=False)
    [self.idx, self.idy] = np.meshgrid(x, x)
    self.phiinvx = self.idx.copy()
    self.phiinvy = self.idy.copy()
    self.psiinvx = self.idx.copy()
    self.psiinvy = self.idy.copy()
    if self.compute_phi:
      self.phix = self.idx.copy()
      self.phiy = self.idy.copy() 
      self.psix = self.idx.copy()
      self.psiy = self.idy.copy() 
    self.tmpx = self.idx.copy()
    self.tmpy = self.idy.copy()

    # test case
    #self.phiinvy += 5.e-8*self.phiinvy**2*(self.s-1-self.phiinvy)**2 + 5.e-8*self.phiinvx**2*(self.s-1-self.phiinvx)**2# compare with += 3.e-7*(...)
    #self.phiinvx += 1.e-7*self.phiinvx**2*(self.s-1-self.phiinvx)**2


    # Allocate and initialize the metrics
    self.g = np.array([[np.ones_like(I1),np.zeros_like(I1)],[np.zeros_like(I1),np.ones_like(I1)]])
    self.h = np.array([[np.ones_like(I1),np.zeros_like(I1)],[np.zeros_like(I1),np.ones_like(I1)]])
    self.hdet = np.zeros_like(I1)
    self.dhaadx = np.zeros_like(I1)
    self.dhbadx = np.zeros_like(I1)
    self.dhabdx = np.zeros_like(I1)
    self.dhbbdx = np.zeros_like(I1)
    self.dhaady = np.zeros_like(I1)
    self.dhbady = np.zeros_like(I1)
    self.dhabdy = np.zeros_like(I1)
    self.dhbbdy = np.zeros_like(I1)
    self.yddy = np.zeros_like(I1)
    self.yddx = np.zeros_like(I1)
    self.xddy = np.zeros_like(I1)
    self.xddx = np.zeros_like(I1)
    self.G = np.zeros_like(np.array([self.g,self.g]))
    self.Jmap = np.zeros_like(np.array([I1,I1]))

    # Create wavenumber vectors
    k = [np.hstack((np.arange(n//2),np.arange(-n//2,0))) for n in self.I0.shape]

    # Create wavenumber tensors
    K = np.meshgrid(*k, sparse=False, indexing='ij')

    # Create Fourier multiplicator
    self.multipliers = np.ones_like(K[0])
    self.multipliers = self.multipliers*self.alpha
    for Ki in K:
      Ki = Ki*self.beta
      self.multipliers = self.multipliers+Ki**2
    if self.alpha == 0:
      self.multipliers[0,0]=1.0#self.multipliers[(0 for _ in self.s)] = 1. # Avoid division by zero
      self.Linv = 1./self.multipliers
      self.multipliers[0,0]=0.
    else:
      self.Linv = 1./self.multipliers

    
  def run(self, niter=300, epsilon=0.1):
    """
    Carry out the matching process.

    Implements to algorithm in the paper by Modin and Karlsson (to appear).

    Parameters
    ----------
    niter : int
      Number of iterations to take.
    epsilon : float
      The stepsize in the gradient descent method.
    yielditer : bool
      If `True`, then a yield statement is executed at the start of
      each iterations. This is useful for example when animating 
      the warp in real-time.

    Returns
    -------
    None or Python generator

    See also
    --------
    N/A

    Examples
    --------
    N/A
    """   

    kE = len(self.E)
    self.E = np.hstack((self.E,np.zeros(niter)))

    for k in range(niter):
      
      # OUTPUT
      np.copyto(self.tmpx, self.I)
      np.copyto(self.tmpy, self.I0)
      self.tmpx = self.tmpx-self.tmpy
      self.tmpx **= 2
      self.E[k+kE] = self.tmpx.sum()
      np.copyto(self.tmpx, (self.h[0,0]-self.g[0,0])**2+(self.h[1,0]-self.g[1,0])**2+\
        (self.h[0,1]-self.g[0,1])**2+(self.h[1,1]-self.g[1,1])**2)
      self.E[k+kE] += self.sigma*self.tmpx.sum()

      self.image_compose(self.I1, self.phiinvx, self.phiinvy, self.I)
      
      self.diffeo_gradient_y(self.phiinvy, self.yddx, self.yddy)
      self.diffeo_gradient_x(self.phiinvx, self.xddx, self.xddy)
      np.copyto(self.h[0,0], self.yddy*self.yddy+self.xddy*self.xddy)
      np.copyto(self.h[1,0], self.yddx*self.yddy+self.xddx*self.xddy)
      np.copyto(self.h[0,1], self.yddy*self.yddx+self.xddy*self.xddx)
      np.copyto(self.h[1,1], self.yddx*self.yddx+self.xddx*self.xddx)

      self.image_gradient(self.h[0,0], self.dhaadx, self.dhaady)
      self.image_gradient(self.h[0,1], self.dhabdx, self.dhabdy)
      self.image_gradient(self.h[1,0], self.dhbadx, self.dhbady)
      self.image_gradient(self.h[1,1], self.dhbbdx, self.dhbbdy)

      np.copyto(self.Jmap[0], -(self.h[0,0]-self.g[0,0])*self.dhaady -(self.h[0,1]-self.g[0,1])*self.dhabdy -(self.h[1,0]-self.g[1,0])*self.dhbady -(self.h[1,1]-self.g[1,1])*self.dhbbdy +\
        2*self.dhaady*self.h[0,0] + 2*self.dhabdx*self.h[0,0] + 2*self.dhbady*self.h[1,0] + 2*self.dhbbdx*self.h[1,0] +\
        2*(self.h[0,0]-self.g[0,0])*self.dhaady + 2*(self.h[1,0]-self.g[1,0])*self.dhbady + 2*(self.h[0,1]-self.g[0,1])*self.dhaadx + 2*(self.h[1,1]-self.g[1,1])*self.dhbadx)
      
      np.copyto(self.Jmap[1], -(self.h[0,0]-self.g[0,0])*self.dhaadx -(self.h[0,1]-self.g[0,1])*self.dhabdx -(self.h[1,0]-self.g[1,0])*self.dhbadx -(self.h[1,1]-self.g[1,1])*self.dhbbdx +\
        2*self.dhaady*self.h[0,1] + 2*self.dhabdx*self.h[0,1] + 2*self.dhbady*self.h[1,1] + 2*self.dhbbdx*self.h[1,1] +\
        2*(self.h[0,0]-self.g[0,0])*self.dhabdy + 2*(self.h[1,0]-self.g[1,0])*self.dhbbdy + 2*(self.h[0,1]-self.g[0,1])*self.dhabdx + 2*(self.h[1,1]-self.g[1,1])*self.dhbbdx)

      self.image_gradient(self.I, self.dIdx, self.dIdy)
      self.vx = -(self.I-self.I0)*self.dIdx + 2*self.sigma*self.Jmap[1]# axis: [1]
      self.vy = -(self.I-self.I0)*self.dIdy + 2*self.sigma*self.Jmap[0]# axis: [0]
      fftx = np.fft.fftn(self.vx)
      ffty = np.fft.fftn(self.vy)
      fftx *= self.Linv
      ffty *= self.Linv
      self.vx[:] = -np.fft.ifftn(fftx).real # vx[:]=smth will copy while vx=smth directs a pointer
      self.vy[:] = -np.fft.ifftn(ffty).real

      # STEP 4 (v = -grad E, so to compute the inverse we solve \psiinv' = -epsilon*v o \psiinv)
      np.copyto(self.tmpx, self.vx)
      self.tmpx *= epsilon
      np.copyto(self.psiinvx, self.idx)
      self.psiinvx -= self.tmpx
      if self.compute_phi: # Compute forward phi also (only for output purposes)
        np.copyto(self.psix, self.idx)
        self.psix += self.tmpx

      np.copyto(self.tmpy, self.vy)
      self.tmpy *= epsilon
      np.copyto(self.psiinvy, self.idy)
      self.psiinvy -= self.tmpy
      if self.compute_phi: # Compute forward phi also (only for output purposes)
        np.copyto(self.psiy, self.idy)
        self.psiy += self.tmpy

      self.diffeo_compose(self.phiinvx, self.phiinvy, self.psiinvx, self.psiinvy, \
                self.tmpx, self.tmpy) # Compute composition phi o psi = phi o (1-eps*v)
      np.copyto(self.phiinvx, self.tmpx)
      np.copyto(self.phiinvy, self.tmpy)
      if self.compute_phi: # Compute forward phi also (only for output purposes)
        self.diffeo_compose(self.psix, self.psiy, \
                  self.phix, self.phiy, \
                  self.tmpx, self.tmpy)
        np.copyto(self.phix, self.tmpx)
        np.copyto(self.phiy, self.tmpy)

if __name__ == '__main__':
  pass


# ON COPYTO
# https://stackoverflow.com/questions/6431973/how-to-copy-data-from-a-numpy-array-to-another
*/