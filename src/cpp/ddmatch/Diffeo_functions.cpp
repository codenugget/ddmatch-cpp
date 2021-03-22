#include "Diffeo_functions.h"

// returns v0_idx, v1_idx, frac_dv
std::tuple<int, int, double> periodic_1d(const double v, const int s) {
  // NOTE: what should we do when v is much larger than int allows?
  assert(v <= std::numeric_limits<int>::max());
  assert(v >= std::numeric_limits<int>::lowest());
  int v0 = int(floor(v)); // floor(-2.1) = -3, floor(3.9) = 3
  int v1 = v0 + 1;
  double dv = v - double(v0);

  // If dv is negative it means that v is negative, so v0
  // is larger than v. We then reduce v0 and v1 by 1 and
  // after that impose the periodic boundary conditions.
  //if (dv < 0 or v0 < 0) // dv > 0
  if (v0 < 0)
  {
    //dv += 1; // dv > 0
    //v0 -= 1;
    v0 %= s;
    if (v0 < 0)
      v0 += s; // modulo works differently in c++ vs python for negative numbers
    //v1 -= 1;
    if (v1 < 0) {
      v1 %= s;
      if (v1 < 0)
        v1 += s; // modulo works differently in c++ vs python for negative numbers
    }
  }
  else if (v0 >= s) {
    v0 %= s;
    v1 %= s;
  }
  else if (v1 >= s) {
    v1 %= s;
  }
  return std::make_tuple(v0, v1, dv);
}

// returns v0_idx, v1_idx, v0_shift, v1_shift, frac_dv
std::tuple<int, int, double, double, double> periodic_1d_shift(const double v, const int s) {
  // NOTE: what should we do when v is much larger than int allows?
  assert(v <= std::numeric_limits<int>::max());
  assert(v >= std::numeric_limits<int>::lowest());
  int v0 = int(floor(v)); // floor(2.1) = -3, floor(3.9) = 3
  int v1 = v0 + 1;
  double dv = v - double(v0); // dv is always >= 0

  double v0_shift = 0.0;
  double v1_shift = 0.0;

  // If dv is negative it means that v is negative, so v0
  // is larger than v. We then reduce v0 and v1 by 1 and
  // after that impose the periodic boundary conditions.
  if (v0 < 0) {
    v0 %= s;
    if (v0 < 0)
      v0 += s; // modulo differs between c++ and python
    v0_shift = -double(s);

    if (v1 < 0) {
      v1 %= s;
      if (v1 < 0)
        v1 += s;
      v1_shift = -double(s);
    }
  }
  else if (v0 >= s) {
    v0 %= s;
    v1 %= s;
    v0_shift = double(s);
    v1_shift = double(s);
  }
  else if (v1 >= s) {
    v1 %= s;
    v1_shift = double(s);
  }

  return std::make_tuple(v0, v1, v0_shift, v1_shift, dv);
}

bool image_compose_2d(const dGrid& I, const dGrid& xphi, const dGrid& yphi, dGrid& Iout) {
  if (!I.is_same_shape(xphi) or !I.is_same_shape(yphi) or !I.is_same_shape(Iout))
    return false;
  int w = I.cols();
  int h = I.rows();
  if (w != h) // for now assume width == height
    return false;
  int s = w;
  for(int py = 0; py < h; ++py) {
    for(int px = 0; px < w; ++px) {
      const auto [x0, x1, dx] = periodic_1d(xphi[py][px], w);
      const auto [y0, y1, dy] = periodic_1d(yphi[py][px], h);

      double val = 0;
      val += I[y0][x0] * (1-dx) * (1-dy);
      val += I[y0][x1] * dx     * (1-dy);
      val += I[y1][x0] * (1-dx) * dy;
      val += I[y1][x1] * dx     * dy;
      Iout[py][px] = val;

      /*
      int x0_idx = int(xphi[py][px]);
      int y0_idx = int(yphi[py][px]);
      // QUESTION: why add to index here and not int x1_idx = int(xphi->get(px+1, py  , 0))?
      int x1_idx = x0_idx + 1;
      // QUESTION: why add to index here and not int y1_idx = int(xphi->get(px  , py+1, 0))?
      int y1_idx = y0_idx + 1;

      double frac_dx = xphi[py][px] - double(x0_idx);
      double frac_dy = yphi[py][px] - double(y0_idx);

      // Id frac_dx is negative it means that xphi is negative, so x0_idx
      // is larger than xphi. We then reduce x0_idx and x1_idx by 1 and
      // after that impose the periodic boundary conditions.
      if (frac_dx < 0 or x0_idx < 0) {
        frac_dx  += 1;
        x0_idx -= 1;
        x0_idx %= s; // QUESTION: replace s with 'w'?
        x1_idx -= 1;
        x1_idx %= s;
      } else if (x0_idx >= s) {
        x0_idx %= s;
        x1_idx %= s;
      } else if (x1_idx >= s) {
        x1_idx %= s;
      }

      if (frac_dy < 0 or y0_idx < 0) {
        frac_dy  += 1;
        y0_idx -= 1;
        y0_idx %= s; // QUESTION: replace s with 'h'?
        y1_idx -= 1;
        y1_idx %= s;
      } else if (y0_idx >= s) {
        y0_idx %= s;
        y1_idx %= s;
      } else if (y1_idx >= s) {
        y1_idx %= s;
      }

      double val = 0;
      val += I[y0_idx][x0_idx] * (1.-frac_dx) * (1.-frac_dy);
      val += I[y0_idx][x1_idx] * frac_dx      * (1.-frac_dy);
      val += I[y1_idx][x0_idx] * (1.-frac_dx) * frac_dy;
      val += I[y1_idx][x1_idx] * frac_dx      * frac_dy;
      Iout[py][px] = val;
      */
    }
  }
  return true;
/*
  def image_compose_2d(I,xphi,yphi,Iout):
    for i in range(s):
      for j in range(s):
        xind = int(xphi[i,j])
        yind = int(yphi[i,j])
        xindp1 = xind+1
        yindp1 = yind+1
        deltax = xphi[i,j]-float(xind)
        deltay = yphi[i,j]-float(yind)
        
        # Id xdelta is negative it means that xphi is negative, so xind
        # is larger than xphi. We then reduce xind and xindp1 by 1 and
        # after that impose the periodic boundary conditions.
        if (deltax < 0 or xind < 0):
          deltax += 1.0
          xind -= 1
          xind %= s
          xindp1 -= 1
          xindp1 %= s
        elif (xind >= s):
          xind %= s
          xindp1 %= s
        elif (xindp1 >= s):
          xindp1 %= s

        if (deltay < 0 or xind < 0):
          deltay += 1.0
          yind -= 1
          yind %= s
          yindp1 -= 1
          yindp1 %= s
        elif (yind >= s):
          yind %= s
          yindp1 %= s
        elif (yindp1 >= s):
          yindp1 %= s
          
        onemdeltax = 1.-deltax
        onemdeltay = 1.-deltay

        Iout[i,j] = I[yind,xind]     * onemdeltax * onemdeltay+\
                    I[yind,xindp1]   * deltax     * onemdeltay+\
                    I[yindp1,xind]   * deltay     * onemdeltax+\
                    I[yindp1,xindp1] * deltay     * deltax
*/
}

bool eval_diffeo_2d(
  const dGrid& xpsi,                const dGrid& ypsi,
  const std::vector<double>& xvect, const std::vector<double>& yvect,
        std::vector<double>& xout,        std::vector<double>& yout) {
  if (!xpsi.is_same_shape(ypsi))
    return false;
  int w = xpsi.cols();
  int h = xpsi.rows();
  if (w != h) // for now assume width == height
    return false;
  int s = w;
  // d = xvect.shape[0]
  int n = (int) xvect.size();
  for(int i = 0; i < n; ++i) {
    const auto [x0_idx, x1_idx, x0_shift, x1_shift, frac_dx] = periodic_1d_shift(xvect[i], w);
    const auto [y0_idx, y1_idx, y0_shift, y1_shift, frac_dy] = periodic_1d_shift(yvect[i], h);
    double val = 0;
    val += (xpsi[y0_idx][x0_idx] + x0_shift) * (1.-frac_dx)* (1.-frac_dy);
    val += (xpsi[y0_idx][x1_idx] + x1_shift) * frac_dx     * (1.-frac_dy);
    val += (xpsi[y1_idx][x0_idx] + x0_shift) * (1.-frac_dx)* frac_dy;
    val += (xpsi[y1_idx][x1_idx] + x1_shift) * frac_dx     * frac_dy;
    xout[i] = val;

    val = 0;
    val += (ypsi[y0_idx][x0_idx] + y0_shift) * (1.-frac_dx)* (1.-frac_dy);
    val += (ypsi[y0_idx][x1_idx] + y0_shift) * frac_dx     * (1.-frac_dy);
    val += (ypsi[y1_idx][x0_idx] + y1_shift) * (1.-frac_dx)* frac_dy;
    val += (ypsi[y1_idx][x1_idx] + y1_shift) * frac_dx     * frac_dy;
    yout[i] = val;

    /*
    int x0_idx = int(xvect[i]);
    int y0_idx = int(yvect[i]);
    int x1_idx = x0_idx + 1;
    int y1_idx = y0_idx + 1;
    double frac_dx = xvect[i] - double(x0_idx);
    double frac_dy = yvect[i] - double(y0_idx);
    double x0_shift = 0.0;
    double x1_shift = 0.0;
    double y0_shift = 0.0;
    double y1_shift = 0.0;

    // Id frac_dx is negative it means that xphi is negative, so x0_idx
    // is larger than xphi. We then reduce x0_idx and x1_idx by 1 and
    // after that impose the periodic boundary conditions.
    if (frac_dx < 0 or x0_idx < 0) {
      frac_dx  += 1;
      x0_idx -= 1;
      x0_idx %= s; // QUESTION: replace s with 'w'?
      x0_shift -= double(s);

      x1_idx -= 1;
      if (x1_idx < 0) {
        x1_idx %= s;
        x1_shift = -double(s);
      }
    } else if (x0_idx >= s) {
      x0_idx %= s;
      x1_idx %= s;
      x0_shift = double(s);
      x1_shift = double(s);
    } else if (x1_idx >= s) {
      x1_idx %= s;
      x1_shift = double(s);
    }

    if (frac_dy < 0 or y0_idx < 0) {
      frac_dy  += 1;
      y0_idx -= 1;
      y0_idx %= s; // QUESTION: replace s with 'h'?
      y0_shift -= double(s);

      y1_idx -= 1;
      if (y1_idx < 0) {
        y1_idx %= s;
        y1_shift = -double(s);
      }
    } else if (y0_idx >= s) {
      y0_idx %= s;
      y1_idx %= s;
      y0_shift = double(s);
      y1_shift = double(s);
    } else if (y1_idx >= s) {
      y1_idx %= s;
      y1_shift = double(s);
    }


    double val = 0;
    val += (xpsi[y0_idx][x0_idx] + x0_shift)* (1.-frac_dx)* (1.-frac_dy);
    val += (xpsi[y0_idx][x1_idx] + x1_shift)* frac_dx     * (1.-frac_dy);
    val += (xpsi[y1_idx][x0_idx] + x0_shift)* (1.-frac_dx)* frac_dy;
    val += (xpsi[y1_idx][x1_idx] + x1_shift)* frac_dx     * frac_dy;
    xout[i] = val;

    val = 0;
    val += (ypsi[y0_idx][x0_idx] + y0_shift)* (1.-frac_dx)* (1.-frac_dy);
    val += (ypsi[y0_idx][x1_idx] + y0_shift)* frac_dx     * (1.-frac_dy);
    val += (ypsi[y1_idx][x0_idx] + y1_shift)* (1.-frac_dx)* frac_dy;
    val += (ypsi[y1_idx][x1_idx] + y1_shift)* frac_dx     * frac_dy;
    yout[i] = val;
    */
  }
  return true;
/*
  # Evaluate diffeo psi(x,y) for each pair in xvect, yvect. 
  # Assuming psi is periodic.

    d = xvect.shape[0]

    for i in range(d):
      xind = int(xvect[i])
      yind = int(yvect[i])
      xindp1 = xind+1
      yindp1 = yind+1
      deltax = xvect[i]-float(xind)
      deltay = yvect[i]-float(yind)
      xshift = 0.0
      xshiftp1 = 0.0
      yshift = 0.0
      yshiftp1 = 0.0
      
      # Id xdelta is negative it means that xphi is negative, so xind
      # is larger than xphi. We then reduce xind and xindp1 by 1 and
      # after that impose the periodic boundary conditions.
      if (deltax < 0 or xind < 0):
        deltax += 1.0
        xind -= 1
        xindp1 -= 1
        xind %= s
        xshift = -float(s) # Should use floor_divide here instead.
        if (xindp1 < 0):
          xindp1 %= s
          xshiftp1 = -float(s)
      elif (xind >= s):
        xind %= s
        xindp1 %= s
        xshift = float(s)
        xshiftp1 = float(s)
      elif (xindp1 >= s):
        xindp1 %= s
        xshiftp1 = float(s)

      if (deltay < 0 or yind < 0):
        deltay += 1.0
        yind -= 1
        yindp1 -= 1
        yind %= s
        yshift = -float(s) # Should use floor_divide here instead.
        if (yindp1 < 0):
          yindp1 %= s
          yshiftp1 = -float(s)
      elif (yind >= s):
        yind %= s
        yindp1 %= s
        yshift = float(s)
        yshiftp1 = float(s)
      elif (yindp1 >= s):
        yindp1 %= s
        yshiftp1 = float(s)
        
      xout[i] = (xpsi[yind,xind]+xshift)*(1.-deltax)*(1.-deltay)+\
        (xpsi[yind,xindp1]+xshiftp1)*deltax*(1.-deltay)+\
        (xpsi[yindp1,xind]+xshift)*deltay*(1.-deltax)+\
        (xpsi[yindp1,xindp1]+xshiftp1)*deltay*deltax
      
      yout[i] = (ypsi[yind,xind]+yshift)*(1.-deltax)*(1.-deltay)+\
        (ypsi[yind,xindp1]+yshift)*deltax*(1.-deltay)+\
        (ypsi[yindp1,xind]+yshiftp1)*deltay*(1.-deltax)+\
        (ypsi[yindp1,xindp1]+yshiftp1)*deltay*deltax

  return eval_diffeo_2d
*/
}

//#define FIND_BUG
#ifdef FIND_BUG
#include <cstring>
#include <fstream>
void dump_data(const dGrid& g, const char* filename) {
  std::ofstream fp(filename, std::ios::binary);

  int w = g.cols();
  int h = g.rows();
  for(int y = 0; y < h; ++y) {
    for(int x = 0; x < w; ++x) {
      char buf[128];
      sprintf(buf, "%.4f, ", g[y][x]);
      auto l = strlen(buf);
      fp.write(buf, l);
    }
    fp.write("\n", 1);
  }
}
#endif

// NOTE: boundary_conditions1(-2.16e-16, 64) will return invalid numbers... add to gtest
// returns (i0,i1, i0shift, i1shift, fraction)
/*std::tuple<int,int, double,double,double> boundary_conditions1(const double value, const int sz) {
  int i0 = static_cast<int>(value);
  int i1 = i0 + 1;
  double i0_shift = 0.0;
  double i1_shift = 0.0;
  double frac = value - static_cast<double>(i0);
  // TODO: check that this holds for c++!
  // If frac_dx is negative it means that xphi is negative, so x0_idx
  // is larger than xphi. We then reduce x0_idx and x1_idx by 1 and
  // after that impose the periodic boundary conditions.
  if (frac < 0 or i0 < 0) {
    frac += 1.0;
    i0 -= 1;
    i1 -= 1;
    i0 = ((i0 % sz) + sz) % sz;
    i0_shift = -float(sz); // # Should use floor_divide here instead.
    if (i1 < 0) {
      i1 %= sz;
      i1_shift = -float(sz);
    }
  }
  else if (i0 >= sz) {
    i0 %= sz;
    i1 %= sz;
    i0_shift = float(sz);
    i1_shift = float(sz);
  }
  else if (i1 >= sz) {
    i1 %= sz;
    i1_shift = float(sz);
  }
  return { i0, i1, i0_shift, i1_shift, frac };
}*/

bool diffeo_compose_2d(
  const dGrid& xpsi, const dGrid& ypsi,
  const dGrid& xphi, const dGrid& yphi,
  dGrid& xout, dGrid& yout) {
  // Compute composition psi o phi. 
  // Assuming psi and phi are periodic.
  if (!xpsi.is_same_shape(ypsi) or !xpsi.is_same_shape(xphi) or !xpsi.is_same_shape(yphi))
    return false;

  int w = xpsi.cols();
  int h = xpsi.rows();
  if (w != h) // only allow square sizes now
    return false;

#ifdef FIND_BUG
  dump_data(xpsi, "xpsi.dat");
  dump_data(ypsi, "ypsi.dat");
  dump_data(xphi, "xphi.dat");
  dump_data(yphi, "yphi.dat");
  dump_data(xout, "xout.dat");
  dump_data(yout, "yout.dat");
#endif

#if 1
  for(int py = 0; py < h; ++py) {
    for(int px = 0; px < w; ++px) {
      if (px == 21) {
        int dbg = 0;
      }
      // NOTE: xphi seems to grow uncontrolled around edges (when writing xphi[0][62] is very large and causes an assertion failure, however neigboring values seems to be suspiciously large as well)
      const auto [x0_idx, x1_idx, x0_shift, x1_shift, frac_dx] = periodic_1d_shift(xphi[py][px], w);
      const auto [y0_idx, y1_idx, y0_shift, y1_shift, frac_dy] = periodic_1d_shift(yphi[py][px], h);
      double val = 0;
      val += (xpsi[y0_idx][x0_idx] + x0_shift) * (1.-frac_dx) * (1.-frac_dy);
      val += (xpsi[y0_idx][x1_idx] + x1_shift) * frac_dx      * (1.-frac_dy);
      val += (xpsi[y1_idx][x0_idx] + x0_shift) * (1.-frac_dx) * frac_dy;
      val += (xpsi[y1_idx][x1_idx] + x1_shift) * frac_dx      * frac_dy;
      xout[py][px] = val;

      val = 0;
      val += (ypsi[y0_idx][x0_idx] + y0_shift) * (1.-frac_dx) * (1.-frac_dy);
      val += (ypsi[y0_idx][x1_idx] + y0_shift) * frac_dx      * (1.-frac_dy);
      val += (ypsi[y1_idx][x0_idx] + y1_shift) * (1.-frac_dx) * frac_dy;
      val += (ypsi[y1_idx][x1_idx] + y1_shift) * frac_dx      * frac_dy;
      yout[py][px] = val;
    }
  }
  return true;
#else
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j) {
        int xind = int(xphi[i][j]);
        int yind = int(yphi[i][j]);
        int xindp1 = xind+1;
        int yindp1 = yind+1;
        double deltax = xphi[i][j]-float(xind);
        double deltay = yphi[i][j]-float(yind);
        double xshift = 0.0;
        double xshiftp1 = 0.0;
        double yshift = 0.0;
        double yshiftp1 = 0.0;
        // If xdelta is negative it means that xphi is negative, so xind
        // is larger than xphi. We then reduce xind and xindp1 by 1 and
        // after that impose the periodic boundary conditions.
        if (deltax < 0 or xind < 0) {
          deltax += 1.0;
          xind -= 1;
          xindp1 -= 1;
          xind %= s;
          xshift = -float(s);// # Should use floor_divide here instead.
          if (xindp1 < 0) {
            xindp1 %= s;
            xshiftp1 = -float(s);
          }
        }
        else if (xind >= s) {
          xind %= s;
          xindp1 %= s;
          xshift = float(s);
          xshiftp1 = float(s);
        }
        else if (xindp1 >= s) {
          xindp1 %= s;
          xshiftp1 = float(s);
        }

        if (deltay < 0 or yind < 0) {
          deltay += 1.0;
          yind -= 1;
          yindp1 -= 1;
          yind %= s;
          yshift = -float(s);// # Should use floor_divide here instead.
          if (yindp1 < 0) {
            yindp1 %= s;
            yshiftp1 = -float(s);
          }
        }
        else if (yind >= s) {
          yind %= s;
          yindp1 %= s;
          yshift = float(s);
          yshiftp1 = float(s);
        }
        else if (yindp1 >= s) {
          yindp1 %= s;
          yshiftp1 = float(s);
        }
        xout[i][j] =(xpsi[yind][xind]    +xshift  )*(1.-deltax)*(1.-deltay)+\
                    (xpsi[yind][xindp1]  +xshiftp1)*deltax     *(1.-deltay)+\
                    (xpsi[yindp1][xind]  +xshift  )*deltay     *(1.-deltax)+\
                    (xpsi[yindp1][xindp1]+xshiftp1)*deltay     *deltax;
        
        yout[i][j] =(ypsi[yind][xind]    +yshift  )*(1.-deltax)*(1.-deltay)+\
                    (ypsi[yind][xindp1]  +yshift  )*deltax     *(1.-deltay)+\
                    (ypsi[yindp1][xind]  +yshiftp1)*deltay     *(1.-deltax)+\
                    (ypsi[yindp1][xindp1]+yshiftp1)*deltay     *deltax;
    }
  }
  return true;
#endif
}
/*
  def diffeo_compose_2d(xpsi,ypsi,xphi,yphi,xout,yout):
  # Compute composition psi o phi. 
  # Assuming psi and phi are periodic.

    for i in range(s):
      for j in range(s):
        xind = int(xphi[i,j])
        yind = int(yphi[i,j])
        xindp1 = xind+1
        yindp1 = yind+1
        deltax = xphi[i,j]-float(xind)
        deltay = yphi[i,j]-float(yind)
        xshift = 0.0
        xshiftp1 = 0.0
        yshift = 0.0
        yshiftp1 = 0.0
        
        # Id xdelta is negative it means that xphi is negative, so xind
        # is larger than xphi. We then reduce xind and xindp1 by 1 and
        # after that impose the periodic boundary conditions.
        if (deltax < 0 or xind < 0):
          deltax += 1.0
          xind -= 1
          xindp1 -= 1
          xind %= s
          xshift = -float(s) # Should use floor_divide here instead.
          if (xindp1 < 0):
            xindp1 %= s
            xshiftp1 = -float(s)
        elif (xind >= s):
          xind %= s
          xindp1 %= s
          xshift = float(s)
          xshiftp1 = float(s)
        elif (xindp1 >= s):
          xindp1 %= s
          xshiftp1 = float(s)

        if (deltay < 0 or yind < 0):
          deltay += 1.0
          yind -= 1
          yindp1 -= 1
          yind %= s
          yshift = -float(s) # Should use floor_divide here instead.
          if (yindp1 < 0):
            yindp1 %= s
            yshiftp1 = -float(s)
        elif (yind >= s):
          yind %= s
          yindp1 %= s
          yshift = float(s)
          yshiftp1 = float(s)
        elif (yindp1 >= s):
          yindp1 %= s
          yshiftp1 = float(s)

        xout[i,j] = (xpsi[yind,xind]    +xshift  )*(1.-deltax)*(1.-deltay)+\
                    (xpsi[yind,xindp1]  +xshiftp1)*deltax     *(1.-deltay)+\
                    (xpsi[yindp1,xind]  +xshift  )*deltay     *(1.-deltax)+\
                    (xpsi[yindp1,xindp1]+xshiftp1)*deltay     *deltax
        
        yout[i,j] = (ypsi[yind,xind]    +yshift  )*(1.-deltax)*(1.-deltay)+\
                    (ypsi[yind,xindp1]  +yshift  )*deltax     *(1.-deltay)+\
                    (ypsi[yindp1,xind]  +yshiftp1)*deltay     *(1.-deltax)+\
                    (ypsi[yindp1,xindp1]+yshiftp1)*deltay     *deltax

  return diffeo_compose_2d
*/

bool diffeo_gradient_y_2d(const dGrid& I, dGrid& dIdx, dGrid& dIdy) {
  if (!I.is_same_shape(dIdx) or !I.is_same_shape(dIdy))
    return false;

  int w = I.cols();
  int h = I.rows();
  if (w != h)
    return false;

  int s = w;

  int i = 0;
  int j = 0;
  int ip1 = i+1;
  int jp1 = j+1;
  int im1 = s-1;
  int jm1 = s-1;
  int im2 = s-2;
  int jm2 = s-2;
  for (int j = 0; j < s; ++j) {
    dIdy[  0][j] = (I[ip1][j] - I[im1][j] + s)/2.0;
    dIdy[im1][j] = (I[  i][j] - I[im2][j] + s)/2.0;
  }
  for (int i = 0; i < s; ++i) {
    dIdx[i][  0] = (I[i][jp1] - I[i][jm1])/2.0;
    dIdx[i][jm1] = (I[i][  j] - I[i][jm2])/2.0;
  }
  im1 = 0;
  jm1 = 0;
  for(int i = 1; i < s-1; ++i) {
    ip1 = i+1;
    for(int j = 0; j < s; ++j) {
      dIdy[i][j] = (I[ip1][j] - I[im1][j])/2.0;
    }
    im1 = i;
  }

  for(int j = 1; j < s-1; ++j) {
    jp1 = j+1;
    for(int i = 0; i < s; ++i) {
      dIdx[i][j] = (I[i][jp1] - I[i][jm1])/2.0;
    }
    jm1 = j;
  }
  return true;
  /*
    def diffeo_gradient_y_2d(I,dIdx,dIdy):
    i = 0
    j = 0
    ip1 = i+1
    jp1 = j+1
    im1 = s-1
    jm1 = s-1
    im2 = s-2
    jm2 = s-2
    for j in range(s):
      dIdy[0,j] = (I[ip1,j]-I[im1,j]+s)/2.0
      dIdy[im1,j] = (I[i,j]-I[im2,j]+s)/2.0
    j = 0
    for i in range(s):
      dIdx[i,0] = (I[i,jp1]-I[i,jm1])/2.0
      dIdx[i,jm1] = (I[i,j]-I[i,jm2])/2.0
    im1 = 0
    jm1 = 0
    for i in range(1,s-1):
      ip1 = i+1
      for j in range(s):
        dIdy[i,j] = (I[ip1,j]-I[im1,j])/2.0
      im1 = i
    for j in range(1,s-1):
      jp1 = j+1
      for i in range(s):
        dIdx[i,j] = (I[i,jp1]-I[i,jm1])/2.0
      jm1 = j

  return diffeo_gradient_y_2d
*/
}

bool diffeo_gradient_x_2d(const dGrid& I, dGrid& dIdx, dGrid& dIdy) {
  if (!I.is_same_shape(dIdx) or !I.is_same_shape(dIdy))
    return false;

  int w = I.cols();
  int h = I.rows();
  if (w != h)
    return false;

  int s = w;

  int i = 0;
  int j = 0;
  int ip1 = i+1;
  int jp1 = j+1;
  int im1 = s-1;
  int jm1 = s-1;
  int im2 = s-2;
  int jm2 = s-2;
  for (int j = 0; j < s; ++j) {
    dIdy[  0][j] = (I[ip1][j] - I[im1][j])/2.0;
    dIdy[im1][j] = (I[  i][j] - I[im2][j])/2.0;
  }
  j = 0;
  for (int i = 0; i < s; ++i) {
    dIdx[i][  0] = (I[i][jp1] - I[i][jm1] + s)/2.0;
    dIdx[i][jm1] = (I[i][  j] - I[i][jm2] + s)/2.0;
  }
  im1 = 0;
  jm1 = 0;
  for(int i = 1; i < s-1; ++i) {
    ip1 = i+1;
    for(int j = 0; j < s; ++j) {
      dIdy[i][j] = (I[jp1][j] - I[im1][j])/2.0;
    }
    im1 = i;
  }
  for(int j = 1; j < s-1; ++j) {
    jp1 = j+1;
    for(int i = 0; i < s; ++i) {
      dIdx[i][j] = (I[i][jp1] - I[i][jm1])/2.0;
    }
    jm1 = j;
  }
  return true;
  /*
    def diffeo_gradient_x_2d(I,dIdx,dIdy):
    i = 0
    j = 0
    ip1 = i+1
    jp1 = j+1
    im1 = s-1
    jm1 = s-1
    im2 = s-2
    jm2 = s-2
    for j in range(s):
      dIdy[0,j] = (I[ip1,j]-I[im1,j])/2.0
      dIdy[im1,j] = (I[i,j]-I[im2,j])/2.0
    j = 0
    for i in range(s):
      dIdx[i,0] = (I[i,jp1]-I[i,jm1]+s)/2.0
      dIdx[i,jm1] = (I[i,j]-I[i,jm2]+s)/2.0
    im1 = 0
    jm1 = 0
    for i in range(1,s-1):
      ip1 = i+1
      for j in range(s):
        dIdy[i,j] = (I[ip1,j]-I[im1,j])/2.0
      im1 = i
    for j in range(1,s-1):
      jp1 = j+1
      for i in range(s):
        dIdx[i,j] = (I[i,jp1]-I[i,jm1])/2.0
      jm1 = j

  return diffeo_gradient_x_2d*/
}


bool image_gradient_2d(const dGrid& I, dGrid& dIdx, dGrid& dIdy) {
  if (!I.is_same_shape(dIdx) or !I.is_same_shape(dIdy))
    return false;

  int w = I.cols();
  int h = I.rows();
  if (w != h)
    return false;
  int s = w;

  int im1 = s-1;
  int jm1 = s-1;
  int jp1;
  for (int i = 0; i< s-1; ++i) {
    int ip1 = i+1;
    for (int j = 0; j < s-1; ++j) {
      jp1 = j+1;
      double val = (I[ip1][j] - I[im1][j])/2.0;
      dIdy[i][j] = val;
      val = (I[i][jp1] - I[jm1][i])/2.0;
      dIdx[i][j] = val;
      jm1 = j;
    }
    double val = (I[ip1][s-1] - I[im1][s-1])/2.0;
    dIdy[i][s-1] = val;
    val = (I[i][0] - I[i][s-2])/2.0;
    dIdx[i][s-1] = val;
    jm1 = s-1;
    im1 = i;
  }
  for(int j = 0; j < s-1; ++j) {
    jp1 = j+1;
    double val = (I[0][j] - I[im1][j])/2.0;
    dIdy[s-1][j] = val;
    val = (I[s-1][jp1] - I[s-1][jm1])/2.0;
    dIdx[s-1][j] = val;
    jm1 = j;
  }
  double val = (I[0][s-1] - I[s-2][s-1])/2.0;
  dIdy[s-1][s-1] = val;
  //dIdy[s-1,s-1] = (I[0,s-1]-I[s-2,s-1])/2.0
  val = (I[s-1][0] - I[s-1][s-2])/2.0;
  dIdx[s-1][s-1] = val;
  //dIdx[s-1,s-1] = (I[s-1,0]-I[s-1,s-2])/2.0
  return true;
/*
  def image_gradient_2d(I,dIdx,dIdy):
    im1 = s-1
    jm1 = s-1
    for i in range(s-1):
      ip1 = i+1
      for j in range(s-1):
        jp1 = j+1
        dIdy[i,j] = (I[ip1,j]-I[im1,j])/2.0
        dIdx[i,j] = (I[i,jp1]-I[i,jm1])/2.0
        jm1 = j
      dIdy[i,s-1] = (I[ip1,s-1]-I[im1,s-1])/2.0
      dIdx[i,s-1] = (I[i,0]-I[i,s-2])/2.0
      jm1 = s-1
      im1 = i
    for j in range(s-1):
      jp1 = j+1
      dIdy[s-1,j] = (I[0,j]-I[im1,j])/2.0
      dIdx[s-1,j] = (I[s-1,jp1]-I[s-1,jm1])/2.0
      jm1 = j
    dIdy[s-1,s-1] = (I[0,s-1]-I[s-2,s-1])/2.0
    dIdx[s-1,s-1] = (I[s-1,0]-I[s-1,s-2])/2.0
*/
}

bool image_gradient_2d_forward(const dGrid& I, dGrid& dIdx, dGrid& dIdy) {
  if (!I.is_same_shape(dIdx) or !I.is_same_shape(dIdy))
    return false;

  int w = I.cols();
  int h = I.rows();
  if (w != h)
    return false;
  int s = w;

  int im1 = s-1;
  int jm1 = s-1;
  int ip1;
  int jp1;
  for (int i = 0; i < s-1; ++i) {
    ip1 = i+1;
    for (int j = 0; j < s-1; ++j) {
      jp1 = j+1;
      double val = (I[ip1][j] - I[im1][j]) / 2.0;
      dIdy[i][j] = val;
      val = (I[i][jp1] - I[i][jm1]) / 2.0;
      dIdx[i][j] = val;
      jm1 = j;
    }
    double val = (I[ip1][s-1] - I[im1][s-1]) / 2.0;
    dIdy[i][s-1] = val;
    val = (I[i][0] - I[i][s-2]) / 2.0;
    dIdx[i][s-1] = val;
    jm1 = s-1;
    im1 = i;
  }
  for(int j = 0; j < s-1; ++j) {
    jp1 = j+1;
    double val = (I[0][j] - I[im1][j])/2.0;
    dIdy[s-1][j] = val;
    val = (I[s-1][jp1] - I[s-1][jm1])/2.0;
    dIdx[s-1][j] = val;
    jm1 = j;
  }
  double val = (I[0][s-1] - I[s-2][s-1])/2.0;
  dIdy[s-1][s-1] = val;
  val = (I[s-1][0] - I[s-1][s-2])/2.0;
  dIdx[s-1][s-1] = val;
  return true;
/*
  def image_gradient_2d_forward(I,dIdx,dIdy):
    im1 = s-1
    jm1 = s-1
    for i in range(s-1):
      ip1 = i+1
      for j in range(s-1):
        jp1 = j+1
        dIdy[i,j] = (I[ip1,j]-I[im1,j])/2.0
        dIdx[i,j] = (I[i,jp1]-I[i,jm1])/2.0
        jm1 = j
      dIdy[i,s-1] = (I[ip1,s-1]-I[im1,s-1])/2.0
      dIdx[i,s-1] = (I[i,0]-I[i,s-2])/2.0
      jm1 = s-1
      im1 = i
    for j in range(s-1):
      jp1 = j+1
      dIdy[s-1,j] = (I[0,j]-I[im1,j])/2.0
      dIdx[s-1,j] = (I[s-1,jp1]-I[s-1,jm1])/2.0
      jm1 = j
    dIdy[s-1,s-1] = (I[0,s-1]-I[s-2,s-1])/2.0
    dIdx[s-1,s-1] = (I[s-1,0]-I[s-1,s-2])/2.0

  return image_gradient_2d_forward
*/
}


bool divergence_2d(const dGrid& vx, const dGrid& vy, dGrid& divv) {
  if (!vx.is_same_shape(vy) or !vx.is_same_shape(divv))
    return false;

  int w = vx.cols();
  int h = vx.rows();
  if (w != h)
    return false;
  int s = w;

  int im1 = s-1;
  int jm1 = s-1;
  int ip1;
  int jp1;
  for (int i = 0; i < s-1; ++i) {
    ip1 = i+1;
    for (int j = 0; i < s-1; ++j) {
      jp1 = j+1;
      double dy = vy[ip1][  j] - vy[im1][  j];
      double dx = vx[  i][jp1] - vx[  i][jm1];
      divv[i][j] = (dx + dy)/2.0;
      //divv[i,j] = (vy[ip1,j]-vy[im1,j] + vx[i,jp1]-vx[i,jm1])/2.0
      jm1 = j;
    }
    double dy = vy[ip1][s-1] - vy[im1][s-1];
    double dx = vx[  i][  0] - vx[  i][s-2];
    divv[i][s-1] = (dx + dy)/2.0;
    //divv[i,s-1] = (vy[ip1,s-1]-vy[im1,s-1] + vx[i,0]-vx[i,s-2])/2.0
    jm1 = s-1;
    im1 = i;
  }
  for(int j = 0; j < s-1; ++j) {
    jp1 = j+1;
    //divv[s-1,j] = (vy[0,j]-vy[im1,j] + vx[s-1,jp1]-vx[s-1,jm1])/2.0
    double dy = vy[  0][  j] - vy[im1][  j];
    double dx = vx[s-1][jp1] - vx[s-1][jm1];
    divv[s-1][j] = (dx + dy)/2.0;
    jm1 = j;
  }
  //divv[s-1,s-1] = (vy[0,s-1]-vy[s-2,s-1] + vx[s-1,0]-vx[s-1,s-2])/2.0
  double dy = vy[s-1][  0] - vy[s-1][s-2];
  double dx = vx[  0][s-1] - vx[s-2][s-1];
  divv[s-1][s-1] = (dx + dy)/2.0;
  return true;
/*
  def divergence_2d(vx,vy,divv):
    im1 = s-1
    jm1 = s-1
    for i in range(s-1):
      ip1 = i+1
      for j in range(s-1):
        jp1 = j+1
        divv[i,j] = (vy[ip1,j]-vy[im1,j] + vx[i,jp1]-vx[i,jm1])/2.0
        jm1 = j
      divv[i,s-1] = (vy[ip1,s-1]-vy[im1,s-1] + vx[i,0]-vx[i,s-2])/2.0
      jm1 = s-1
      im1 = i
    for j in range(s-1):
      jp1 = j+1
      divv[s-1,j] = (vy[0,j]-vy[im1,j] + vx[s-1,jp1]-vx[s-1,jm1])/2.0
      jm1 = j
    divv[s-1,s-1] = (vy[0,s-1]-vy[s-2,s-1] + vx[s-1,0]-vx[s-1,s-2])/2.0

  return divergence_2d
*/
}

constexpr double det_2d(const double a11, const double a21, const double a12, const double a22) {
    return a11*a22 - a12*a21;
}


bool jacobian_2d_forward(const dGrid& xphi, const dGrid& yphi, dGrid& jac) {
  if (!xphi.is_same_shape(yphi) or !xphi.is_same_shape(jac))
    return false;

  int w = xphi.cols();
  int h = xphi.rows();
  if (w != h)
    return false;
  int s = w;

  for (int i = 0; i < s-1; ++i) {
    for (int j = 0; j < s-1; ++j) {
      double dxphi_dx = xphi[  i][j+1] - xphi[i][j];
      double dxphi_dy = xphi[i+1][  j] - xphi[i][j];
      double dyphi_dx = yphi[  i][j+1] - yphi[i][j];
      double dyphi_dy = yphi[i+1][  j] - yphi[i][j];
      double det = det_2d(dxphi_dx,dyphi_dx,
                          dxphi_dy,dyphi_dy);
      jac[i][j] = det;
    }
    // TODO: why +s here?
    double dxphi_dx = xphi[  i][  0] + s - xphi[i][s-1];
    double dxphi_dy = xphi[i+1][s-1]     - xphi[i][s-1];
    double dyphi_dx = yphi[  i][  0]     - yphi[i][s-1];
    double dyphi_dy = yphi[i+1][s-1]     - yphi[i][s-1];
    double det = det_2d(dxphi_dx,dyphi_dx,
                        dxphi_dy,dyphi_dy);

    jac[i][s-1] = det;
  }
  for (int j = 0; j < s-1; ++j) {
    double dxphi_dx = xphi[s-1][j+1] - xphi[s-1][j];
    double dxphi_dy = xphi[  0][  j] - xphi[s-1][j];
    double dyphi_dx = yphi[s-1][j+1] - yphi[s-1][j];
    double dyphi_dy = yphi[  0][  j] - yphi[s-1][j];
    double det = det_2d(dxphi_dx,dyphi_dx,
                        dxphi_dy,dyphi_dy);
    jac[s-1][j] = det;
  }

  // TODO: why +s here?
  double dxphi_dx = xphi[s-1][  0] + s - xphi[s-1][s-1];
  double dxphi_dy = xphi[  0][s-1]     - xphi[s-1][s-1];
  double dyphi_dx = yphi[s-1][  0]     - yphi[s-1][s-1];
  double dyphi_dy = yphi[  0][s-1] + s - yphi[s-1][s-1];
  double det = det_2d(dxphi_dx,dyphi_dx,
                      dxphi_dy,dyphi_dy);
  jac[s-1][s-1] = det;
  return true;
}
/*
  def jacobian_2d_forward(xphi,yphi,jac):
    for i in range(s-1):
      for j in range(s-1):
        dxphi_dx = xphi[i,j+1]-xphi[i,j]
        dxphi_dy = xphi[i+1,j]-xphi[i,j]
        dyphi_dx = yphi[i,j+1]-yphi[i,j]
        dyphi_dy = yphi[i+1,j]-yphi[i,j]
        jac[i,j] = det_2d(dxphi_dx,dyphi_dx,\
                  dxphi_dy,dyphi_dy)

      dxphi_dx = xphi[i,0]+s-xphi[i,s-1]
      dxphi_dy = xphi[i+1,s-1]-xphi[i,s-1]
      dyphi_dx = yphi[i,0]-yphi[i,s-1]
      dyphi_dy = yphi[i+1,s-1]-yphi[i,s-1]
      jac[i,s-1] = det_2d(dxphi_dx,dyphi_dx,\
                dxphi_dy,dyphi_dy)
    for j in range(s-1):
      dxphi_dx = xphi[s-1,j+1]-xphi[s-1,j]
      dxphi_dy = xphi[0,j]-xphi[s-1,j]
      dyphi_dx = yphi[s-1,j+1]-yphi[s-1,j]
      dyphi_dy = yphi[0,j]+s-yphi[s-1,j]
      jac[s-1,j] = det_2d(dxphi_dx,dyphi_dx,\
                dxphi_dy,dyphi_dy)

    dxphi_dx = xphi[s-1,0]+s-xphi[s-1,s-1]
    dxphi_dy = xphi[0,s-1]-xphi[s-1,s-1]
    dyphi_dx = yphi[s-1,0]-yphi[s-1,s-1]
    dyphi_dy = yphi[0,s-1]+s-yphi[s-1,s-1]
    jac[s-1,s-1] = det_2d(dxphi_dx,dyphi_dx,\
              dxphi_dy,dyphi_dy)

  return jacobian_2d_forward
*/