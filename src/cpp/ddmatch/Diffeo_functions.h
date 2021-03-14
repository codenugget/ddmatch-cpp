#pragma once

#include <memory>
#include <vector>

#include "core/MyArrays.h"

bool image_compose_2d(const dGrid& I, const dGrid& xphi, const dGrid& yphi, dGrid& Iout);

bool eval_diffeo_2d(
  const dGrid& xpsi,                const dGrid& ypsi,
  const std::vector<double>& xvect, const std::vector<double>& yvect,
        std::vector<double>& xout,        std::vector<double>& yout);

bool diffeo_compose_2d(
  const dGrid& xpsi, const dGrid& ypsi,
  const dGrid& xphi, const dGrid& yphi,
  dGrid& xout, dGrid& yout);


bool diffeo_gradient_y_2d(const dGrid& I, dGrid& dIdx, dGrid& dIdy);
bool diffeo_gradient_x_2d(const dGrid& I, dGrid& dIdx, dGrid& dIdy);

bool image_gradient_2d(const dGrid& I, dGrid& dIdx, dGrid& dIdy);

bool image_gradient_2d_forward(const dGrid& I, dGrid& dIdx, dGrid& dIdy);

bool divergence_2d(const dGrid&  vx, const dGrid& vy, dGrid& divv);

bool jacobian_2d_forward(const dGrid& xphi, const dGrid& yphi, dGrid& jac);

// returns (i0,i1, i0shift, i1shift, fraction)
std::tuple<int, int, double, double, double> boundary_conditions1(const double value, const int sz);
