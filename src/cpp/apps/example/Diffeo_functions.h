#pragma once

#include <memory>
#include <vector>

#include "image/Image.h"

bool image_compose_2d(const ImageLib::dImage* I, const ImageLib::dImage* xphi, const ImageLib::dImage* yphi, ImageLib::dImage* Iout);

bool eval_diffeo_2d(
  const ImageLib::dImage* xpsi,     const ImageLib::dImage* ypsi,
  const std::vector<double>& xvect, const std::vector<double>& yvect,
        std::vector<double>& xout,        std::vector<double>& yout);

bool diffeo_compose_2d(
  const ImageLib::dImage* xpsi, const ImageLib::dImage* ypsi,
  const ImageLib::dImage* xphi, const ImageLib::dImage* yphi,
        ImageLib::dImage* xout,       ImageLib::dImage* yout);


bool diffeo_gradient_y_2d(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool diffeo_gradient_x_2d(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool image_gradient_2d(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool image_gradient_2d_forward(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool divergence_2d(const ImageLib::dImage* vx, const ImageLib::dImage* vy, ImageLib::dImage* divv);
bool jacobian_2d_forward(const ImageLib::dImage* xphi, const ImageLib::dImage* yphi, ImageLib::dImage* jac);
