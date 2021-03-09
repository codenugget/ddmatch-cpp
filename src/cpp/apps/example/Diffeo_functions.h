#pragma once

#include <memory>
#include <vector>

#include "core/MyArrays.h"

//bool image_compose_2d(const dGrid& I, const dGrid& xphi, const dGrid& yphi, dGrid& Iout);
bool image_compose_2d(dGrid& I, dGrid& xphi, dGrid& yphi, dGrid& Iout);

//bool eval_diffeo_2d(
//  const ImageLib::dImage* xpsi,     const ImageLib::dImage* ypsi,
//  const std::vector<double>& xvect, const std::vector<double>& yvect,
//        std::vector<double>& xout,        std::vector<double>& yout);
bool eval_diffeo_2d(
        dGrid& xpsi,                      dGrid& ypsi,
  const std::vector<double>& xvect, const std::vector<double>& yvect,
        std::vector<double>& xout,        std::vector<double>& yout);

//bool diffeo_compose_2d(
//  const ImageLib::dImage* xpsi, const ImageLib::dImage* ypsi,
//  const ImageLib::dImage* xphi, const ImageLib::dImage* yphi,
//        ImageLib::dImage* xout,       ImageLib::dImage* yout);
bool diffeo_compose_2d(
  dGrid& xpsi, dGrid& ypsi,
  dGrid& xphi, dGrid& yphi,
  dGrid& xout, dGrid& yout);


//bool diffeo_gradient_y_2d(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool diffeo_gradient_y_2d(dGrid& I, dGrid& dIdx, dGrid& dIdy);

//bool diffeo_gradient_x_2d(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool diffeo_gradient_x_2d(dGrid& I, dGrid& dIdx, dGrid& dIdy);

//bool image_gradient_2d(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool image_gradient_2d(dGrid& I, dGrid& dIdx, dGrid& dIdy);

//bool image_gradient_2d_forward(const ImageLib::dImage* I, ImageLib::dImage* dIdx, ImageLib::dImage* dIdy);
bool image_gradient_2d_forward(dGrid& I, dGrid& dIdx, dGrid& dIdy);

//bool divergence_2d(const ImageLib::dImage* vx, const ImageLib::dImage* vy, ImageLib::dImage* divv);
bool divergence_2d(dGrid&  vx, dGrid& vy, dGrid& divv);

bool jacobian_2d_forward(dGrid& xphi, dGrid& yphi, dGrid& jac);
//bool jacobian_2d_forward(const ImageLib::dImage* xphi, const ImageLib::dImage* yphi, ImageLib::dImage* jac);
