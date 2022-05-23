/*
 * This file is part of the Air Quality Ensemble Data Assimilation suite AQ.
 *
 * (C) Copyright 2022 CERFACS
 *
 * AQ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * AQ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * A copy of the GNU Lesser General Public License is distributed
 * along with AQ (files LICENSE.md, COPYING and COPYING.LESSER).
 */

#include "aq/ObsAuxControl.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "aq/ObsAuxIncrement.h"
#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
ObsAuxControl::ObsAuxControl(const ObsSpace &, const Parameters_ & params)
  : active_(false), geovars_(), hdiags_() {
  oops::Log::info() << "ObsAuxControl: conf = " << params << std::endl;
  bias_.fill(0.0);
  active_ = params.insitu.value() != boost::none;
  if (active_) {
    if (params.insitu.value() != boost::none) bias_[0] = *params.insitu.value();
    std::string strn = "";
    for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
      if (jj > 0) strn += ", ";
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    }
    oops::Log::info() << "ObsAuxControl::ObsAuxControl created, bias = " << strn << std::endl;
  }
}
// -----------------------------------------------------------------------------
ObsAuxControl::ObsAuxControl(const ObsAuxControl & other, const bool copy)
  : active_(other.active_),
    geovars_(other.geovars_), hdiags_(other.hdiags_)
{
  if (active_ && copy) {
    for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] = other.bias_[jj];
  } else {
    bias_.fill(0.0);
  }
}
// -----------------------------------------------------------------------------
ObsAuxControl & ObsAuxControl::operator+=(const ObsAuxIncrement & dx) {
  if (active_) {
    for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] += dx[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsAuxControl & ObsAuxControl::operator=(const ObsAuxControl & rhs) {
  if (active_) {
    for (unsigned int jj = 0; jj < ntypes; ++jj) bias_[jj] = rhs.bias_[jj];
  }
  return *this;
}
// -----------------------------------------------------------------------------
double ObsAuxControl::norm() const {
  double zz = 0.0;
  if (active_) {
    double ztmp = 0.0;
    std::size_t ii = 0;
    for (unsigned int jj = 0; jj < ntypes; ++jj) {
      ztmp = bias_[jj]*bias_[jj];
      zz += ztmp;
      if (ztmp > 0.0) ++ii;
    }
    zz = std::sqrt(zz/ii);
  }
  return zz;
}
// -----------------------------------------------------------------------------
void ObsAuxControl::print(std::ostream & os) const {
  if (active_) {
    std::string strn = "";
    for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
      if (jj > 0) strn += ", ";
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    }
    os << "ObsAuxControl = " << strn;
  }
}
// -----------------------------------------------------------------------------
}  // namespace aq
