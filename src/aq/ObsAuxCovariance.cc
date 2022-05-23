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

#include "aq/ObsAuxCovariance.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>

#include "aq/ObsAuxIncrement.h"
#include "aq/ObsAuxPreconditioner.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
ObsAuxCovariance::ObsAuxCovariance(const ObsSpace &, const Parameters_ & params)
{
  std::array<double, ObsAuxControl::ntypes> zz;
  zz.fill(0.0);
  if (params.covariance.value() != boost::none) {
    const ObsAuxCovarianceParameters& covparams = *params.covariance.value();
    if (covparams.insitu.value() != boost::none) zz[0] = *covparams.insitu.value();
  }
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (std::abs(zz[jj]) > 1.0e-8) {
      variance_[jj] = zz[jj] * zz[jj];
      std::ostringstream strs;
      strs << variance_[jj];
      strn += strs.str();
    } else {
      variance_[jj] = 0.0;
      strn += "0.0";
    }
  }
  oops::Log::info() << "ObsAuxCovariance created, variances = " << strn << std::endl;
}
// -----------------------------------------------------------------------------
void ObsAuxCovariance::multiply(const ObsAuxIncrement & dxin,
                                 ObsAuxIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] * variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsAuxCovariance::inverseMultiply(const ObsAuxIncrement & dxin,
                                        ObsAuxIncrement & dxout) const {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dxout[jj] = dxin[jj] / variance_[jj];
    } else {
      dxout[jj] = 0.0;
    }
  }
}
// -----------------------------------------------------------------------------
void ObsAuxCovariance::randomize(ObsAuxIncrement & dx) const {
  static util::NormalDistribution<double> dist(ObsAuxControl::ntypes, 0.0, 1.0, 4);
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (variance_[jj] > 0.0) {
      dx[jj] = dist[jj] * std::sqrt(variance_[jj]);
    } else {
      dx[jj] = 0.0;
    }
  }
}

// -----------------------------------------------------------------------------
std::unique_ptr<ObsAuxPreconditioner> ObsAuxCovariance::preconditioner() const {
  return std::make_unique<ObsAuxPreconditioner> (variance_);
}

// -----------------------------------------------------------------------------
void ObsAuxCovariance::print(std::ostream & os) const {
  os << "ObsAuxCovariance::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace aq
