/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
