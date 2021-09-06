/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "aq/ObsAuxIncrement.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxCovariance.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace aq {
// -----------------------------------------------------------------------------
ObsAuxIncrement::ObsAuxIncrement(const ObsSpace &, const Parameters_ & params)
  : bias_(ObsAuxControl::ntypes, 0.0), active_(ObsAuxControl::ntypes, false)
{
  if (params.covariance.value() != boost::none) {
    const ObsAuxCovarianceParameters& covparams = *params.covariance.value();
    active_[0] = (covparams.stream.value() != boost::none);
    active_[1] = (covparams.uwind.value() != boost::none);
    active_[2] = (covparams.vwind.value() != boost::none);
    active_[3] = (covparams.wspeed.value() != boost::none);
  }
  bool on = false;
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (active_[jj]) {
      strn += "on";
      on = true;
    } else {
      strn += "off";
    }
  }
  if (on) {oops::Log::trace() << "ObsAuxIncrement created : " << strn << std::endl;}
}
// -----------------------------------------------------------------------------
ObsAuxIncrement::ObsAuxIncrement(const ObsAuxIncrement & other,
                                   const bool copy)
  : bias_(ObsAuxControl::ntypes, 0.0), active_(other.active_)
{
  if (copy) {
    for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] = other.bias_[jj];
  }
  this->makePassive();
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::makePassive() {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (!active_[jj]) bias_[jj] = 0.0;
  }
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::diff(const ObsAuxControl & b1, const ObsAuxControl & b2) {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    bias_[jj] = b1[jj] - b2[jj];
  }
  this->makePassive();
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::zero() {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] = 0.0;
}
// -----------------------------------------------------------------------------
ObsAuxIncrement & ObsAuxIncrement::operator=(const ObsAuxIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] = rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsAuxIncrement & ObsAuxIncrement::operator+=(const ObsAuxIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] += rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsAuxIncrement & ObsAuxIncrement::operator-=(const ObsAuxIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] -= rhs.bias_[jj];
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
ObsAuxIncrement & ObsAuxIncrement::operator*=(const double fact) {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] *= fact;
  this->makePassive();
  return *this;
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::axpy(const double fact, const ObsAuxIncrement & rhs) {
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) bias_[jj] += fact * rhs.bias_[jj];
  this->makePassive();
}
// -----------------------------------------------------------------------------
double ObsAuxIncrement::dot_product_with(const ObsAuxIncrement & rhs) const {
  double zz = 0.0;
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (active_[jj]) zz += bias_[jj] * rhs.bias_[jj];
  }
  return zz;
}
// -----------------------------------------------------------------------------
double ObsAuxIncrement::norm() const {
  double zz = 0.0;
  int ii = 0;
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (active_[jj]) {
      zz += bias_[jj] * bias_[jj];
      ++ii;
    }
  }
  if (ii > 0) zz = std::sqrt(zz/ii);
  return zz;
}
// -----------------------------------------------------------------------------
size_t ObsAuxIncrement::serialSize() const {
  size_t nn = bias_.size();
  return nn;
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::serialize(std::vector<double> & vect) const {
  vect.insert(vect.end(), bias_.begin(), bias_.end());
  oops::Log::trace() << "ObsAuxIncrement::serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::deserialize(const std::vector<double> & vect, size_t & index) {
  for (unsigned int jj = 0; jj < bias_.size(); ++jj) {
    bias_[jj] = vect[index];
    ++index;
  }
  oops::Log::trace() << "ObsAuxIncrement::deserialize done" << std::endl;
}
// -----------------------------------------------------------------------------
void ObsAuxIncrement::print(std::ostream & os) const {
  bool on = false;
  std::string strn = "";
  for (unsigned int jj = 0; jj < ObsAuxControl::ntypes; ++jj) {
    if (jj > 0) strn += ", ";
    if (active_[jj]) {
      on = true;
      std::ostringstream strs;
      strs << bias_[jj];
      strn += strs.str();
    } else {
      strn += "0.0";
    }
  }
  if (on) os << std::endl << "ObsAuxIncrement = " << strn;
}
// -----------------------------------------------------------------------------
}  // namespace aq
