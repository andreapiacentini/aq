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

#ifndef AQ_OBSAUXCOVARIANCE_H_
#define AQ_OBSAUXCOVARIANCE_H_

#include <array>
#include <memory>
#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/Printable.h"

#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxParameters.h"
#include "aq/ObsAuxPreconditioner.h"

namespace aq {
  class ObsAuxControl;
  class ObsAuxIncrement;
  class ObsSpace;

// -----------------------------------------------------------------------------

class ObsAuxCovariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsAuxCovariance> {
 public:
  typedef ObsAuxControlParameters Parameters_;

  static const std::string classname() {return "aq::ObsAuxCovariance";}

/// Constructor, destructor
  ObsAuxCovariance(const ObsSpace &, const Parameters_ &);
  ~ObsAuxCovariance() {}

/// Linear algebra operators
  void linearize(const ObsAuxControl &, const eckit::Configuration &) {}
  void multiply(const ObsAuxIncrement &, ObsAuxIncrement &) const;
  void inverseMultiply(const ObsAuxIncrement &, ObsAuxIncrement &) const;
  void randomize(ObsAuxIncrement &) const;

  std::unique_ptr<ObsAuxPreconditioner> preconditioner() const;

  /// I/O and diagnostics
  void write(const Parameters_ &) const {}
  bool active(const unsigned int ii) const {return variance_[ii] > 0.0;}

 private:
  void print(std::ostream &) const;
  std::array<double, ObsAuxControl::ntypes> variance_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSAUXCOVARIANCE_H_
