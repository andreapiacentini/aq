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

#ifndef AQ_COVARIANCE_H_
#define AQ_COVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/Geometry.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace aq {
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for AQ model.

class Covariance : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<Covariance> {
 public:
  static const std::string classname() {return "aq::Covariance";}

  Covariance(const Geometry &, const oops::Variables &,
                    const eckit::Configuration &, const State &, const State &);
  ~Covariance();

  void multiply(const Increment &, Increment &) const;
  void inverseMultiply(const Increment &, Increment &) const;
  void randomize(Increment &) const;

 private:
  void print(std::ostream &) const;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_COVARIANCE_H_
