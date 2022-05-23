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

#ifndef AQ_MODELAUXCOVARIANCE_H_
#define AQ_MODELAUXCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace aq {
  class ModelAuxControl;
  class ModelAuxIncrement;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelAuxCovariance : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<ModelAuxCovariance> {
 public:
  static const std::string classname() {return "aq::ModelAuxCovariance";}

/// Constructor, destructor
  ModelAuxCovariance(const eckit::Configuration & conf, const Geometry &): conf_(conf) {}
  ~ModelAuxCovariance() {}

/// Linear algebra operators
  void linearize(const ModelAuxControl &, const Geometry &) {}
  void multiply(const ModelAuxIncrement &, ModelAuxIncrement &) const {}
  void inverseMultiply(const ModelAuxIncrement &, ModelAuxIncrement &) const {}
  void randomize(ModelAuxIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODELAUXCOVARIANCE_H_
