/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
