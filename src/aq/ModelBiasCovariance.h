/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODELBIASCOVARIANCE_H_
#define AQ_MODELBIASCOVARIANCE_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace aq {
  class ModelBias;
  class ModelBiasIncrement;
  class GeometryAQ;

// -----------------------------------------------------------------------------

class ModelBiasCovariance : public util::Printable,
                            private boost::noncopyable,
                            private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname() {return "aq::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration & conf, const GeometryAQ &): conf_(conf) {}
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const GeometryAQ &) {}
  void multiply(const ModelBiasIncrement &, ModelBiasIncrement &) const {}
  void inverseMultiply(const ModelBiasIncrement &, ModelBiasIncrement &) const {}
  void randomize(ModelBiasIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODELBIASCOVARIANCE_H_
