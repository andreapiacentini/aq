/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_OBSAUXCOVARIANCE_H_
#define AQ_OBSAUXCOVARIANCE_H_

#include <array>
#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/Printable.h"

#include "aq/ObsAuxControl.h"
#include "aq/ObsAuxParameters.h"

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
