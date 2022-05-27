/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 CERFACS.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
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
