/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_ERRORCOVARIANCEAQ_H_
#define AQ_ERRORCOVARIANCEAQ_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/Configuration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "aq/AqFortran.h"
#include "aq/GeometryAQ.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace aq {
  class IncrementAQ;
  class StateAQ;

// -----------------------------------------------------------------------------
/// Background error covariance matrix for AQ model.

class ErrorCovarianceAQ : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ErrorCovarianceAQ> {
 public:
  static const std::string classname() {return "aq::ErrorCovarianceAQ";}

  ErrorCovarianceAQ(const GeometryAQ &, const oops::Variables &,
                    const eckit::Configuration &, const StateAQ &, const StateAQ &);
  ~ErrorCovarianceAQ();

  void multiply(const IncrementAQ &, IncrementAQ &) const;
  void inverseMultiply(const IncrementAQ &, IncrementAQ &) const;
  void randomize(IncrementAQ &) const;

 private:
  void print(std::ostream &) const;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_ERRORCOVARIANCEAQ_H_
