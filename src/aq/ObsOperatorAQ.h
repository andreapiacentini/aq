/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_OBSOPERATORAQ_H_
#define AQ_MODEL_OBSOPERATORAQ_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GomAQ;
  class LocationsAQ;
  class ObsBias;
  class ObsDiagsAQ;
  class ObsOpBaseAQ;
  class ObsSpaceAQ;
  class ObsVecAQ;

// -----------------------------------------------------------------------------

class ObsOperatorAQ : public util::Printable,
                      private boost::noncopyable {
 public:
  ObsOperatorAQ(const ObsSpaceAQ &, const eckit::Configuration &);
  ~ObsOperatorAQ();

/// Obs Operator
  void simulateObs(const GomAQ &, ObsVecAQ &, const ObsBias &, ObsDiagsAQ &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required input requiredVars from Model
  std::unique_ptr<LocationsAQ> locations() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOpBaseAQ> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_MODEL_OBSOPERATORAQ_H_
