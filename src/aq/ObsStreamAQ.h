/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_OBSSTREAMAQ_H_
#define AQ_OBSSTREAMAQ_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "aq/AqTraits.h"
#include "aq/ObsOpBaseAQ.h"
#include "aq/ObsSpaceAQ.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GomAQ;
  class ObsBias;
  class ObsVecAQ;

// -----------------------------------------------------------------------------
/// Streamfunction observation for AQ model.

class ObsStreamAQ : public ObsOpBaseAQ,
                    private util::ObjectCounter<ObsStreamAQ> {
 public:
  static const std::string classname() {return "aq::ObsStreamAQ";}

  ObsStreamAQ(const ObsSpaceAQ &, const eckit::Configuration &);

// Obs Operator
  void simulateObs(const GomAQ &, ObsVecAQ &, const ObsBias &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}
  std::unique_ptr<LocationsAQ> locations() const override;

 private:
  void print(std::ostream &) const override;
  const ObsSpaceAQ & obsdb_;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_OBSSTREAMAQ_H_
