/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_MODEL_OBSSTREAMTLAD_H_
#define AQ_MODEL_OBSSTREAMTLAD_H_

#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "oops/aq/ObsOpBaseTLAD.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GomAQ;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsSpaceAQ;
  class ObsVecAQ;

// -----------------------------------------------------------------------------
/// Streamfunction TL/AD observation operator for AQ model.

class ObsStreamTLAD : public ObsOpBaseTLAD,
                      private util::ObjectCounter<ObsStreamTLAD> {
 public:
  static const std::string classname() {return "aq::ObsStreamTLAD";}

  ObsStreamTLAD(const ObsSpaceAQ &, const eckit::Configuration &);

// Obs Operators
  void setTrajectory(const GomAQ &, const ObsBias &) override;
  void simulateObsTL(const GomAQ &, ObsVecAQ &, const ObsBiasIncrement &) const override;
  void simulateObsAD(GomAQ &, const ObsVecAQ &, ObsBiasIncrement &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_MODEL_OBSSTREAMTLAD_H_
