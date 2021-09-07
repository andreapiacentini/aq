/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_LINEAROBSINSITU_H_
#define AQ_LINEAROBSINSITU_H_

#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "aq/LinearObsOpBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GeoVals;
  class ObsAuxControl;
  class ObsAuxIncrement;
  class ObsSpace;
  class ObsVec;

// -----------------------------------------------------------------------------
/// Streamfunction TL/AD observation operator for AQ model.

class LinearObsInSitu : public LinearObsOpBase,
                        private util::ObjectCounter<LinearObsInSitu> {
 public:
  static const std::string classname() {return "aq::LinearObsInSitu";}

  LinearObsInSitu(const ObsSpace &, const eckit::Configuration &);

// Obs Operators
  void setTrajectory(const GeoVals &, const ObsAuxControl &) override;
  void simulateObsTL(const GeoVals &, ObsVec &, const ObsAuxIncrement &) const override;
  void simulateObsAD(GeoVals &, const ObsVec &, ObsAuxIncrement &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_LINEAROBSINSITU_H_
