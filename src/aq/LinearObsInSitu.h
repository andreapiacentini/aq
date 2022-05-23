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
