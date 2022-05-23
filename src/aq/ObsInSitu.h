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

#ifndef AQ_OBSINSITU_H_
#define AQ_OBSINSITU_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "aq/ObsOpBase.h"
#include "aq/ObsSpace.h"
#include "aq/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace aq {
  class GeoVals;
  class ObsAuxControl;
  class ObsVec;

// -----------------------------------------------------------------------------
/// Streamfunction observation for AQ model.

class InSitu : public ObsOpBase,
                    private util::ObjectCounter<InSitu> {
 public:
  static const std::string classname() {return "aq::InSitu";}

  InSitu(const ObsSpace &, const eckit::Configuration &);

// Obs Operator
  void simulateObs(const GeoVals &, ObsVec &, const ObsAuxControl &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}
  std::unique_ptr<Locations> locations() const override;

 private:
  void print(std::ostream &) const override;
  const ObsSpace & obsdb_;
  const oops::Variables varin_;
};
// -----------------------------------------------------------------------------

}  // namespace aq
#endif  // AQ_OBSINSITU_H_
