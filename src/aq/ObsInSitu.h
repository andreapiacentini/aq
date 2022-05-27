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
