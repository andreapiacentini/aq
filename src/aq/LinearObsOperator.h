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

#ifndef AQ_LINEAROBSOPERATOR_H_
#define AQ_LINEAROBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "aq/ObsOperatorParameters.h"

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace aq {
  class GeoVals;
  class LinearObsOpBase;
  class ObsAuxControl;
  class ObsAuxIncrement;
  class ObsSpace;
  class ObsVec;

// -----------------------------------------------------------------------------

class LinearObsOperator : public util::Printable,
                       private boost::noncopyable {
 public:
  typedef ObservationParameters Parameters_;

  LinearObsOperator(const ObsSpace &, const Parameters_ &);
  ~LinearObsOperator();

/// Obs Operator
  void setTrajectory(const GeoVals &, const ObsAuxControl &);
  void simulateObsTL(const GeoVals &, ObsVec &, const ObsAuxIncrement &) const;
  void simulateObsAD(GeoVals &, const ObsVec &, ObsAuxIncrement &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required inputs requiredVars from Model

 private:
  void print(std::ostream &) const;
  std::unique_ptr<LinearObsOpBase> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_LINEAROBSOPERATOR_H_
