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
