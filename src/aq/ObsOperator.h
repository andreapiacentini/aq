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

#ifndef AQ_OBSOPERATOR_H_
#define AQ_OBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "aq/ObsOperatorParameters.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"


namespace aq {
  class GeoVals;
  class Locations;
  class ObsAuxControl;
  class ObsDiagnostics;
  class ObsOpBase;
  class ObsSpace;
  class ObsVec;

// -----------------------------------------------------------------------------

class ObsOperator : public util::Printable,
                      private boost::noncopyable {
 public:
  typedef ObservationParameters Parameters_;

  ObsOperator(const ObsSpace &, const Parameters_ &);
  ~ObsOperator();

/// Obs Operator
  void simulateObs(const GeoVals &, ObsVec &, const ObsAuxControl &, ObsVec &,
                   ObsDiagnostics &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required input requiredVars from Model
  std::unique_ptr<Locations> locations() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOpBase> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSOPERATOR_H_
