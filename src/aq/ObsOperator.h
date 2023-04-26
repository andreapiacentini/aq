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

#ifndef AQ_OBSOPERATOR_H_
#define AQ_OBSOPERATOR_H_

#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "aq/ObsOperatorParameters.h"

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace oops {
  template <typename OBS> class Locations;
}

namespace aq {
  class GeoVals;
  class ObsAuxControl;
  class ObsDiagnostics;
  class ObsOpBase;
  class ObsSpace;
  class ObsVec;
  struct ObsTraits;

// -----------------------------------------------------------------------------

class ObsOperator : public util::Printable,
                      private boost::noncopyable {
 public:
  typedef ObservationParameters Parameters_;
  typedef oops::Locations<ObsTraits> Locations_;

  ObsOperator(const ObsSpace &, const Parameters_ &);
  ~ObsOperator();

/// Obs Operator
  void simulateObs(const GeoVals &, ObsVec &, const ObsAuxControl &, ObsVec &,
                   ObsDiagnostics &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required input requiredVars from Model
  Locations_ locations() const;

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOpBase> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSOPERATOR_H_
