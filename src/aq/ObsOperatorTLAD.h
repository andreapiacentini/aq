/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef AQ_OBSOPERATORTLAD_H_
#define AQ_OBSOPERATORTLAD_H_

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
  class GomAQ;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsOpBaseTLAD;
  class ObsSpaceAQ;
  class ObsVecAQ;

// -----------------------------------------------------------------------------

class ObsOperatorTLAD : public util::Printable,
                       private boost::noncopyable {
 public:
  typedef ObservationParameters Parameters_;

  ObsOperatorTLAD(const ObsSpaceAQ &, const Parameters_ &);
  ~ObsOperatorTLAD();

/// Obs Operator
  void setTrajectory(const GomAQ &, const ObsBias &);
  void simulateObsTL(const GomAQ &, ObsVecAQ &, const ObsBiasIncrement &) const;
  void simulateObsAD(GomAQ &, const ObsVecAQ &, ObsBiasIncrement &) const;

/// Other
  const oops::Variables & requiredVars() const;  // Required inputs requiredVars from Model

 private:
  void print(std::ostream &) const;
  std::unique_ptr<ObsOpBaseTLAD> oper_;
};

// -----------------------------------------------------------------------------

}  // namespace aq

#endif  // AQ_OBSOPERATORTLAD_H_
